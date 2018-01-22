=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut

#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $va = $reg->get_adaptor("human","variation","variation");
  $pa = $reg->get_adaptor("human","variation","population");

  $va = $db->get_VariationAdaptor();
  $pa = $db->get_PopulationAdaptor();

  # Get a Variation by its internal identifier
  $var = $va->fetch_by_dbID(145);

  # fetch a variation by its name
  $var = $va->fetch_by_name('rs100');

  # check if the variation is failed and if so, get the failed description
  if ($var->is_failed()) {
    $desc = $var->failed_description();
  }

  # fetch all variations from a population
  $pop = $pa->fetch_by_name('PACIFIC');
  @vars = {$va->fetch_all_by_Population($pop)};
  
  # Modify the include_failed_variations flag in DBAdaptor to also return variations that have been flagged as failed
  $va->db->include_failed_variations(1);
  
  # Fetch all variations from a population, including variations flagged as failed
  @vars = {$va->fetch_all_by_Population($pop)};
  
=head1 DESCRIPTION

This adaptor provides database connectivity for Variation objects.
Variations (SNPs, etc.) may be retrieved from the Ensembl variation database by
several means using this module.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor;

use Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref check_ref wrap_array);

use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Variation::Allele;
use Bio::EnsEMBL::Utils::Iterator;

use DBI qw(:sql_types);

our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor');

my $DEFAULT_ITERATOR_CACHE_SIZE = 10000;

sub store {
    my ($self, $var) = @_;
    
	my $dbh = $self->dbc->db_handle;
    
    # look up source_id
    if(!defined($var->{_source_id})) {
        my $sth = $dbh->prepare(q{
            SELECT source_id FROM source WHERE name = ?
        });
        $sth->execute($var->{source}->name);
        
        my $source_id;
		$sth->bind_columns(\$source_id);
		$sth->fetch();
		$sth->finish();
		$var->{_source_id} = $source_id;
    }
    if( defined $var->{evidence}){
	## store these by attrib id to allow different values in different species
	my $aa = $self->db->get_AttributeAdaptor;

	foreach my $ev_term( @{$var->{evidence}} ){

	    my $ev_class_id = $aa->attrib_id_for_type_value('evidence',$ev_term);
	    push @{$var->{evidence_attribs}},  $ev_class_id;
	}
    }
    
    throw("No source ID found for source name ", $var->source_name) unless defined($var->{_source_id});
    
    my $sth = $dbh->prepare(q{
        INSERT INTO variation (
            source_id,
            name,
            ancestral_allele,
            flipped,
            class_attrib_id,
            somatic,
            minor_allele,
            minor_allele_freq,
            minor_allele_count,
            clinical_significance,
            evidence_attribs
        ) VALUES (?,?,?,?,?,?,?,?,?,?,?)
    });
    
    $sth->execute(
        $var->{source} ? $var->{source}->dbID : $var->{_source_id},
        $var->name,
        $var->ancestral_allele,
        $var->{flipped},
        $var->{class_attrib_id} || ( $var->{class_SO_term} && $self->db->get_AttributeAdaptor->attrib_id_for_type_value('SO_term', $var->{class_SO_term})) || 18,
        $var->is_somatic,
        $var->minor_allele,
        $var->minor_allele_frequency,
        $var->minor_allele_count,
        $var->{clinical_significance} ? (join ",", @{$var->{clinical_significance}}) : undef,
        $var->{evidence_attribs} ? (join ",", @{$var->{evidence_attribs}}) : undef  
    );
    
    $sth->finish;
    
    # get dbID
	my $dbID = $dbh->last_insert_id(undef, undef, 'variation', 'variation_id');
    $var->{dbID}    = $dbID;
    $var->{adaptor} = $self;

    if(defined $var->{synonyms}){
        $self->store_synonyms($var);
    }
}


sub store_synonyms{

    my $self = shift;
    my $var  = shift;

    my $dbh = $self->dbc->db_handle;
    my $sth = $dbh->prepare(q{
        INSERT IGNORE INTO variation_synonym (
        variation_id,
        name,
        source_id
     ) VALUES (?,?,?)
    });

    foreach my $source_name (keys %{$var->{synonyms}} ){  ### update this to allow storage of id

      my $source_adaptor = $self->db->get_SourceAdaptor();
      my $source = $source_adaptor->fetch_by_name($source_name);

      throw("No source found for name $source_name") unless defined $source;

      foreach my $name (keys %{$var->{synonyms}->{$source_name}}){
        $sth->execute($var->dbID(), $name, $source->dbID()) ;
      }
    }
    return $var;
}

sub update {
    my ($self, $var) = @_;
    
	my $dbh = $self->dbc->db_handle;
    
    # look up source_id
    if(!defined($var->{_source_id})) {
        my $sth = $dbh->prepare(q{
            SELECT source_id FROM source WHERE name = ?
        });
        $sth->execute($var->{source}->name);
        
        my $source_id;
		$sth->bind_columns(\$source_id);
		$sth->fetch();
		$sth->finish();
		$var->{_source_id} = $source_id;
    }
    
    throw("No source ID found for source name ", $var->{source}->name) unless defined($var->{_source_id});

 if( defined $var->{evidence}){
	## store these by attrib id to allow differnt values in different species
	my $aa = $self->db->get_AttributeAdaptor;

	foreach my $ev_term( @{$var->{evidence}} ){

	    my $ev_class_id = $aa->attrib_id_for_type_value('evidence',$ev_term);
	    push @{$var->{evidence_attribs}},  $ev_class_id;
	}
    }
    
    my $sth = $dbh->prepare(q{
        UPDATE variation
           SET source_id = ?,
               name = ?,
               ancestral_allele = ?,
               flipped = ?,
               class_attrib_id = ?,
               somatic = ?,
               minor_allele = ?,
               minor_allele_freq = ?,
               minor_allele_count = ?,
               clinical_significance = ?,
               evidence_attribs = ?
         WHERE variation_id = ?
    });
    
    $sth->execute(
        $var->{source} ? $var->{source}->dbID : $var->{_source_id},
        $var->name,
        $var->ancestral_allele,
        $var->{flipped},
        $var->{class_attrib_id} || $var->adaptor->db->get_AttributeAdaptor->attrib_id_for_type_value('SO_term', $var->{class_SO_term}) || 18,
        $var->is_somatic,
        $var->minor_allele,
        $var->minor_allele_frequency,
        $var->minor_allele_count,
        join(",",@{$var->{clinical_significance}}) || undef,
        join(",",@{$var->{evidence_attribs}}) || undef,   ### HERE
        $var->dbID
    );
    
    $sth->finish;
}


sub store_attributes{

    my ($self, $var) = @_;

    my $dbh = $self->dbc->db_handle;

    my $sth = $dbh->prepare(q{
        INSERT IGNORE INTO variation_attrib
        ( variation_id, attrib_id, value)
        VALUES (?,?,?)
    });


    foreach my $attrib ( keys %{$var->{attribs}}){

        my $attrib_id = $var->adaptor->db->get_AttributeAdaptor->attrib_id_for_type_value('var_att', $attrib);
        throw("No ID found for attrib_type ", $attrib) unless defined($attrib_id);


        $sth->execute(
            $var->{dbID},
            $attrib_id,
            $var->{attribs}->{$attrib}
        );
    }
    $sth->finish;
}

=head2 fetch_all

  Description: Returns a listref of all germline variations
  Returntype : listref of Variations
  Status     : Stable

=cut

sub fetch_all {
    my $self = shift;
    my $constraint = 'v.somatic = 0';
    return $self->generic_fetch($constraint);
}

=head2 fetch_all_somatic

  Description: Returns a listref of all somatic variations
  Returntype : listref of Variations
  Status     : Stable

=cut

sub fetch_all_somatic {
    my $self = shift;
    my $constraint = 'v.somatic = 1';
    return $self->generic_fetch($constraint);
}

=head2 fetch_Iterator

  Arg [1]    : int $cache_size (optional)
  Example    : $var_iterator = $var_adaptor->fetch_Iterator;
  Description: returns an iterator over all germline variations in the database
  Returntype : Bio::EnsEMBL::Utils::Iterator
  Exceptions : none
  Caller     : general
  Status     : Experimental

=cut

sub fetch_Iterator {
    my ($self, $cache_size) = @_;
    return $self->_generic_fetch_Iterator($cache_size, 'v.somatic = 0');
}

=head2 fetch_Iterator_somatic

  Arg [1]    : int $cache_size (optional)
  Example    : $var_iterator = $var_adaptor->fetch_Iterator;
  Description: returns an iterator over all somatic variations in the database
  Returntype : Bio::EnsEMBL::Utils::Iterator
  Exceptions : none
  Caller     : general
  Status     : Experimental

=cut

sub fetch_Iterator_somatic {
    my ($self, $cache_size) = @_;
    return $self->_generic_fetch_Iterator($cache_size, 'v.somatic = 1');
}

sub _generic_fetch_Iterator {

    my ($self, $cache_size, $constraint) = @_;

    my $full_constraint = $constraint ?
        $constraint.' AND '.$self->db->_exclude_failed_variations_constraint :
        $self->db->_exclude_failed_variations_constraint;

    # prepare and execute a query to fetch all dbIDs

    my $sth = $self->prepare(qq{
        SELECT      v.variation_id 
        FROM        (variation v, source s)
        LEFT JOIN   failed_variation fv on v.variation_id = fv.variation_id
        WHERE       v.source_id = s.source_id
        AND         $full_constraint
    });

    $sth->execute;

    my $var_id;

    $sth->bind_columns(\$var_id);

    # we probably can't fit all of these into memory at once though, 
    # so create an iterator that fetches $cache_size dbIDs from the
    # statement handle at a time and then fetches these objects, 
    # storing them in a cache. We then return variation objects 
    # from this cache one by one, before filling it again if 
    # necessary

    $cache_size ||= $DEFAULT_ITERATOR_CACHE_SIZE;
    
    my @cache;

    my $items_to_fetch = 1;

    return Bio::EnsEMBL::Utils::Iterator->new(sub{

        if (@cache == 0 && $items_to_fetch) {
            
            # our cache is empty, and there are still items to fetch, so
            # fetch the next chunk of dbIDs and create objects from them

            my @dbIDs;

            my $item_count = 0;

            while( $sth->fetch ) {

                push @dbIDs, $var_id;
                
                if (++$item_count == $cache_size) {
                    # we have fetched a cache's worth of dbIDs, so flag that
                    # there are still items to fetch and last out of the loop
                    $items_to_fetch = 1;
                    last;
                }
                
                # if this is the last row, this flag will be 0 outside the loop
                $items_to_fetch = 0;
            }

            $sth->finish unless $items_to_fetch;
            
            @cache = @{ $self->fetch_all_by_dbID_list(\@dbIDs) } if @dbIDs;
        }

        return shift @cache;
    });
}

=head2 fetch_by_dbID

  Arg [1]    : int $dbID
  Example    : $var = $var_adaptor->fetch_by_dbID(5526);
  Description: Retrieves a Variation object via its internal identifier.
               If no such variation exists undef is returned.
  Returntype : Bio::EnsEMBL::Variation::Variation
  Exceptions : throw if dbID arg is not defined
  Caller     : general, IndividualAdaptor
  Status     : Stable

=cut

sub fetch_by_dbID {
    my $self = shift;
    my $dbID = shift;
    
    # Now, it should be fine to just use the superclass method
    my $result = $self->SUPER::fetch_by_dbID($dbID);
    
    return $result;
    
}

sub _columns {
    my $self = shift;
  
    my @cols = (
        "v.variation_id",
        "v.name AS v_name",
        "v.class_attrib_id AS v_class_attrib_id",
        "v.source_id AS v_source_id",
        "v.somatic AS v_somatic",
        "v.flipped AS v_flipped",
        "v.ancestral_allele AS v_ancestral_allele",
        "vs.name AS vs_name",
        "s2.name AS vs_source_name",
        "v.minor_allele",
        "v.minor_allele_freq",
        "v.minor_allele_count",
        "v.clinical_significance",
        "v.evidence_attribs"
    );
    
    
    return @cols;
}

sub _tables {
    my $self = shift;
    
    my @tables = (
        ['variation', 'v'],
        ['source', 's1'],
        ['variation_synonym', 'vs'],
        ['source', 's2']
    );
    
    # If we are constraining on population_id, add the allele table
    push(@tables,['allele', 'a']) if ($self->{'_constrain_population'});

    
    return @tables;
}
 
sub _left_join {
    my $self = shift;
    
    my @left_join = (
        ['variation_synonym', 'v.variation_id = vs.variation_id'],
        ['source s2', 'vs.source_id = s2.source_id']
    );
    
 
    return @left_join;
}

sub _default_where_clause {
    my $self = shift;
    
    my $constraint = qq{
        s1.source_id = v.source_id
    };
    
    # If we are constraining on population_id, we should have a constraint on the allele tables as well
    $constraint .= qq{ AND a.variation_id = v.variation_id } if ($self->{'_constrain_population'});

    # constrain to variants passing QC filters or those which are cited, unless include_failed_variants set
    $constraint .= qq{ AND v.display = 1 }  unless $self->db->include_failed_variations() ;
    
    return $constraint;
}
 
=head2 fetch_by_name

  Arg [1]    : string $name
  Arg [2]    : string $source (optional)
  Example    : $var = $var_adaptor->fetch_by_name('rs1453');
  Description: Retrieves a variation object via its name
  Returntype : Bio::EnsEMBL::Variation::Variation
  Exceptions : throw if name argument is not defined
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_name {
    my $self = shift;
    my $name = shift;
    my $source = shift;

    throw('name argument expected') if(!defined($name));
    
    # hack fix for phencode names with '+' in getting encoded/decoded by the webcode
    $name =~ s/\s/\+/ if $self->db->species =~ /homo_sapiens/i;
      
    # Add a constraint on the name
    my $constraint = qq{v.name = ?};

    # If a source is given, add a constraint on the source  
    $constraint .= qq{ AND s1.name = ?} if ( defined $source );
    
    # Bind the parameters
    $self->bind_param_generic_fetch($name,SQL_VARCHAR);
    $self->bind_param_generic_fetch($source,SQL_VARCHAR) if ( defined $source );
    
    # Get the results from generic fetch method
    my $result = $self->generic_fetch($constraint);
        
    # We need to check the synonym table in case the name was not found in the variation table
    unless (scalar(@{$result})) {
        
        # Call the fetch by synonym method
        $result = wrap_array($self->fetch_by_synonym($name,$source));
        
    }
    
    # If we still can't find any result and the name looks like an ssId, try fetching by subsnp_id instead
    unless (scalar(@{$result}) || $name !~ m/^ss\d+$/) {
        
        $result = wrap_array($self->fetch_by_subsnp_id($name));
        
    }
        
    # Return the result
    return undef unless (scalar(@{$result}));
    return $result->[0];
}

# alias for fetch_by_name
sub fetch_by_stable_id {
		my $self = shift;
		return $self->fetch_by_name(@_);
}

# Internal method for getting the internal dbIDs for a list of names. Will also query the variation_synonym and allele (for subsnp_ids) tables
sub _name_to_dbID {
    my $self = shift;
    my $name_list = shift;
    my $synonym = shift;
    my $subsnp = shift;
    
    $name_list = wrap_array($name_list);
    throw ("A list of names is required") unless (scalar(@{$name_list}));
    
    # Use a hash to store the name to dbID mapping
    my %dbIDs;
    
    # Determine which columns we are returning
    my $cols = ($synonym && $subsnp ? "CONCAT('ss',v.subsnp_id) AS name, v.variation_id" : "v.name, v.variation_id");
    
    # Determine which table we are querying
    my $table = ($synonym ? ($subsnp ? 'allele' : 'variation_synonym') : 'variation');
    
    # Statement to get the dbIDs from variation or variation_synonym
    my $stmt = qq{
        SELECT
            $cols
        FROM
            $table v
        WHERE 
    };
    my $sth;
    
    # Work on batches of $batch_size;
    my $batch_size = 200;
    # Make a local copy of the list to work on
    my $local_list = [@{$name_list}];
    while (scalar(@{$local_list})) {
        
        # Get the next batch and construct the constraint
        my @names = splice(@{$local_list},0,$batch_size);
        my $constraint = "('" . join("','",@names) . "')";
        $constraint = ($synonym && $subsnp ? "v.subsnp_id" : "v.name") . qq{ IN $constraint};
        
        # Prepare a statement
        $sth = $self->prepare($stmt . qq{ $constraint });
        $sth->execute();
        
        # Fetch the results and populate the hash
        my ($name,$dbID);
        $sth->bind_columns(\$name,\$dbID);
        while ($sth->fetch()) {
            $dbIDs{$name} = $dbID;
        }
        
    }
    
    # If we are querying variation and have unmapped names, also query variation_synonym and allele
    if (!$synonym || !$subsnp) {
        
        # Get the unmapped names
        my @unmapped = grep {!exists($dbIDs{$_})} @{$name_list};
        
        # If we are going to query for subsnp_ids, get the names that look like ssIds and strip the ss prefix
        if ($synonym) {
            @unmapped = grep {$_ =~ s/^ss(\d+)$/$1/} @unmapped;
        }
        
        # Get the dbID mapping from the synonym table and add them to the hash if there are unmapped
        if (scalar(@unmapped)) {
            my $names = $self->_name_to_dbID(\@unmapped,1,$synonym); 
            map {$dbIDs{$_} = $names->{$_}} keys(%{$names});
        }
    }
        
    # Return a hashref with the name -> dbID mapping
    return \%dbIDs;
}

sub fetch_by_synonym {
    my $self = shift;
    my $name = shift;
    my $source = shift;
    
    # Do a query to get the current variation for the synonym and call a fetch method on this variation
    my $constraint = qq{vs.name = ?};
    $constraint .= qq{ AND s.name = ?} if (defined($source));
    
    # This statement will only return 1 row which is consistent with the behaviour of fetch_by_name.
    # However, the synonym name is only guaranteed to be unique in combination with the source 
    my $stmt = qq{
        SELECT
            vs.variation_id
        FROM
            variation_synonym vs JOIN
            source s ON (
                s.source_id = vs.source_id
            )
        WHERE
            $constraint
        LIMIT 1
    };
    my $sth = $self->prepare($stmt);
    $sth->bind_param(1,$name,SQL_VARCHAR);
    $sth->bind_param(2,$source,SQL_VARCHAR) if (defined($source));
    $sth->execute();
    
    # Bind the results
    my $dbID;
    $sth->bind_columns(\$dbID);
    # Fetch the results
    $sth->fetch();
    
    # Return undef in case no data could be found
    return undef unless (defined($dbID));
    
    # Call the fetch_by_name method using the updated name for the synonym
    return $self->fetch_by_dbID($dbID);
}


=head2 fetch_by_subsnp_id

  Arg [1]    : string $subsnp_id
  Example    : $var = $var_adaptor->fetch_by_subsnp_id('ss123');
  Description: Retrieves a variation object via a component subsnp ID
  Returntype : Bio::EnsEMBL::Variation::Variation
  Exceptions : throw if name argument is not defined
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_subsnp_id {
    my $self = shift;
    my $name = shift;
  
    throw('name argument expected') if(!defined($name));
  
    # Strip away any ss prefix
    $name =~ s/^ss//gi;

    $self->{'_constrain_population'} = 1;
    
    # Add a constraint on the subsnp_id
    my $constraint = qq{a.subsnp_id = ?};
    
    # Bind the parameters
    $self->bind_param_generic_fetch($name,SQL_INTEGER);
    
    # Get the results from generic fetch method
    my $result = $self->generic_fetch($constraint);
    
    delete($self->{'_constrain_population'});
    
    # Return the result
    return undef unless (scalar(@{$result}));
    return $result->[0];
    
}



=head2 fetch_all_by_source

  Arg [1]    : string $source_name
  Arg [2]    : int $primary
  Example    : $var = $var_adaptor->fetch_all_by_source();
  Description: Retrieves all Variation objects associated with a source. By default ($primary=0)
               returns variations that have the source or variation_synonym that have the source.
               If primary set to 1, it returns only variations where the primary name is associated
               with the source
  Returntype : listref of Bio::EnsEMBL::Variation::Variation
  Exceptions : thrown if source_name not provided
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_source {
    my $self = shift;
    my $source_name = shift;
    my $primary = shift;

    throw('name argument expected') unless (defined($source_name));
    
    # By default, returns ALL variation and variation_synonyms where source = $name
    $primary ||= 0; 

    # Add the constraint on the source name. If primary is true, only the variation source will be queried, otherwise the variation_synonym source will be as well
    my $constraint = qq{s1.name = ?};
    $constraint = qq{($constraint OR s2.name = ?)} unless ($primary);
    
    # If necessary, add the constraint for filtering failed variations
    $constraint .= qq{ AND } . $self->db->_exclude_failed_variations_constraint();
  
    # Bind the source name parameter
    $self->bind_param_generic_fetch($source_name,SQL_VARCHAR);
    $self->bind_param_generic_fetch($source_name,SQL_VARCHAR);

    # Execute the superclass generic fetch
    return $self->generic_fetch($constraint);

}


=head2 fetch_all_by_source_type

  Arg [1]    : string $source_type
  Arg [2]    : int $primary
  Example    : $var = $var_adaptor->fetch_all_by_source_type('chip');
  Description: Retrieves all Variation objects associated with a type of source. By default ($primary=0)
               returns variations that have the source or variation_synonym that have the source.
               If primary set to 1, it returns only variations where the primary name is associated
               with the source
  Returntype : listref of Bio::EnsEMBL::Variation::Variation
  Exceptions : thrown if source_type not provided
  Caller     : general
  Status     : At risk

=cut

sub fetch_all_by_source_type {
    my $self = shift;
    my $source_type = shift;
    my $primary = shift;
  
    throw('name argument expected') unless (defined($source_type));
	
    # Get the source names that match the type
    my $stmt = qq{
        SELECT
            name
        FROM
            source
        WHERE
            type = ?
    };
    my $sth = $self->prepare($stmt);
    $sth->bind_param(1,$source_type,SQL_VARCHAR);
    $sth->execute();
    
    # Fetch the results from one source at a time. This should probably be implemented in a more efficient manner.
    my $name;
    $sth->bind_columns(\$name);
    # Store in a hash to avoid duplicates
    my %out;
    while ($sth->fetch()) {
        
        my $result = $self->fetch_all_by_source($name,$primary);
        next unless (scalar(@{$result}));
        
        # Add the fetched results to the total hash
        map {$out{$_->dbID()} = $_} @{$result};        
    }
    
    # Return the reference to the list of variations
    return [values(%out)];
}


=head2 fetch_all_by_dbID_list

  Arg [1]    : reference to list of ints $list
  Example    : @vars = @{$va->fetch_all_by_dbID_list([124, 56, 90])};
  Description: Retrieves a set of variations via their internal identifiers.
               This is faster than repeatedly calling fetch_by_dbID if there
               are a large number of variations to retrieve
  Returntype : reference to list of Bio::EnsEMBL::Variation::Variation objects
  Exceptions : throw on bad argument
  Caller     : general, IndividualGenotypeAdaptor, PopulationGenotypeAdaptor
  Status     : Stable

=cut

sub fetch_all_by_dbID_list {
    my $self = shift;
    my $list = shift;

    # Get an Iterator for the list
    my $iterator = $self->fetch_Iterator_by_dbID_list($list);
    
    # Get an arrayref representing the contents of the Iterator
    my $result = $iterator->to_arrayref();
    
    return $result;
}

=head2 fetch_Iterator_by_dbID_list

  Arg [1]    : reference to list of ints $list
  Example    : $variation_iterator = $va->fetch_Iterator_by_dbID_list([124, 56, 90]);
  Description: Retrieves an iterator over a set of variations via their internal identifiers.
  Returntype : Bio::EnsEMBL::Utils::Iterator
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Experimental

=cut

sub fetch_Iterator_by_dbID_list {
    my ($self, $dbid_list, $cache_size) = @_;
    
    unless ((defined $dbid_list) && (ref $dbid_list eq 'ARRAY')) {
        throw("list reference argument is required");
    }

    $cache_size ||= $DEFAULT_ITERATOR_CACHE_SIZE;

    # create an iterator that fetches variations in blocks of
    # $cache_size and returns them in turn

    my @object_cache;

    return Bio::EnsEMBL::Utils::Iterator->new(sub {

            if (@object_cache == 0 && @$dbid_list > 0 ) {
                my @dbids = splice @$dbid_list, 0, $cache_size;
                
                # Create a constraint on the dbIDs
                my $id_str = "(" . join(",",@dbids) . ")";
                my $constraint = qq{v.variation_id IN $id_str};
                
                @object_cache = @{ $self->generic_fetch($constraint) };
            }

            return shift @object_cache;
        }
    );
}

=head2 fetch_all_by_name_list

  Arg [1]    : reference to list of names $list
  Example    : @vars = @{$va->fetch_all_by_name_list(["rs3", "rs1333049"])};
  Description: Retrieves a set of variations via their names. This is faster
               than repeatedly calling fetch_by_name if there are a large number
			   of variations to retrieve
  Returntype : reference to list of Bio::EnsEMBL::Variation::Variation objects
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut


sub fetch_all_by_name_list {
    my $self = shift;
    my $list = shift;

    # Get a list of dbIDs for the names
    my $dbIDs = $self->_name_to_dbID($list);
    
    # Then fetch the variations by dbID list instead
    my @dbID_list = values(%{$dbIDs});
    return $self->fetch_all_by_dbID_list(\@dbID_list);
}

=head2 get_all_sources

  Args        : none
  Example     : $sources = $va->get_all_sources();
  Description : Retrieves from the database all sources in the Source table
  ReturnType  : array ref of string
  Exceptions  : none
  Caller      : web
  Status      : Stable

=cut

sub get_all_sources{
    my $self = shift;
    my @sources;

    my $source_name;
    my $sth = $self->prepare(qq{SELECT name from source
				});
    $sth->execute();
    $sth->bind_columns(\$source_name);
    
    while ($sth->fetch()){
			push @sources, $source_name
    }
    $sth->finish();

    return \@sources;
}


=head2 get_source_version

  Arg[1]      : string $name
  Example     : $version = $va->get_source_version('dbSNP');
  Description : Retrieves from the database the version for the source given as an argument
  ReturnType  : int
  Exceptions  : none
  Caller      : general
  Status      : Stable

=cut

sub get_source_version{
    my $self = shift;
    my $name = shift;
    my $version;
    my $sth = $self->prepare(qq{SELECT version from source where name = ?
				});

    $sth->bind_param(1,$name,SQL_VARCHAR);

    $sth->execute();
    $sth->bind_columns(\$version);
    $sth->fetch();
    $sth->finish();

    return $version;
}

=head2 get_flanking_sequence

    Arg[1]      : int $variationID
    Example     : $flankinq_sequence = $va->get_flanking_sequence('652');
    Description : Retrieves the appropriate flanking sequence (five,three) for the variation
                  position from the genomic reference sequence
    ReturnType  : reference to a list containing (three_flank,five_flank)
    Exceptions  : throw when not possible to obtain sequence
    Caller      : general, Variation
    Status      : Stable

=cut


sub get_flanking_sequence{

  my $self = shift;
  my $variationID = shift;
  my $length = shift ||325;


  my $flanking_sequence; #reference to an array for the three_prime and five_prime seqs
  my ($seq_region_id, $seq_region_strand, $seq_region_start,  $seq_region_end, $up_seq_region_end, $down_seq_region_start,
$match );

  ## variant may have multiple mappings, some perfect (alignment_quality = 1) others not (alignment_quality =0)
  my $sth = $self->prepare(qq{
                              SELECT seq_region_id, seq_region_start, seq_region_end, seq_region_strand, alignment_quality
                              FROM variation_feature
                              WHERE variation_id = ?
                              ORDER BY alignment_quality DESC
                             });

  $sth->bind_param(1,$variationID,SQL_INTEGER);
  $sth->execute(); #retrieve the flank from the variation database
  $sth->bind_columns(\($seq_region_id, $seq_region_start,  $seq_region_end, $seq_region_strand,$match ));
  $sth->fetch();
  $sth->finish();


  unless( $seq_region_id){
      warn( "*****[ERROR]: No seq_region_id for SNP with dbID: $variationID. ".
            "Cannot retrieve flanking region******\n" );
      return;
  }
  unless( $seq_region_start &&  $seq_region_end){
      warn( "*****[ERROR]: No mapping for SNP with dbID: $variationID. ".
               "Cannot retrieve flanking region******\n" );
      return;
  }

  if($seq_region_start >  $seq_region_end){
    ##insertion - re-order start and end coordinates
    ($up_seq_region_end,$down_seq_region_start)  =  ($seq_region_end, $seq_region_start );
  }
  else{
    ## deletion or substitution - increment start and end coordinates to avoid variant base(s)
    ($up_seq_region_end,$down_seq_region_start)  =  ($seq_region_start,  $seq_region_end);
    $up_seq_region_end--;
    $down_seq_region_start++;
  }
  my $up_seq_region_start  = $up_seq_region_end - $length;
  my $down_seq_region_end  = $down_seq_region_start + $length;

  my $down_seq = $self->_get_flank_from_core( $seq_region_id,
                                              $down_seq_region_start,
                                              $down_seq_region_end,
                                              $seq_region_strand );

  my $up_seq = $self->_get_flank_from_core( $seq_region_id,
                                            $up_seq_region_start,
                                            $up_seq_region_end,
                                            $seq_region_strand );


  push @{$flanking_sequence},$down_seq,$up_seq; #add to the array the 3 and 5 prime sequences

  return $flanking_sequence;

}


=head2 fetch_all_by_Population

  Arg [1]    : Bio::EnsEMBL::Variation::Population
  Arg [2]	 : $minimum_frequency (optional)
  Example    : $pop = $pop_adaptor->fetch_by_dbID(1345);
               @vars = @{$va_adaptor->fetch_all_by_Population($pop)};
  Description: Retrieves all variations which are stored for a specified
               population. If $minimum_frequency is supplied, only variations
			   with a minor allele frequency (MAF) greater than
			   $minimum_frequency will be returned.
  Returntype : listref of Bio::EnsEMBL::Variation::Variation
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_Population {
    my $self = shift;
    my $pop = shift;
    my $freq = shift;
    
    assert_ref($pop,'Bio::EnsEMBL::Variation::Population');

    if(!defined($pop->dbID())) {
        warning("Cannot retrieve genotypes for population without set dbID");
        return [];
    }
  
    # Constraint the query using the population_id for the population
    my $constraint = qq{a.population_id = ?};
    
    # adjust frequency if given a percentage
    if (defined($freq)) {
	
	   $freq /= 100 if $freq > 1;
	   $constraint .= qq{ AND (IF(a.frequency > 0.5, 1-a.frequency, a.frequency) > ?) };
    }
  
    # Add the constraint for failed variations
    $constraint .= " AND " . $self->db->_exclude_failed_variations_constraint();
  
    # Bind the parameters
    $self->bind_param_generic_fetch($pop->dbID(),SQL_INTEGER);
    $self->bind_param_generic_fetch($freq,SQL_DOUBLE) if (defined($freq));
    
    # Set the flag to indicate that we are constraining on population and should not left join to allele
    $self->{'_constrain_population'} = 1;
    
    # Execute the generic fetch
    my $result = $self->generic_fetch($constraint);
    
    # Unset the flag
    delete($self->{'_constrain_population'});
    
    return $result;
}

=head2 fetch_all_by_VariationSet

  Arg [1]    : Bio::EnsEMBL::Variation::VariationSet
  Example    : @vars = @{$va_adaptor->fetch_all_by_VariationSet($vs)};
  Description: Retrieves all variations which are present in a specified
               variation set and its subsets.
  Returntype : listref of Bio::EnsEMBL::Variation::Variation
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_VariationSet {
    my $self = shift;
    return $self->_generic_fetch_by_VariationSet(0, @_);
}

=head2 fetch_Iterator_by_VariationSet

  Arg [1]    : Bio::EnsEMBL::Variation::VariationSet
  Example    : $var_iterator = $va_adaptor->fetch_Iterator_by_VariationSet($vs);
  Description: Retrieves an iterator for all variations which are present in a specified
               variation set and its subsets.
  Returntype : Bio::EnsEMBL::Utils::Iterator object
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : Experimental

=cut

sub fetch_Iterator_by_VariationSet {
    my $self = shift;
    my $set = shift;
    my $cache_size = shift || $DEFAULT_ITERATOR_CACHE_SIZE;
    
    # First, get ids for all subsets,
    my @var_set_ids = ($set->dbID);
    map {push(@var_set_ids,$_->dbID())} @{$set->adaptor->fetch_all_by_super_VariationSet($set)};
    my $var_set_id = join(",",@var_set_ids);
    
    # Prepare a query for getting the span of variation_ids
    my $stmt = qq{
        FROM
            variation_set_variation vsv LEFT JOIN
            failed_variation fv ON (
    	       fv.variation_id = vsv.variation_id
            )
        WHERE
            vsv.variation_set_id IN ($var_set_id)
    };
      
    # Add the constraint for failed variations 
    my $constraint = " ";
    $constraint .=  " AND fv.variation_id is null "
       unless $self->db->include_failed_variations() ;
    
    my $sth = $self->prepare(qq{SELECT MIN(vsv.variation_id), MAX(vsv.variation_id) $stmt $constraint});
    $sth->execute();
    my ($min_variation_id,$max_variation_id);
    $sth->bind_columns(\$min_variation_id,\$max_variation_id);
    $sth->fetch();
    $max_variation_id ||= 0;
    $min_variation_id ||= 1;
    
    # Prepare a statement for getting the ids in a range
    $sth = $self->prepare(qq{SELECT vsv.variation_id $stmt AND vsv.variation_id BETWEEN ? AND ? $constraint});
    
    # Internally, we keep an Iterator that works on the dbID span we're at
    my $iterator;
        
    return Bio::EnsEMBL::Utils::Iterator->new(sub {

        # If the iterator is empty, get a new chunk of dbIDs, unless we've fetched all dbIDs 
        unless ((defined($iterator) && $iterator->has_next()) || $min_variation_id > $max_variation_id) {
            
            ## check there are ids in the range to return
            my $count_sth = $self->prepare(qq{SELECT count(vsv.variation_id) $stmt AND vsv.variation_id BETWEEN ? AND ? $constraint});

            my $done = 0;
            while( $min_variation_id < $max_variation_id && $done == 0  ){
                $count_sth->execute($min_variation_id, $min_variation_id+$cache_size);
                my $count = $count_sth->fetchall_arrayref();
                if ($count->[0]->[0] > 0){
                    $done =1;
                }
                else{
                    $min_variation_id += ($cache_size + 1);
                }
            }

            # Get the next chunk of dbIDs
            $sth->execute($min_variation_id,$min_variation_id+$cache_size);
            $min_variation_id += ($cache_size + 1);
            
            # Use a hash to keep track of the seen dbIDs
            my %seen;
            
            # Loop over the dbIDs and avoid duplicates
            my $dbID;
            my @dbIDs;
            $sth->bind_columns(\$dbID);
            while ($sth->fetch()) {
                push (@dbIDs,$dbID) unless ($seen{$dbID}++);
            }
    
            # Get a new Iterator based on the new dbID span
            $iterator = $self->fetch_Iterator_by_dbID_list(\@dbIDs);
            
        }
        
        return $iterator->next();
    });
}

sub _generic_fetch_by_VariationSet {
    my $self = shift;
    my $want_iterator = shift;
    my $set = shift;
    
    assert_ref($set,'Bio::EnsEMBL::Variation::VariationSet');

    if(!defined($set->dbID())) {
        warning("Cannot retrieve variations for variation set without a dbID");
        return [];
    }
  
    # Get the unique dbIDs for all variations in this set and all of its subsets
    my $dbid_list = $self->_fetch_all_dbIDs_by_VariationSet($set);
 
    my $num_vars = @$dbid_list;

    if ($num_vars > 100_000 && !$want_iterator) {
        warn "This set contains a large number ($num_vars) of variations, these may not fit".
            "into memory at once, considering using fetch_Iterator_by_VariationSet instead";
    }

    # Use the dbIDs to get all variations and return them
    return $want_iterator ? 
        $self->fetch_Iterator_by_dbID_list($dbid_list) : 
        $self->fetch_all_by_dbID_list($dbid_list);
}

sub _fetch_all_dbIDs_by_VariationSet {
  my $self = shift;
  my $set = shift;
  
  # First, get ids for all subsets,
  
  my @var_set_ids = ($set->dbID);
  
  foreach my $var_set (@{$set->adaptor->fetch_all_by_super_VariationSet($set)}) {
    push @var_set_ids, $var_set->dbID;
  }
  
  my $set_str = "(" . join(",",@var_set_ids) .")";

  # Add the constraint for failed variations
  my $constraint = " ";
  $constraint .= " AND fv.variation_id is null "
       unless $self->db->include_failed_variations() ;
  
  # Then get the dbIDs for all these sets
  my $stmt = qq{
    SELECT DISTINCT
      vsv.variation_id
    FROM
      variation_set_variation vsv LEFT JOIN
      failed_variation fv ON (
	fv.variation_id = vsv.variation_id
      )
    WHERE
      vsv.variation_set_id in $set_str
      $constraint
  };

  my $sth = $self->prepare($stmt);
  
  $sth->execute();
  
  my @result;
  my $dbID;
  
  $sth->bind_columns(\$dbID);
  
  while ($sth->fetch()) {
        push @result, $dbID;
  }

  return \@result;
}

=head2 fetch_all_by_publication

  Arg [1]    : Bio::EnsEMBL::Variation::Publication
  Example    : @vars = @{$va_adaptor->fetch_all_by_publication($pub)};
  Description: Retrieves all variations which are cited in a specific publication
  Returntype : listref of Bio::EnsEMBL::Variation::Variation
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : At Risk
=cut

sub fetch_all_by_publication{

    my $self    = shift;
    my $pub_obj = shift;

    my @var;
    my $variation_id;

    my $stmt = "SELECT  vc.variation_id 
                from variation_citation vc 
                left join variation v on vc.variation_id = v.variation_id 
                where vc.publication_id = ? ";

    # Add the constraint for failed variations
    $stmt .= " AND " .  $self->db->_exclude_failed_variations_constraint()
        unless $self->db->include_failed_variations();

    my $sth = $self->prepare($stmt);
    $sth->execute( $pub_obj->dbID() );
    $sth->bind_columns(\$variation_id);
    while ($sth->fetch()){
      push @var, $self->fetch_by_dbID($variation_id);
    }
    $sth->finish;

    return \@var;

}

# Internal method
# Fetches attributes
sub _fetch_attribs_by_dbID {
    my $self = shift;
    my $id = shift; 
  
    throw("Cannot fetch attributes without dbID") unless defined($id);
    
    my $attribs = {};
 
    my $sth = $self->dbc->prepare(qq{
      SELECT at.value, a.value
      FROM variation_attrib a, attrib at
      WHERE a.attrib_id = at.attrib_id
      AND a.variation_id = ?
    });

    $sth->bind_param(1,$id,SQL_INTEGER);
    $sth->execute();

    my ($key, $value);

    $sth->bind_columns(\$key, \$value);
    $attribs->{$key} .= $value ."," while $sth->fetch;
    $sth->finish;
    $attribs->{$key} =~ s/\,$// if ($key);
    return $attribs;
} 


=head2 get_all_failed_descriptions

  Arg[1]      : Bio::EnsEMBL::Variation::Variation $variation
	               The variation object to get the failed descriptions for
  Example     : 
                my $failed_descriptions = $adaptor->get_all_failed_descriptions($var);
                if (scalar(@{$failed_descriptions})) {
		          print "The variation '" . $var->name() . "' has been flagged as failed because '" . join("' and '",@{$failed_descriptions}) . "'\n";
                }
		
  Description : Gets the unique descriptions for the reasons why the supplied variation has failed.
  ReturnType  : reference to list of strings
  Exceptions  : none
  Caller      : general
  Status      : Stable

=cut

sub get_all_failed_descriptions {
    my $self = shift;
    my $variation = shift;
    
    # Call the internal get method without any constraints
    my $description = $self->_internal_get_failed_descriptions($variation) || [];
    
    return $description;
}


# API-internal method for getting failed descriptions for a variation
sub _internal_get_failed_descriptions {
    my $self = shift;
    my $variation = shift;
    my $constraint = shift;
    
    # Assert that the object passed is a Variation
    assert_ref($variation,'Bio::EnsEMBL::Variation::Variation');
    
    my $stmt = qq{
        SELECT DISTINCT
            fd.description
        FROM
            failed_variation fv JOIN
            failed_description fd ON (
                fd.failed_description_id = fv.failed_description_id
            )
        WHERE
            fv.variation_id = ?
    };
    $stmt .= qq{ AND $constraint } if (defined($constraint));
    
    my $sth = $self->prepare($stmt);
    $sth->execute($variation->dbID());
    
    return [map {$_->[0]} @{$sth->fetchall_arrayref([0])}];
}

sub _get_flank_from_core{
    my $self = shift;
    my $seq_region_id = shift;
    my $seq_region_start = shift;
    my $seq_region_end = shift;
    my $seq_region_strand = shift;

    my $flanking_sequence;
    if (defined $self->db()->dnadb()){
	my $slice_adaptor = $self->db()->dnadb()->get_SliceAdaptor();
	my $slice = $slice_adaptor->fetch_by_seq_region_id($seq_region_id);
	if (!$slice){
	    throw("Not possible to obtain slice for seq_region_id \"$seq_region_id\"\n");
	}
	my $flank = $slice->subseq($seq_region_start,$seq_region_end,$seq_region_strand);
	return $slice->subseq($seq_region_start,$seq_region_end,$seq_region_strand);
    }
    return '';
}


sub _objs_from_sth {
    my ($self, $sth) = @_;
  
    my %row;

    # Create the row hash using column names as keys
    $sth->bind_columns( \( @row{ @{$sth->{NAME_lc} } } ));

    while ($sth->fetch) {

        # we don't actually store the returned object because
        # the _obj_from_row method stores them in a temporary
        # hash _temp_objs in $self 

        $self->_obj_from_row(\%row);
    }

    # Get the created objects from the temporary hash
    my @objs = values %{ $self->{_temp_objs} };
    delete $self->{_temp_objs};
   
    # If the flag for attaching allele objects was set, trigger loading of them by calling the getter method on the variation object
    if ($self->{_load_alleles}) {
        map {$_->get_all_Alleles()} @objs;
    }
    
    # Return the created objects 
    return \@objs;
}


sub _obj_from_row {

    my ($self, $row) = @_;

    return undef unless $row->{variation_id};
    
    # If the variation for this variation_id hasn't already been created, do that
    my $obj = $self->{_temp_objs}{$row->{variation_id}}; 
    
    unless (defined($obj)) {

        ## convert attrib ids to values for evidence
        my @evidence;
        if (defined($row->{evidence_attribs})) {
            my $aa  = $self->db->get_AttributeAdaptor;
            my @attrib_ids = split(/,/,$row->{evidence_attribs});
            foreach my $attrib_id (@attrib_ids){
                my $evidence_value = $aa->attrib_value_for_id($attrib_id);
                push @evidence, $evidence_value ; 
            }
        } 
	      my @clin_sig;
        if (defined($row->{clinical_significance})) {
            @clin_sig = split(/,/,$row->{clinical_significance}); 
        } 
        
        # Create the variation object
        $obj = Bio::EnsEMBL::Variation::Variation->new(
            -dbID   => $row->{variation_id},
            -ADAPTOR => $self,
            -NAME   => $row->{v_name},
            -_SOURCE_ID => $row->{v_source_id},
            -IS_SOMATIC => $row->{v_somatic},
            -FLIPPED => $row->{v_flipped},
            -ANCESTRAL_ALLELE => $row->{v_ancestral_allele},
            -FLANK_FLAG => $row->{fs_flank_flag},
            -CLASS_SO_TERM => $self->AttributeAdaptor()->attrib_value_for_id($row->{v_class_attrib_id}),
            -CLINICAL_SIGNIFICANCE => \@clin_sig,
            -MINOR_ALLELE => $row->{minor_allele},
            -MINOR_ALLELE_FREQUENCY => $row->{minor_allele_freq},
            -MINOR_ALLELE_COUNT => $row->{minor_allele_count},
            -EVIDENCE => \@evidence,
        );
        
        $self->{_temp_objs}{$row->{variation_id}} = $obj;
    }
    
    # Add a synonym if available
    if (defined($row->{vs_source_name}) && defined($row->{vs_name})) {
        $obj->add_synonym($row->{vs_source_name},$row->{vs_name});
    }
    
}

=head2 load_alleles

  Arg[1]      : boolean $load
	               A flag determining whether to load alleles at the same time as Variations are created (1) or to lazy-load them (0 = default)
  Example     : # Tell the VariationAdaptor to load alleles when creating variations
                my $variation_adaptor->load_alleles(1);
                my $variation = $variation_adaptor->fetch_by_dbID(1);
		
  Description : Sets the behaviour when it comes to fetching alleles at the time of object creation or on demand. This setting will be in effect for
                all fetch methods for the lifespan of the adaptor or until reset.
  ReturnType  : none
  Exceptions  : none
  Caller      : general
  Status      : At Risk

=cut

sub load_alleles {
    my $self = shift;
    
    $self->{_load_alleles} = defined(shift);
    
}

1;
