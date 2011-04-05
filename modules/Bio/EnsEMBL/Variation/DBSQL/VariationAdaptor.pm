=head1 LICENSE

 Copyright (c) 1999-2011 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
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

  #Êcheck if the variation is failed and if so, get the failed description
  if ($var->is_failed()) {
    $desc = $var->failed_description();
  }

  # fetch all variations from a population
  $pop = $pa->fetch_by_name('PACIFIC');
  @vars = {$va->fetch_all_by_Population($pop)};
  
  #ÊModify the include_failed_variations flag in DBAdaptor to also return variations that have been flagged as failed
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
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);

use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Variation::Allele;
use Bio::EnsEMBL::Utils::Iterator;

use DBI qw(:sql_types);

our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor');

my $DEFAULT_ITERATOR_CACHE_SIZE = 1000;

=head2 fetch_all

  Description: Returns a listref of all germline variations
  Returntype : listref of Variations
  Status     : At risk

=cut

sub fetch_all {
    my $self = shift;
    my $constraint = 'v.somatic = 0';
    return $self->generic_fetch($constraint);
}

=head2 fetch_all_somatic

  Description: Returns a listref of all somatic variations
  Returntype : listref of Variations
  Status     : At risk

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

  throw('dbID argument expected') if(!defined($dbID));

  # NB flanking sequence table is joined here and in fetch_by_name
  # so we can set flank_flag where the flanking sequence differs
  # from the reference - it is just set to 0 when fetching variations
  # by other methods since otherwise the join takes too long  
  my $sth = $self->prepare
    (q{SELECT v.variation_id, v.name, v.validation_status, v.class_attrib_id, s1.name, s1.description, s1.url, s1.type, v.somatic, v.ancestral_allele,
              a.allele_id, a.subsnp_id, a.allele, a.frequency, a.count, a.sample_id, vs.moltype,
              vs.name, s2.name, (fs.up_seq is not null OR fs.down_seq is not null)
       FROM   (variation v, source s1)
	       LEFT JOIN allele a ON v.variation_id = a.variation_id
		   LEFT JOIN flanking_sequence fs ON v.variation_id = fs.variation_id
                LEFT JOIN variation_synonym vs on v.variation_id = vs.variation_id 
   	            LEFT JOIN source s2 on  vs.source_id = s2.source_id
		    LEFT JOIN failed_variation fv on v.variation_id = fv.variation_id
       WHERE    v.source_id = s1.source_id
       AND    v.variation_id = ?});
  $sth->bind_param(1,$dbID,SQL_INTEGER);
  $sth->execute();

  my $result = $self->_objs_from_sth($sth);
  $sth->finish();
  return undef if(!@$result);

  return $result->[0];
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
  my $extra_sql = "";

  if ( defined $source ) {
    $extra_sql = qq(AND    s1.name = ?);
  }

  # NB flanking sequence table is joined here and in fetch_by_dbID
  # so we can set flank_flag where the flanking sequence differs
  # from the reference - it is just set to 0 when fetching variations
  # by other methods since otherwise the join takes too long
  my $sth = $self->prepare
    (qq{SELECT v.variation_id, v.name, v.validation_status, v.class_attrib_id, s1.name, s1.description, s1.url, s1.type, v.somatic, v.ancestral_allele, 
               a.allele_id, a.subsnp_id, a.allele, a.frequency, a.count, a.sample_id, vs.moltype,
              vs.name, s2.name, (fs.up_seq is not null OR fs.down_seq is not null)
	  FROM   (variation v, source s1)
	     LEFT JOIN allele a on v.variation_id = a.variation_id
		 LEFT JOIN flanking_sequence fs on v.variation_id = fs.variation_id
	      LEFT JOIN variation_synonym vs on v.variation_id = vs.variation_id 
              LEFT JOIN source s2 on  vs.source_id = s2.source_id
	      LEFT JOIN failed_variation fv on v.variation_id = fv.variation_id
       WHERE    v.source_id = s1.source_id
       AND    v.name = ?
       $extra_sql  
       ORDER BY a.allele_id});

  $sth->bind_param(1,$name,SQL_VARCHAR);
  $sth->bind_param(2,$source,SQL_VARCHAR) if defined $source;
  $sth->execute();

  my $result = $self->_objs_from_sth($sth);
  $sth->finish();

  if(!@$result) {
    # try again if nothing found, but check synonym table instead
    $sth = $self->prepare
      (qq{SELECT v.variation_id, v.name, v.validation_status, v.class_attrib_id, s1.name, s1.description, s1.url, s1.type, v.somatic, v.ancestral_allele, 
                 a.allele_id, a.subsnp_id, a.allele, a.frequency, a.count, a.sample_id, vs1.moltype,
                vs2.name, s2.name, (fs.up_seq is not null OR fs.down_seq is not null)
         FROM variation v, source s1, source s2, allele a,
                variation_synonym vs1, variation_synonym vs2, flanking_sequence fs
         WHERE  v.variation_id = a.variation_id
         AND    v.variation_id = vs1.variation_id
         AND    v.variation_id = vs2.variation_id
         AND    v.source_id = s1.source_id
         AND    vs2.source_id = s2.source_id
		 AND    v.variation_id = fs.variation_id
         AND    vs1.name = ?
                $extra_sql
         ORDER BY a.allele_id});
    $sth->bind_param(1,$name,SQL_VARCHAR);
    $sth->bind_param(2,$source,SQL_VARCHAR) if defined $source;
    $sth->execute();
    $result = $self->_objs_from_sth($sth);

    $sth->finish();
  }
  
  if(!@$result && $name =~ /^ss\d+/i) {
	# try again if the ID looks like a subsnp ID
	push @{$result}, $self->fetch_by_subsnp_id($name);
  }
  
  return undef if(!@$result);
  
  return $result->[0];
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
  
  $name =~ s/^ss//gi;

  throw('name argument expected') if(!defined($name));

  # NB flanking sequence table is joined here and in fetch_by_dbID
  # so we can set flank_flag where the flanking sequence differs
  # from the reference - it is just set to 0 when fetching variations
  # by other methods since otherwise the join takes too long
  my $sth = $self->prepare
    (qq{SELECT v.variation_id, v.name, v.validation_status, v.class_attrib_id, s1.name, s1.description, s1.url, s1.type, v.somatic, v.ancestral_allele, 
               a.allele_id, a.subsnp_id, a.allele, a.frequency, a.count, a.sample_id, vs.moltype,
              vs.name, s2.name, (fs.up_seq is not null OR fs.down_seq is not null)
	  FROM   (variation v, source s1)
	     LEFT JOIN allele a on v.variation_id = a.variation_id
		 LEFT JOIN flanking_sequence fs on v.variation_id = fs.variation_id
	     LEFT JOIN variation_synonym vs on v.variation_id = vs.variation_id 
         LEFT JOIN source s2 on  vs.source_id = s2.source_id
	     LEFT JOIN failed_variation fv on v.variation_id = fv.variation_id
       WHERE    v.source_id = s1.source_id
       AND    a.subsnp_id = ?
       ORDER BY a.allele_id});

  $sth->bind_param(1,$name,SQL_VARCHAR);
  $sth->execute();

  my $result = $self->_objs_from_sth($sth);
  $sth->finish();
  
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

  my @out;

  $primary ||= 0; #by default, returns ALL variation and variation_synonyms where source = $name

  throw('name argument expected') if(!defined($source_name));

  #ÊAdd the constraint for failed variations
  my $constraint = " AND " . $self->db->_exclude_failed_variations_constraint();
  
  my $sth = $self->prepare
    (qq{SELECT v.variation_id, v.name, v.validation_status, v.class_attrib_id, s1.name, s1.description, s1.url, s1.type, v.somatic, v.ancestral_allele, 
               a.allele_id, a.subsnp_id, a.allele, a.frequency, a.count, a.sample_id, vs.moltype,
              vs.name, s2.name, 0
	  FROM   (variation v, source s1)
	     LEFT JOIN allele a on v.variation_id = a.variation_id 
	      LEFT JOIN variation_synonym vs on v.variation_id = vs.variation_id 
              LEFT JOIN source s2 on  vs.source_id = s2.source_id
	      LEFT JOIN failed_variation fv on v.variation_id = fv.variation_id
       WHERE    v.source_id = s1.source_id
       AND    s1.name = ?
       $constraint
   });

  $sth->bind_param(1,$source_name,SQL_VARCHAR);
  $sth->execute();

  my $result = $self->_objs_from_sth($sth);
  $sth->finish();
  
  push @out, @{$result};

  #we need to include variation_synonym as well, where the variation was merged with dbSNP
  if (!$primary){
      $sth = $self->prepare
	  (qq{SELECT v.variation_id, v.name, v.validation_status, v.class_attrib_id, s1.name, s1.description, s1.url, s1.type, v.somatic, v.ancestral_allele,
	      a.allele_id, a.subsnp_id, a.allele, a.frequency, a.count, a.sample_id, vs1.moltype,
	      vs1.name, s2.name, 0
		  FROM   (variation v, source s1, source s2,  variation_synonym vs1)
		  LEFT JOIN allele a ON v.variation_id = a.variation_id
		  LEFT JOIN failed_variation fv on v.variation_id = fv.variation_id
		  WHERE  v.variation_id = vs1.variation_id
		  AND    v.source_id = s1.source_id
		  AND    vs1.source_id = s2.source_id
		  AND    s2.name = ?
		  $constraint
		  ORDER BY v.variation_id
	      });
    $sth->bind_param(1,$source_name,SQL_VARCHAR);
    $sth->execute();
    $result = $self->_objs_from_sth($sth);
    $sth->finish();      
     #need to merge both lists, trying to avoid duplicates
      push @out, @{$result};
  }
  return \@out;
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
  
  my @out;
  $primary ||= 0; #by default, returns ALL variation and variation_synonyms where source_type = $source_type
  
  throw('name argument expected') if(!defined($source_type));
	
  my $sth = $self->prepare (qq{SELECT name FROM source WHERE type = ?});

  $sth->bind_param(1,$source_type,SQL_VARCHAR);
  $sth->execute();

  my $result = $self->_objs_from_sth($sth);
  $sth->finish();
  
	# Loop to get all the variations by source
  foreach my $source_name (@{$result}) {
		my $var_list = $self->fetch_all_by_source($source_name,$primary);
		push @out, @{$var_list};
	}
	return \@out;
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
  Status     : At Risk

=cut

sub fetch_all_by_dbID_list {
  my $self = shift;
  my $list = shift;

  if(!defined($list) || ref($list) ne 'ARRAY') {
    throw("list reference argument is required");
  }

  return [] if(!@$list);

  my @out;

  # mysql is faster and we ensure that we do not exceed the max query size by
  # splitting large queries into smaller queries of 200 ids
  my $max = 200;

  while(@$list) {
    my @ids = (@$list > $max) ? splice(@$list, 0, $max) : splice(@$list, 0);

    my $id_str = (@ids > 1)  ? " IN (".join(',',@ids).")"   :   ' = \''.$ids[0].'\'';

    my $sth = $self->prepare
      (qq{SELECT v.variation_id, v.name, v.validation_status, v.class_attrib_id, s1.name, s1.description, s1.url, s1.type, v.somatic, v.ancestral_allele,
                 a.allele_id, a.subsnp_id, a.allele, a.frequency, a.count, a.sample_id, vs.moltype,
                 vs.name, s2.name, 0
	     FROM   (variation v, source s1)
		  LEFT JOIN allele a on v.variation_id = a.variation_id
		     LEFT JOIN variation_synonym vs on v.variation_id = vs.variation_id 
		       LEFT JOIN source s2 on  vs.source_id = s2.source_id
		         LEFT JOIN failed_variation fv on v.variation_id = fv.variation_id
          WHERE    v.source_id = s1.source_id
          AND    v.variation_id $id_str});
    $sth->execute();

    my $result = $self->_objs_from_sth($sth);

    $sth->finish();

    push @out, @$result if(@$result);
  }

  return \@out;
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
                @object_cache = @{ $self->fetch_all_by_dbID_list(\@dbids) };
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
  Status     : At Risk

=cut


sub fetch_all_by_name_list {
  my $self = shift;
  my $list = shift;

  if(!defined($list) || ref($list) ne 'ARRAY') {
    throw("list reference argument is required");
  }

  return [] if(!@$list);

  my @out;

  # mysql is faster and we ensure that we do not exceed the max query size by
  # splitting large queries into smaller queries of 200 ids
  my $max = 200;

  while(@$list) {
    my @ids = (@$list > $max) ? splice(@$list, 0, $max) : splice(@$list, 0);

    my $id_str = (@ids > 1)  ? " IN (".join(',', map {"'$_'"} @ids).")"   :   ' = \''.$ids[0].'\'';
    my $sth = $self->prepare
      (qq{SELECT v.variation_id, v.name, v.validation_status, v.class_attrib_id, s1.name, s1.description, s1.url, s1.type, v.somatic, v.ancestral_allele,
                 a.allele_id, a.subsnp_id, a.allele, a.frequency, a.count, a.sample_id, vs.moltype,
                 vs.name, s2.name, 0
	     FROM   (variation v, source s1)
		  LEFT JOIN allele a on v.variation_id = a.variation_id
		     LEFT JOIN variation_synonym vs on v.variation_id = vs.variation_id 
		       LEFT JOIN source s2 on  vs.source_id = s2.source_id
		         LEFT JOIN failed_variation fv on v.variation_id = fv.variation_id
          WHERE    v.source_id = s1.source_id
          AND    v.name $id_str});
    $sth->execute();

    my $result = $self->_objs_from_sth($sth);

    $sth->finish();

    # Implementing the alias search functionality in fetch_by_name
    
    ## Finding which elements in $list yielded no results above
    my %missing_variations_temp=();
    @missing_variations_temp{@ids}=();
    foreach(@$result){
      delete $missing_variations_temp{$_->name};
  }
    my @missing_names=keys(%missing_variations_temp);

    ## A near-verbatim copy of the code in fetch_by_name, changing only for execution method
    if(@missing_names) {
      my $missing_id_str = (@missing_names > 1)  ? " IN (".join(',', map {"'$_'"} @missing_names).")"   :   ' = \''.$missing_names[0].'\'';
      $sth = $self->prepare
        (qq{SELECT v.variation_id, v.name, v.validation_status, v.class_attrib_id, s1.name, s1.description, s1.url, s1.type, v.somatic, v.ancestral_allele, 
                   a.allele_id, a.subsnp_id, a.allele, a.frequency, a.count, a.sample_id, vs1.moltype,
                  vs2.name, s2.name, 0
           FROM variation v, source s1, source s2, allele a,
                  variation_synonym vs1, variation_synonym vs2, flanking_sequence fs
           WHERE  v.variation_id = a.variation_id
           AND    v.variation_id = vs1.variation_id
           AND    v.variation_id = vs2.variation_id
           AND    v.source_id = s1.source_id
           AND    vs2.source_id = s2.source_id
	  AND    v.variation_id = fs.variation_id
	    AND    vs1.name $missing_id_str
	    ORDER BY a.allele_id
        });
      $sth->execute();
      my $missing_result = $self->_objs_from_sth($sth);
      push @out, @$missing_result if (@$missing_result);
      $sth->finish();
      }
    push @out, @$result if(@$result);
  }
  return \@out;
}


=head2 get_all_sources

  Args        : none
  Example     : $sources = $va->get_all_sources();
  Description : Retrieves from the database all sources in the Source table
  ReturnType  : array ref of string
  Exceptions  : none
  Caller      : web
  Status      : At Risk

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


=head2 get_default_source

  Args      : none
  Example     : $default_source = $va->get_default_source();
  Description : Retrieves from the database the default source used for display purposes
  ReturnType  : string
  Exceptions  : none
  Caller      : web
  Status      : At Risk

=cut

sub get_default_source{
    my $self = shift;

    my $source_name;
    my $sth = $self->prepare(qq{SELECT meta_value from meta where meta_key = ?
				});
    $sth->bind_param(1,'source.default_source',SQL_VARCHAR);
    $sth->execute();
    $sth->bind_columns(\$source_name);
    $sth->fetch();
    $sth->finish();

    return $source_name;
}


=head2 get_source_version

  Arg[1]      : string $name
  Example     : $version = $va->get_source_version('dbSNP');
  Description : Retrieves from the database the version for the source given as an argument
  ReturnType  : int
  Exceptions  : none
  Caller      : general
  Status      : At Risk

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
    Description : Retrieves from the database the appropriate flanking sequence (five,three) for the variation. If the flanking sequence is not in
                  the Flankinq_sequence table, access the core database with the coordinates
    ReturnType  : reference to a list containing (three_flank,five_flank)
    Exceptions  : throw when not possible to obtain sequence
    Caller      : general, Variation
    Status      : Stable

=cut

sub get_flanking_sequence{
  my $self = shift;
  my $variationID = shift;

  my $flanking_sequence; #reference to an array for the three_prime and five_prime seqs
  my ($seq_region_id, $seq_region_strand, $up_seq, $down_seq, $up_seq_region_start, $up_seq_region_end, $down_seq_region_start, $down_seq_region_end);
  
  my $sth = $self->prepare(qq{
			      SELECT seq_region_id, seq_region_strand, up_seq, 
			      down_seq, up_seq_region_start, up_seq_region_end, 
			      down_seq_region_start, down_seq_region_end
			      FROM flanking_sequence
			      WHERE variation_id = ?
			     });

  $sth->bind_param(1,$variationID,SQL_INTEGER);
  $sth->execute(); #retrieve the flank from the variation database
  $sth->bind_columns(\($seq_region_id, $seq_region_strand, $up_seq, $down_seq, $up_seq_region_start, $up_seq_region_end, $down_seq_region_start, $down_seq_region_end));
  $sth->fetch();
  $sth->finish();
 
  if (!defined $down_seq){
	if( $seq_region_id){
	  $down_seq = $self->_get_flank_from_core($seq_region_id, 
											  $down_seq_region_start, 
											  $down_seq_region_end, 
											  $seq_region_strand);
	} else {
	  warn( "*****[ERROR]: No seq_region_id for SNP with dbID: $variationID. ".
			"Cannot retrieve flanking region******\n" );    
	}
  }
  if (!defined $up_seq){
	if( $seq_region_id){
	  $up_seq = $self->_get_flank_from_core($seq_region_id, 
											$up_seq_region_start, 
											$up_seq_region_end, 
											$seq_region_strand);
	} else {
	  warn( "*****[ERROR]: No seq_region_id for SNP with dbID: $variationID. ".
			"Cannot retrieve flanking region******\n" );
	}
  }
  
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

  if(!ref($pop) || !$pop->isa('Bio::EnsEMBL::Variation::Population')) {
    throw('Bio::EnsEMBL::Variation::Population argument expected');
  }

  if(!defined($pop->dbID())) {
    warning("Cannot retrieve genotypes for population without set dbID");
    return [];
  }
  
  my $freq = shift;
  my $extra_sql = '';
  
  if(defined $freq) {
	
	# adjust frequency if given a percentage
	$freq /= 100 if $freq > 1;
	$extra_sql = qq{ AND (IF(a.frequency > 0.5, 1-a.frequency, a.frequency) > $freq) }
  }
  
  #ÊAdd the constraint for failed variations
  $extra_sql .= " AND " . $self->db->_exclude_failed_variations_constraint();
  
  my $sth = $self->prepare
    (qq{SELECT v.variation_id, v.name, v.validation_status, v.class_attrib_id, s1.name, s1.description, s1.url, s1.type, v.somatic, v.ancestral_allele,
               a.allele_id, a.subsnp_id, a.allele, a.frequency, a.count, a.sample_id, vs.moltype,
              vs.name, s2.name, 0
	    FROM   (variation v, source s1, allele a)
	      LEFT JOIN variation_synonym vs on v.variation_id = vs.variation_id 
              LEFT JOIN source s2 on  vs.source_id = s2.source_id
	      LEFT JOIN failed_variation fv on v.variation_id = fv.variation_id
	      WHERE  v.variation_id = a.variation_id
	      AND    v.source_id = s1.source_id
	      AND    a.sample_id = ?
		  $extra_sql});
  $sth->bind_param(1,$pop->dbID,SQL_INTEGER);
  $sth->execute();

  my $results = $self->_objs_from_sth($sth);
  $sth->finish();

  return $results;
}

=head2 fetch_all_by_VariationSet

  Arg [1]    : Bio::EnsEMBL::Variation::VariationSet
  Example    : @vars = @{$va_adaptor->fetch_all_by_VariationSet($vs)};
  Description: Retrieves all variations which are present in a specified
               variation set and its subsets.
  Returntype : listref of Bio::EnsEMBL::Variation::Variation
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : At Risk

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
    return $self->_generic_fetch_by_VariationSet(1, @_);
}

sub _generic_fetch_by_VariationSet {
    my $self = shift;
    my $want_iterator = shift;
    my $set = shift;
    
    if(!ref($set) || !$set->isa('Bio::EnsEMBL::Variation::VariationSet')) {
        throw('Bio::EnsEMBL::Variation::VariationSet argument expected');
    }

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
  
  my $set_str = '('.join(',',@var_set_ids).')';

  #ÊAdd the constraint for failed variations
  my $constraint = " AND " . $self->db->_exclude_failed_variations_constraint();
  
  # Then get the dbIDs for all these sets
  my $stmt = qq{
    SELECT
      DISTINCT(vsv.variation_id)
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

=head2 is_failed

  Description : DEPRECATED. The appropriate subroutine on the Variation/Allele object should be used instead.
  Exceptions  : Thrown on invocation.
  Status      : DEPRECATED

=cut

sub is_failed {
  my $self = shift;
  
    throw("The is_failed subroutine in VariationAdaptor is deprecated. Use the appropriate subroutine on the Variation/Allele object instead");
}

=head2 has_failed_subsnps

  Description : DEPRECATED. Use has_failed_alleles on the Variation object instead.
  Exceptions  : Thrown on invocation.
  Status      : DEPRECATED

=cut

sub has_failed_subsnps {
    my $self = shift;
  
    throw("The has_failed_subsnps subroutine in VariationAdaptor is deprecated. Use has_failed_alleles on the Variation object instead");
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
  ReturnType  : reference to a list of strings
  Exceptions  : none
  Caller      : general
  Status      : At Risk

=cut

sub get_all_failed_descriptions {
    my $self = shift;
    my $variation = shift;
    
    #ÊCall the internal get method without any constraints
    my $description = $self->_internal_get_failed_descriptions($variation) || [];
    
    return $description;
}


=head2 get_failed_description

  Description : DEPRECATED. The appropriate subroutine on the Variation/Allele object should be used instead.
  Exceptions  : Thrown on invocation.
  Status      : DEPRECATED

=cut

sub get_failed_description {
    my $self = shift;
    
    throw("The get_failed_description subroutine in VariationAdaptor is deprecated. Use the appropriate subroutine on the Variation/Allele object instead");
}

#ÊAPI-internal method for getting failed descriptions for a variation
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
  my $self = shift;
  my $sth = shift;

  my ($var_id, $name, $vstatus, $class_attrib_id, $source, $source_desc, $source_url, $source_type, $is_somatic, 
      $ancestral_allele, $allele_id, $allele_ss_id, $allele, $allele_freq, $allele_count, $allele_sample_id, 
      $moltype, $syn_name, $syn_source, $cur_allele_id, $cur_var, $cur_var_id, $flank_flag);

  $sth->bind_columns(\$var_id, \$name, \$vstatus, \$class_attrib_id, \$source, \$source_desc, \$source_url, \$source_type, 
										 \$is_somatic, \$ancestral_allele, \$allele_id, \$allele_ss_id, \$allele, \$allele_freq, \$allele_count,
                     \$allele_sample_id, \$moltype, \$syn_name, \$syn_source, \$flank_flag);

  my @vars;

  my %seen_syns;
  my %seen_pops;

  my $pa = $self->db()->get_PopulationAdaptor();

  my $aa = $self->db->get_AttributeAdaptor;

  while($sth->fetch()) {
    if(!defined($cur_var) || $cur_var_id != $var_id) {
	$vstatus = 0 if (!defined $vstatus);
      my @states = split(',',$vstatus);
      $cur_var = Bio::EnsEMBL::Variation::Variation->new
        (-dbID   => $var_id,
         -ADAPTOR => $self,
         -NAME   => $name,
         -SOURCE => $source,
		 -SOURCE_DESCRIPTION => $source_desc,
		 -SOURCE_URL => $source_url,
		 -SOURCE_TYPE => $source_type,
		 -IS_SOMATIC => $is_somatic,
	     -ANCESTRAL_ALLELE => $ancestral_allele,
	     -MOLTYPE => $moltype,
         -VALIDATION_STATES => \@states,
	     -FLANK_FLAG => $flank_flag,
	     -CLASS_SO_TERM => $aa->attrib_value_for_id($class_attrib_id));
      push @vars, $cur_var;
      $cur_var_id = $var_id;
    }

    if(!defined($cur_allele_id) || $cur_allele_id != $allele_id) {
	  my $pop;
	  if($allele_sample_id) { 
		$pop = $seen_pops{$allele_sample_id} ||=
		  $pa->fetch_by_dbID($allele_sample_id);
	  }
	  
      if (defined $allele_id){
		my $allele = Bio::EnsEMBL::Variation::Allele->new
			(-dbID      => $allele_id,
			 -ADAPTOR   => $self,
			 -ALLELE    => $allele,
			 -FREQUENCY => $allele_freq,
			 -COUNT     => $allele_count,
			 -POPULATION => $pop,
			 -SUBSNP    => $allele_ss_id);
			
		$cur_var->_add_Allele($allele);
		
		$cur_allele_id = $allele_id;
      }
    }

    if(defined ($syn_source) && !$seen_syns{"$syn_source:$syn_name"}) {
      $seen_syns{"$syn_source:$syn_name"} = 1;
      $cur_var->add_synonym($syn_source, $syn_name);
    }
  }

	
  return \@vars;
}

1;
