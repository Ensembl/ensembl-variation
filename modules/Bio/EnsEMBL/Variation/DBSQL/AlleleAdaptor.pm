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
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::AlleleAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::AlleleAdaptor

=head1 SYNOPSIS

  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_registry_from_db(-host => 'ensembldb.ensembl.org', -user => 'anonymous');

  my $aa = $registry->get_adaptor('human', 'variation', 'allele');
  my $va = $registry->get_adaptor('human', 'variation', 'variation');

  my $variation = $va->fetch_by_name('rs678');
  my $alleles = $aa->fetch_all_by_Variation($variation);
  foreach my $allele (@$alleles) {
    my $allele_string = $allele->allele();
    my $population_name = ($allele->population()) ? $allele->population()->name : 'population is NA';
    my $frequency = ($allele->frequency()) ? $allele->frequency : 'frequency is NA';
    print join(' ', $allele_string, $population_name, $frequency), "\n";
  }

=head1 DESCRIPTION

This adaptor provides database connectivity for Allele objects.
Alleles may be retrieved from the Ensembl variation database by
several means using this module.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::AlleleAdaptor;

use Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);

use Bio::EnsEMBL::Variation::Allele;
use Bio::EnsEMBL::Utils::Iterator;

use DBI qw(:sql_types);

our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor');

my $DEFAULT_ITERATOR_CACHE_SIZE = 1000;

our $CACHE_SIZE = 5;


sub store {
	my ($self, $allele) = @_;
	
	my $dbh = $self->dbc->db_handle;
	
	# look up or store allele code
	my $allele_code = $self->_allele_code($allele->allele);
	
	my $sth = $dbh->prepare_cached(q{
		INSERT DELAYED INTO allele (
			variation_id,
			subsnp_id,
			allele_code_id,
			population_id,
			frequency,
			count			
		) VALUES (?,?,?,?,?,?)
	});
	
	$sth->execute(
		$allele->{_variation_id} || $allele->variation->dbID,
		$allele->{subsnp},
		$allele_code,
		$allele->population ? $allele->population->dbID : undef,
		$allele->frequency,
		$allele->count
	);
  
  # reset cache
  delete $self->{_cache} if $self->{_cache};
	
	$sth->finish;
}

sub store_multiple {
	my ($self, $alleles) = @_;
	
	my $dbh = $self->dbc->db_handle;
	
	my $q_string = join ",", map {'(?,?,?,?,?,?)'} @$alleles;
	
	my @args = map {
		$_->{_variation_id} || $_->variation->dbID,
		$_->{subsnp},
		$self->_allele_code($_->allele),
		$_->population ? $_->population->dbID : undef,
		$_->frequency,
		$_->count
	} @$alleles;
	
	my $sth = $dbh->prepare_cached(qq{
		INSERT INTO allele (
			variation_id,
			subsnp_id,
			allele_code_id,
			population_id,
			frequency,
			count			
		) VALUES $q_string
	});
	
	$sth->execute(@args);
  
  # reset cache
  delete $self->{_cache} if $self->{_cache};
	
	$sth->finish;
}

sub store_to_file_handle {
	my ($self, $allele, $file_handle) = @_;
	
	my $dbh = $self->dbc->db_handle;
	
	print $file_handle join("\t",
		$allele->{_variation_id} || $allele->variation->dbID || '\N',
		$allele->{subsnp} || '\N',
		$self->_allele_code($allele->allele),
		$allele->population ? $allele->population->dbID : '\N',
		defined($allele->frequency) ? $allele->frequency :  '\N',
		defined($allele->count) ? $allele->count : '\N',
	)."\n";
}

=head2 fetch_all

  Description: fetch_all should not be used for Alleles.
  Exceptions : thrown on invocation
  Status     : Stable

=cut

sub fetch_all {
    my $self = shift;

    throw("fetch_all cannot be used for Allele objects");
}

=head2 fetch_all_by_subsnp_id

  Arg [1]    : string $subsnp_id
  Example    : $alleles = $allele_adaptor->fetch_all_by_subsnp_id('ss123');
  Description: Retrieves all allele objects via a component subsnp ID
  Returntype : listref of Bio::EnsEMBL::Variation::Allele objects
  Exceptions : throw if name argument is not defined
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_subsnp_id {
  my $self = shift;
  my $name = shift;

  throw('name argument expected') if(!defined($name));
  $name =~ s/^ss//gi;

  # Add the constraint on the subsnp_id column and pass to generic_fetch
  my $constraint = qq{ a.subsnp_id = $name };

  return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_Variation

  Arg [1]    : Bio::EnsEMBL::Variation::Variation
  Arg [2]    : Bio::EnsEMBL::Variation::Population (optional)
  Example    : @alleles = @{$allele_adaptor->fetch_all_by_Variation($var)};
  Description: Retrieves all alleles which are associated with a specified
               variation. If the optional population argument is specified, only retrieve
               alleles for that population.
  Returntype : listref of Bio::EnsEMBL::Variation::Allele
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Variation {
  my $self = shift;
  my $variation = shift;
  my $population = shift;
  
  # Make sure that we are passed a Variation object
  assert_ref($variation,'Bio::EnsEMBL::Variation::Variation');
  
  # If we got a population argument, make sure that it is a Population object
  assert_ref($population,'Bio::EnsEMBL::Variation::Population') if (defined($population));
  
  my $variation_id = $variation->dbID();
  
  my $cached;
  my $return = [];

  if(defined($self->{_cache})) {
    foreach my $stored(@{$self->{_cache}}) {
      my @keys = keys %{$stored};
      $cached = $stored->{$keys[0]} if $keys[0] eq $variation_id;
      last if defined($cached);
    }
  }

  if(!defined($cached)) {

    my $use_vcf = $self->db->use_vcf() || 0;
    my @from_vcf;

    if($use_vcf) {
      my $vfs = $variation->get_all_VariationFeatures;

      if($vfs && @$vfs) {
        @from_vcf =
          map {$_->{adaptor} = $self; $_}
          map {@{$_->get_all_Alleles_by_VariationFeature($vfs->[0])}}
          @{$self->db->get_VCFCollectionAdaptor->fetch_all() || []};
      }
    }

    push @$cached, @from_vcf;
    # return $cached;

    # use_vcf == 2 specifies we don't want anything from the DB
    unless($use_vcf == 2) {
      
      # Add a constraint on the variation_id column and pass to generic fetch
      my $constraint = qq{ a.variation_id = $variation_id };
      
      # If required, add a constraint on the population id
      if (defined($population)) {
        my $population_id = $population->dbID();
        $constraint .= qq{ AND a.population_id = $population_id };
      }
      
      # Add the constraint for failed alleles
      $constraint .= " AND " . $self->db->_exclude_failed_alleles_constraint();
    
      push @$cached, @{$self->generic_fetch($constraint)};
      
      # If a population was specified, attach the population to the object
      map {$_->population($population)} @{$cached} if (defined($population));
    }

    # add freqs from genotypes for human (1KG data), could be from VCF so don't exclude here
    push @$cached, @{$self->_fetch_all_by_Variation_from_Genotypes($variation, $population)};
		
    # don't store if population specified
    return $cached if defined($population);
    
    # add genotypes for this variant to the cache
    push @{$self->{_cache}}, {$variation_id => $cached};
	
    # shift off first element to keep cache within size limit
    shift @{$self->{_cache}} if scalar @{$self->{_cache}} > $CACHE_SIZE;
  }
  
  if (defined($population)) {
    my $population_id = $population->dbID;
    foreach (@{$cached}) {
      my $allele_population = $_->population;
      if ($allele_population) {
        push @$return, $_ if ($allele_population->dbID eq $population_id);
      }
    }
  }
  else {
    $return = $cached;
  }
	
  return $return;
}

sub _fetch_all_by_Variation_from_Genotypes {
  my $self = shift;
  return $self->_generic_fetch_all_by_Variation_from_Genotypes(@_, 'Allele');
}

=head2 get_all_failed_descriptions

  Arg[1]      : Bio::EnsEMBL::Variation::Allele
	               The allele object to get the failed descriptions for
  Example     : 
                my $failed_descriptions = $adaptor->get_all_failed_descriptions($allele);
                if (scalar(@{$failed_descriptions})) {
		          print "The allele '" . $allele->allele() . "' has been flagged as failed because '" . join("' and '",@{$failed_descriptions}) . "'\n";
                }
		
  Description : Gets the unique descriptions for the reasons why the supplied allele has failed.
  ReturnType  : reference to list of strings
  Exceptions  : thrown on incorrect argument
  Caller      : general
  Status      : At Risk

=cut

sub get_all_failed_descriptions {
    my $self = shift;
    my $allele = shift;
    
    # Call the internal get method without any constraints
    my $description = $self->_internal_get_failed_descriptions($allele) || [];
    
    return $description;
}

=head2 get_subsnp_handle

  Arg[1]      : Bio::EnsEMBL::Variation::Allele
	               The allele object to get the observation submitter subsnp handle for
  Arg[2]      : Bio::EnsEMBL::Variation::Population (optional)
                       The population object to get the frequency submitter subsnp handle for                           
  Example     : 
                my $handle = $adaptor->get_subsnp_handle($allele);
		print "The allele '" . $allele->allele() . "' of subsnp 'ss" . $allele->subsnp() . "' was submitted by '$handle'\n";

                my $handle = $adaptor->get_subsnp_handle($allele, $population);
		print "Frequency data for subsnp 'ss" . $allele->subsnp() . "' in population " . $population->name() . " was submitted by '$handle'\n";


		
  Description : Gets the submitter handle for the specified allele
                The initial observation and later allele frequencies may be submitted for the same ss record by different groups
                supplying a population returns the submitter handle for data obtained within this population  
  ReturnType  : String
  Exceptions  : thrown on incorrect argument
  Caller      : general
  Status      : Stable

=cut

sub get_subsnp_handle {
    my $self = shift;
    my $allele = shift;
    my $population = shift ; ## return submitter for specific allele frequency

    # Assert that the object passed is an Allele
    assert_ref($allele,'Bio::EnsEMBL::Variation::Allele');
    
    # If we got a population argument, make sure that it is a Population object
    assert_ref($population,'Bio::EnsEMBL::Variation::Population') if (defined($population));


    # Get the subsnp id and get rid of any 'ss' prefix
    my $ssid = $allele->subsnp() || "";
    $ssid =~ s/^ss//;
    
    my ($stmt, $sth);
    if(defined $population ){
       $stmt = qq{
        SELECT
            submitter_handle.handle
        FROM
            submitter_handle, allele
        WHERE
            allele.subsnp_id = ?
        AND allele.frequency_submitter_handle = submitter_handle.handle_id
        AND allele.population_id = ?
        LIMIT 1
        };
       $sth = $self->prepare($stmt);
       $sth->execute($ssid, $population->dbID());       
       
    }
    else{
        $stmt = qq{
          SELECT
              handle
          FROM
              subsnp_handle
          WHERE
              subsnp_id = ?
          LIMIT 1
        };


    $sth = $self->prepare($stmt);
    $sth->execute($ssid);

    }
    my $handle;
    $sth->bind_columns(\$handle);
    $sth->fetch();
    return $handle;
}


# API-internal method for getting failed descriptions for an Allele
sub _internal_get_failed_descriptions {
    my $self = shift;
    my $allele = shift;
    my $constraint = shift;
    
    # Assert that the object passed is an Allele
    assert_ref($allele,'Bio::EnsEMBL::Variation::Allele');
    
    my $stmt = qq{
        SELECT DISTINCT
            fd.description
        FROM
            failed_allele fa JOIN
            failed_description fd ON (
                fd.failed_description_id = fa.failed_description_id
            )
        WHERE
            fa.allele_id = ?
    };
    $stmt .= qq{ AND $constraint } if (defined($constraint));
    
    my $sth = $self->prepare($stmt);
    $sth->execute($allele->dbID());
    my @descriptions;
    my $description;
    $sth->bind_columns(\$description);
    while ($sth->fetch()) {
        push(@descriptions,$description);
    }
    return \@descriptions;
}

sub _objs_from_sth {
    my $self = shift;
    my $sth = shift;

    my ($allele_id, $variation_id, $subsnp_id, $allele, $frequency, $population_id, $count, $last_allele_id);
    my @alleles;
    
    $sth->bind_columns(\$allele_id, \$variation_id, \$subsnp_id, \$allele, \$frequency, \$population_id, \$count);
    
    while($sth->fetch()) {
    
        # The left join to failed allele can create duplicate rows, so check that we've got a new Allele before creating the object
        unless (defined($last_allele_id) && $last_allele_id == $allele_id) {
        
            my $obj = Bio::EnsEMBL::Variation::Allele->new(
                -dbID           => $allele_id,
                -VARIATION_ID   => $variation_id,
                -SUBSNP         => $subsnp_id,
                -ALLELE         => $allele,
                -FREQUENCY      => $frequency,
                -POPULATION_ID  => $population_id,
                -COUNT          => $count,
                -ADAPTOR        => $self
            );
              
            push(@alleles,$obj);
            $last_allele_id = $allele_id;
        }
        
    }
    
    return \@alleles;
}

# method used by superclass to construct SQL
sub _tables { 
    my $self = shift;
    
    my @tables = (
        ['allele', 'a'], ['allele_code', 'ac']
    );
    
	# If we are excluding failed_alleles, add that table
	push(@tables,['failed_allele', 'fa']) unless ($self->db->include_failed_variations());
	
	return @tables;
}

# Add a left join to the failed_variation table
sub _left_join { 
    my $self = shift;
    
    # If we are including failed variations, skip the left join
    return () if ($self->db->include_failed_variations());
    return ([ 'failed_allele', 'fa.allele_id = a.allele_id']); 
}

sub _columns {
  return qw( a.allele_id a.variation_id a.subsnp_id ac.allele a.frequency a.population_id a.count );
}

sub _write_columns {
	return qw(variation_id subsnp_id allele_code_id population_id frequency count);
}

sub _default_where_clause  {
  return 'a.allele_code_id = ac.allele_code_id';
}

sub _cache_allele_codes {
	my $self = shift;
	
	my $dbh = $self->dbc->db_handle;
	
	my $sth = $dbh->prepare(qq{SELECT allele_code_id, allele FROM allele_code});
	$sth->execute;
	
	my ($code, $allele);
	$sth->bind_columns(\$code, \$allele);
	my %allele_codes;
	$allele_codes{defined($allele) ? $allele : ''} = $code while $sth->fetch;
	$sth->finish();
	
	$self->db->{_allele_codes} = \%allele_codes;
	
	return $self->db->{_allele_codes};
}

sub _allele_code {
	my ($self, $allele) = @_;
	
	# check if cache is loaded
	my $just_loaded = 0;
	
	if(!exists($self->db->{_allele_codes})) {
		$self->_cache_allele_codes;
		$just_loaded = 1;
	}
	
	if(!exists($self->db->{_allele_codes}->{$allele})) {
		
		# check another process hasn't created it by reloading the cache
		$self->_cache_allele_codes unless $just_loaded;
		
		# if it still doesn't exist
		if(!exists($self->db->{_allele_codes}->{$allele})) {
			my $dbh = $self->dbc->db_handle;
			
			my $sth = $dbh->prepare(q{
				INSERT INTO allele_code (
					allele
				)
				VALUES (?)
			});
			eval {
				$sth->execute($allele);
			};
			$sth->finish;
			
			my $allele_code;
			
			# insert failed, another process did it maybe?
			if($@) {
				my $sth2 = $dbh->prepare(qq{
					SELECT allele_code_id
					FROM allele_code
					WHERE allele = ?
				});
				$sth2->execute($allele);
				$sth2->bind_columns(\$allele_code);
				$sth2->fetch();
				
				throw("ERROR: Failed to insert allele '$allele' into allele_code") unless defined($allele_code);
				
				$sth2->finish();
			}
			else {
				$allele_code = $dbh->last_insert_id(undef, undef, 'allele_code', 'allele_code_id');
			}
			
			$self->db->{_allele_codes}->{$allele} = $allele_code;
		}
	}
	
	return $self->db->{_allele_codes}->{$allele};
}
1;
