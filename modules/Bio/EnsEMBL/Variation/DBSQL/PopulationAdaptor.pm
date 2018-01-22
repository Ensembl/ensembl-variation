
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
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::PopulationAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::PopulationAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $pa = $reg->get_adaptor("human","variation","population");

  # Get a Population by its internal identifier
  $pop = $pa->fetch_by_dbID(145);

  # fetch a population by its name
  $pop = $pa->fetch_by_name('PACIFIC'); 

  # fetch all sub populations of a population
  foreach $sub_pop (@{$pa->fetch_all_by_super_Population($pop)}) {
    print $sub_pop->name(), " is a sub population of ", $pop->name(), "\n";
  }

  # fetch all super populations
  foreach $super_pop (@{$pa->fetch_all_by_sub_Population($pop)}) {
    print $pop->name(), " is a sub population of ", $super_pop->name(), "\n";
  }


=head1 DESCRIPTION

This adaptor provides database connectivity for Population objects.
Populations may be retrieved from the Ensembl variation database by
several means using this module.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::PopulationAdaptor;

use Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);
use Bio::EnsEMBL::Utils::Scalar qw(wrap_array);

use Bio::EnsEMBL::Variation::Population;
use DBI qw(:sql_types);
our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor');

sub store {
  my ($self, $pop) = @_;

  my $dbh = $self->dbc->db_handle;

  my $sth = $dbh->prepare(q{
    INSERT INTO population (
      name,
      size,
      description,
      collection,
      freqs_from_gts,
      display
    ) VALUES (?,?,?,?,?,?)
  });

  $sth->execute(
    $pop->name,
    $pop->size,
    $pop->description,
    $pop->collection || 0,
    $pop->_freqs_from_gts || 0,
    $pop->display,
  );
  $sth->finish;

  # get the population_id inserted
  my $dbID = $dbh->last_insert_id(undef, undef, 'population', 'population_id');

  $pop->{dbID}    = $dbID;
  $pop->{adaptor} = $self;
	
}

=head2 fetch_population_by_synonym

    Arg [1]              : String $population_synonym
    Arg [2]              : String $source (optional)
    Example              : my $pop = $pop_adaptor->fetch_population_by_synonym($population_synonym, $source);
    Description          : Retrieves populations for the synonym given in the source. If no source is provided, retrieves all the synonyms
    Returntype           : listref of Bio::EnsEMBL::Variation::Population objects
    Exceptions           : none
    Caller               : general
    Status               : Stable

=cut

sub fetch_population_by_synonym {
  my $self = shift;
  my $synonym_name = shift;
  my $source = shift;
  my ($populations, $population_ids, $population_id);
  my $sql;
  if (defined $source) {
    $sql = qq{
      SELECT ps.population_id 
      FROM population_synonym ps, source s
      WHERE ps.name = ? and ps.source_id = s.source_id AND s.name = "$source"
    };
  }
  else {
    $sql = qq{SELECT population_id FROM population_synonym WHERE name = ?};
  }
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $synonym_name, SQL_VARCHAR);
  $sth->execute();    
  $sth->bind_columns(\$population_id);
  while ($sth->fetch()) {
    push @{$population_ids}, $population_id;
  }

  foreach my $population_id (@{$population_ids}) {
    my $population = $self->fetch_by_dbID($population_id);
    push @{$populations}, $population;
  }
  return $populations;
}

=head2 fetch_synonyms

    Arg [1]              : $population_id
    Arg [2] (optional)   : $source
    Example              : my $dbSNP_synonyms = $pop_adaptor->fetch_synonyms($population_id, $dbSNP);
                           my $all_synonyms = $pop_adaptor->fetch_synonyms($population_id);
    Description: Retrieves synonyms for the source provided. Otherwise, return all the synonyms for the population_id
    Returntype : listref of strings
    Exceptions : none
    Caller     : general
    Status     : Stable

=cut

sub fetch_synonyms {
  my $self = shift;
  my $dbID = shift;
  my $source = shift;
  my $synonym;
  my $synonyms;

  my $sql;
  if (defined $source) {
    $sql = qq{SELECT ps.name FROM population_synonym ps, source s WHERE ps.population_id = ? AND ps.source_id = s.source_id AND s.name = "$source"};
  } else {
    $sql = qq{SELECT name FROM population_synonym WHERE population_id = ?};
  }
  my $sth = $self->prepare($sql);
  $sth->bind_param(1,$dbID,SQL_INTEGER);
  $sth->execute();
  $sth->bind_columns(\$synonym);
  while ($sth->fetch) {
    push @{$synonyms}, $synonym;
  }
  return $synonyms;
}

=head2 fetch_by_name

  Arg [1]    : string $name
  Example    : $pop = $pop_adaptor->fetch_by_name('NUSPAE:Singapore_HDL');
  Description: Retrieves a population object via its name
  Returntype : Bio::EnsEMBL::Variation::Population
  Exceptions : throw if name argument is not defined
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_name {
  my $self = shift;
  my $name = shift;

  throw('Name argument expected.') if (!defined($name));

  my $sth = $self->prepare(q{
    SELECT p.population_id, p.name, p.size, p.description, p.collection, p.freqs_from_gts, p.display, dg.display_name, dg.display_priority
    FROM   population p left outer join display_group dg on dg.display_group_id = p.display_group_id
    WHERE  p.name = ?;
  });

  $sth->bind_param(1, $name, SQL_VARCHAR);
  $sth->execute();
  my $populations = $self->_objs_from_sth($sth);
  $sth->finish();
  return undef if (!@$populations);
  return $populations->[0];
}

=head2 fetch_all_by_dbID_list

  Arg [1]    : listref of dbIDs
  Example    : $pops = $pop_adaptor->fetch_all_by_dbID_list([907, 1132]);
  Description: Retrieves a listref of population objects via a list of internal
               dbID identifiers
  Returntype : listref of Bio::EnsEMBL::Variation::Population objects
  Exceptions : throw if list argument is not defined
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_dbID_list {
  my $self = shift;
  my $list = shift;

  if (!defined($list) || ref($list) ne 'ARRAY') {
    throw("list reference argument is required");
  } 

  return [] unless scalar @$list >= 1;
  my $id_str = (@$list > 1)  ? " IN (".join(',',@$list).")"   :   ' = \''.$list->[0].'\'';
  my $sth = $self->prepare(qq{
    SELECT p.population_id, p.name, p.size, p.description, p.collection, p.freqs_from_gts, p.display, dg.display_name, dg.display_priority
    FROM population p left outer join display_group dg on dg.display_group_id = p.display_group_id
    WHERE  p.population_id $id_str;
  });
  $sth->execute();
  my $populations = $self->_objs_from_sth($sth);
  $sth->finish();
  return undef if(!@$populations);
  return $populations;
}

=head2 fetch_all_by_name_search

  Arg [1]    : string $name
  Example    : $pop = $pop_adaptor->fetch_all_by_name_search('CEU');
  Description: Retrieves a list of population objects whose name matches the
               search term.
  Returntype : Listref of Bio::EnsEMBL::Variation::Population objects
  Exceptions : throw if name argument is not defined
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_name_search {
  my $self = shift;
  my $name = shift;

  throw('Name argument expected.') if(!defined($name));

  my $sth = $self->prepare(q{
    SELECT p.population_id, p.name, p.size, p.description, p.collection, p.freqs_from_gts, p.display, dg.display_name, dg.display_priority
    FROM   population p left outer join display_group dg on dg.display_group_id = p.display_group_id
    WHERE  p.name like concat('%', ?, '%');
  });

  $sth->bind_param(1,$name,SQL_VARCHAR);
  $sth->execute();
  my $populations = $self->_objs_from_sth($sth);
  $sth->finish();
  return $populations;
}


=head2 fetch_all_by_super_Population

  Arg [1]    : Bio::EnsEMBL::Variation::Population $pop
  Example    : foreach $sub_pop (@{$pa->fetch_all_by_super_Population($pop)}) {
                 print $sub_pop->name(), "\n";
               }
  Description: Retrieves all sub populations of a provided population.
  Returntype : listref of Bio::EnsEMBL::Variation::Population objetcs
  Exceptions : throw on bad argument
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_super_Population {
  my $self = shift;
  my $pop  = shift;

  if (!ref($pop) || !$pop->isa('Bio::EnsEMBL::Variation::Population')) {
    throw('Bio::EnsEMBL::Variation::Population argument expected');
  }

  if (!$pop->dbID()) {
    warning("Cannot retrieve sub populations for population without dbID");
    return [];
  }

  my $sth = $self->prepare(q{
    SELECT p.population_id, p.name, p.size, p.description, p.collection, p.freqs_from_gts, p.display, dg.display_name, dg.display_priority
    FROM population_structure ps,  population p
    LEFT OUTER JOIN display_group dg on p.display_group_id = dg.display_group_id
    WHERE  p.population_id = ps.sub_population_id
    AND    ps.super_population_id = ?;
  });

  $sth->bind_param(1,$pop->dbID,SQL_INTEGER);
  $sth->execute();
  my $populations = $self->_objs_from_sth($sth);
  $sth->finish();
  return $populations;
}


=head2 fetch_all_by_sub_Population

  Arg [1]    : Bio::EnsEMBL::Variation::Population $pop
  Example    : foreach $super_pop (@{$pa->fetch_all_by_sub_Population($pop)}) {
                 print $super_pop->name(), "\n";
               }
  Description: Retrieves all super populations for a provided population
  Returntype : listref of Bio::EnsEMBL::Variation::Population objects
  Exceptions : throw on bad argument
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_sub_Population {
  my $self = shift;
  my $pop  = shift;

  if (!ref($pop) || !$pop->isa('Bio::EnsEMBL::Variation::Population')) {
    throw('Bio::EnsEMBL::Variation::Population argument expected');
  }

  if (!$pop->dbID()) {
    warning("Cannot retrieve super populations for population without dbID");
    return [];
  }

  my $sth = $self->prepare(q{
    SELECT p.population_id, p.name, p.size, p.description, p.collection, p.freqs_from_gts, p.display, dg.display_name, dg.display_priority
    FROM   population_structure ps, population p
    LEFT OUTER JOIN display_group dg on dg.display_group_id = p.display_group_id
    WHERE  p.population_id = ps.super_population_id
    AND    ps.sub_population_id = ?;
  });

  $sth->bind_param(1,$pop->dbID,SQL_INTEGER);
  $sth->execute();
  my $populations = $self->_objs_from_sth($sth);
  $sth->finish();
  return $populations;
}


=head2 fetch_default_LDPopulation

    Args        : none
    Example     : $population = $pop_adaptor->fetch_default_LDPopulation();
    Description : Obtains the population it is used as a default in the LD display of the pairwise LD data
    ReturnType  : Bio::EnsEMBL::Variation::Population
    Exceptions  : none
    Caller      : general
    Status      : At Risk

=cut

sub fetch_default_LDPopulation {
  my $self = shift;
  my $population_id;

  my $sth = $self->prepare(qq{SELECT meta_value from meta where meta_key = ?});
  $sth->bind_param(1,'pairwise_ld.default_population',SQL_VARCHAR);
  $sth->execute();
  $sth->bind_columns(\$population_id);
  $sth->fetch();
  $sth->finish;

  if (defined $population_id) {
    return $self->fetch_by_dbID($population_id);
  } else {
    return undef;
  }
}

=head2 fetch_all_vcf_Populations

    Example     : @populations = @{$pop_adaptor->fetch_all_vcf_Populations();
    Description : Gets all populations that are represented in VCF files
    ReturnType  : reference to list of Bio::EnsEMBL::Variation::Population objects
    Exceptions  : none
    Caller      : general
    Status      : stable

=cut

sub fetch_all_vcf_Populations {
  my $self = shift;
  my $use_vcf = $self->db->use_vcf();
  if (!$use_vcf) {
    warning('You need to set use_vcf: $sample_genotype_feature_adaptor->db->use_vcf(1)');
    return [];
  }
  my @vcf_pops = map {@{$_->get_all_Populations}} @{$self->db->get_VCFCollectionAdaptor->fetch_all};
  return \@vcf_pops;
}


=head2 fetch_all_LD_Populations

    Example     : @populations = @{$pop_adaptor->fetch_all_LD Populations();
    Description : Gets all populations that can be used in the LD display
    ReturnType  : reference to list of Bio::EnsEMBL::Variation::Population objects
    Exceptions  : none
    Caller      : general
    Status      : At Risk

=cut

sub fetch_all_LD_Populations {
  my $self = shift;
#   return $self->_fetch_all_db_LD_Populations;
#   my $use_vcf = $self->db->use_vcf;
  
#   if ($use_vcf) {
#     my @vcf_pops = map {@{$_->get_all_Populations}} @{$self->db->get_VCFCollectionAdaptor->fetch_all};
    
#     # value of >1 means use only VCF pops
#     if ($use_vcf > 1) {
#       return \@vcf_pops;
#     }
    
#     # value of 1 means use both VCF and DB pops
#     else {
#       my @db_pops = @{$self->_fetch_all_db_LD_Populations};
      
#       # merge populations
#       my %merged = map {$_->name() => $_} (@vcf_pops, @db_pops);
      
#       return [values %merged];
#     }
#   }
  
#   else {
#     return $self->_fetch_all_db_LD_Populations;
#   }
# }

# sub _fetch_all_db_LD_Populations {
#   my $self = shift;
  return [grep {$_->name !~ /ALL|AFR|AMR|ASN|EUR/} @{$self->generic_fetch(qq{ p.display = 'LD' })}];
}


=head2 fetch_all_HapMap_Populations

    Example     : @populations = @{$pop_adaptor->fetch_all_HapMap_populations();
    Description : Gets all populations from the HapMap project (human only!)
    ReturnType  : reference to list of Bio::EnsEMBL::Variation::Population objects
    Exceptions  : none
    Caller      : general
    Status      : At Risk

=cut

sub fetch_all_HapMap_Populations {
  my $self = shift;
  return $self->generic_fetch(qq{ p.name like 'cshl-hapmap%' });
}


=head2 fetch_all_1KG_Populations

    Example     : @populations = @{$pop_adaptor->fetch_all_1KG_populations();
    Description : Gets all populations from the 1000 genomes project (human only!)
    ReturnType  : reference to list of Bio::EnsEMBL::Variation::Population objects
    Exceptions  : none
    Caller      : general
    Status      : At Risk

=cut

sub fetch_all_1KG_Populations{
  my $self = shift;
  return $self->generic_fetch(qq{ p.name like '1000GENOMES%' });
}


=head2 fetch_all_by_Individual

  Arg [1]     : Bio::EnsEMBL::Variation::Individual $ind
  Example     : my $ind = $ind_adaptor->fetch_by_name('NA12004');
                foreach my $pop (@{$pop_adaptor->fetch_all_by_Individual($ind)}){
		          print $pop->name, "\n";
                }
  Description : Retrieves all populations from a specified individual
  ReturnType  : reference to list of Bio::EnsEMBL::Variation::Population objects
  Exceptions  : throw if incorrect argument is passed
                warning if provided individual does not have a dbID
  Caller      : general
  Status      : Stable

=cut

sub fetch_all_by_Individual {
  my $self = shift;
  my $ind = shift;

  if (!ref($ind) || !$ind->isa('Bio::EnsEMBL::Variation::Individual')) {
    throw("Bio::EnsEMBL::Variation::Individual arg expected");
  }

  if (!$ind->dbID()) {
    warning("Individual does not have dbID, cannot retrieve Individuals");
    return [];
  } 

  my $sth = $self->prepare(qq{
    SELECT p.population_id, p.name, p.size, p.description, p.collection, p.freqs_from_gts, p.display, dg.display_name, dg.display_priority
    FROM sample_population sp, sample s, population p
    LEFT OUTER JOIN display_group dg on p.display_group_id = dg.display_group_id
    WHERE p.population_id = sp.population_id
    AND sp.sample_id = s.sample_id
    AND s.individual_id = ?;
  });
  $sth->bind_param(1,$ind->dbID,SQL_INTEGER);
  $sth->execute();
  my $populations = $self->_objs_from_sth($sth);
  $sth->finish();
  return $populations;
}

=head2 fetch_all_by_Sample

  Arg [1]     : Bio::EnsEMBL::Variation::Sample $sample
  Example     : my $ind = $ind_adaptor->fetch_by_name('NA12004');
                foreach my $pop (@{$pop_adaptor->fetch_all_by_Sample($sample)}){
                  print $pop->name, "\n";
                }
  Description : Retrieves all populations from a specified sample
  ReturnType  : reference to list of Bio::EnsEMBL::Variation::Population objects
  Exceptions  : throw if incorrect argument is passed
                warning if provided sample does not have a dbID
  Caller      : general
  Status      : Stable

=cut

sub fetch_all_by_Sample {
  my $self = shift;
  my $sample = shift;

  if (!ref($sample) || !$sample->isa('Bio::EnsEMBL::Variation::Sample')) {
    throw("Bio::EnsEMBL::Variation::Sample arg expected");
  }

  if (!$sample->dbID()) {
    warning("Sample does not have dbID, cannot retrieve Populations");
    return [];
  } 

  my $sth = $self->prepare(qq{
    SELECT p.population_id, p.name, p.size, p.description, p.collection, p.freqs_from_gts, p.display, dg.display_name, dg.display_priority
    FROM sample_population sp, population p
    LEFT OUTER JOIN display_group dg on p.display_group_id = dg.display_group_id
    WHERE p.population_id = sp.population_id
    AND sp.sample_id = ?;
  });
  $sth->bind_param(1, $sample->dbID, SQL_INTEGER);
  $sth->execute();
  my $populations = $self->_objs_from_sth($sth);
  $sth->finish();
  return $populations;
}

=head2 fetch_all_by_Individual_list

  Arg [1]     : listref of of Bio::EnsEMBL::Variation::Individual objects
  Example     : foreach my $pop (@{$pop_adaptor->fetch_all_by_Individual_list($inds)}){
		          print $pop->name,"\n";
                }
  Description : Retrieves all populations from a specified individual
  ReturnType  : reference to list of Bio::EnsEMBL::Variation::Population objects
  Exceptions  : throw if incorrect argument is passed
                warning if provided individual does not have a dbID
  Caller      : general
  Status      : Stable

=cut

sub fetch_all_by_Individual_list {
	my $self = shift;
	my $list = shift;
	
	if (!ref($list) || !$list->[0]->isa('Bio::EnsEMBL::Variation::Individual')) {
		throw("Listref of Bio::EnsEMBL::Variation::Individual arg expected");
	}
	
	if (!$list->[0]->dbID()) {
		warning("First Individual does not have dbID, cannot retrieve Populations");
		return [];
	}
	
	my $id_str = " IN (" . join(',', map {$_->dbID} @$list). ")";	
	
	my $sth = $self->prepare(qq{
		SELECT p.population_id, p.name, p.size, p.description, p.collection, p.freqs_from_gts, p.display, dg.display_name, dg.display_priority
    FROM sample_population sp, population p
    LEFT OUTER JOIN display_group dg on dg.display_group_id = p.display_group_id
		WHERE p.population_id = sp.population_id
		AND sp.sample_id $id_str
	});
	$sth->execute();
	my $populations = $self->_objs_from_sth($sth);
	$sth->finish();
	return $populations;
}

=head2 fetch_all_by_Sample_list

  Arg [1]     : listref of of Bio::EnsEMBL::Variation::Sample objects
  Example     : foreach my $pop (@{$pop_adaptor->fetch_all_by_Sample_list($samples)}){
                  print $pop->name,"\n";
                }
  Description : Retrieves all populations from a list of samples 
  ReturnType  : reference to list of Bio::EnsEMBL::Variation::Population objects
  Exceptions  : throw if incorrect argument is passed
                warning if provided samples do not have a dbIDs
  Caller      : general
  Status      : Stable

=cut

sub fetch_all_by_Sample_list {
	my $self = shift;
	my $list = shift;
	
	if (!ref($list) || !$list->[0]->isa('Bio::EnsEMBL::Variation::Sample')) {
		throw("Listref of Bio::EnsEMBL::Variation::Sample arg expected");
	}
	
	if (!$list->[0]->dbID()) {
		warning("First Sample does not have a dbID, cannot retrieve Populations");
		return [];
	}
	
	my $id_str = " IN (" . join(',', map {$_->dbID} @$list). ")";	
	
	my $sth = $self->prepare(qq{
		SELECT p.population_id, p.name, p.size, p.description, p.collection, p.freqs_from_gts, p.display, dg.display_name, dg.display_priority
    FROM sample_population sp, population p
    LEFT OUTER JOIN display_group dg on dg.display_group_id = p.display_group_id
		WHERE p.population_id = sp.population_id
		AND sp.sample_id $id_str
	});
	$sth->execute();
	my $populations = $self->_objs_from_sth($sth);
	$sth->finish();
	return $populations;
}


=head2 get_dbIDs_for_population_names

  Arg [1]     : $population_names
                Listref of population names.
  Example     : my $ids = $pop_adaptor->get_dbIDs_for_population_names(['CSHL-HAPMAP:HAPMAP-MEX','1000GENOMES:pilot_1_CHB+JPT_low_coverage_panel']);
                map {printf("Population: \%s has dbID \%d\n",$ids->{$_},$_)} keys(%{$ids});
  Description : Retrieve the dbIDs for a list of population names
  ReturnType  : hashref with dbIDs as keys and population names as values
  Caller      : web
  Status      : At Risk

=cut

sub get_dbIDs_for_population_names {
  my $self = shift;
  my $population_names = shift;

  # Wrap the argument into an arrayref
  $population_names = wrap_array($population_names);

  # Define a statement handle for the lookup query
  my $stmt = qq{
    SELECT population_id, name
    FROM population
    WHERE name = ?
    LIMIT 1
  };
  my $sth = $self->prepare($stmt);

  # Loop over the population names and query the db
  my %dbIDs;
  foreach my $pop_name (@{$population_names}) {
    $sth->execute($pop_name);
    my ($id, $name);
    $sth->bind_columns(\$id,\$name);
    $sth->fetch;
    $dbIDs{$id} = $name if (defined($id));
  }

  return \%dbIDs;
}


sub _get_sample_population_hash {
	my $self = shift;
	my $id_list_ref = shift;

	if(!defined($id_list_ref) || ref($id_list_ref) ne 'ARRAY') {
		throw("id_list list reference argument is required");
	}
	
	return {} if (!@$id_list_ref);
	
	my %ip_hash;
	my $max_size = 200;
	my @id_list = @$id_list_ref;
	
	while (@id_list) {
		my @ids;
		if(@id_list > $max_size) {
			@ids = splice(@id_list, 0, $max_size);
		} else {
			@ids = splice(@id_list, 0);
		}
		
		my $id_str;
		if(@ids > 1) {
			$id_str = " IN (" . join(',', @ids). ")";
		} else {
			$id_str = " = ".$ids[0];
		}
		
		my $sth = $self->prepare(qq/
			SELECT sample_id, population_id
			FROM sample_population
			WHERE sample_id $id_str
		/);
		
		$sth->execute();
		
		my ($sample, $pop);
		$sth->bind_columns(\$sample, \$pop);
		$ip_hash{$pop}{$sample} = 1 while $sth->fetch;
		$sth->finish();
	}
	
	# NB COMMENTED OUT FOR NOW AS IT DOESN'T SEEM TO WORK PROPERLY
	# now get super-populations
	#my @pops = keys %ip_hash;
	#my %new_pops;
	#
	## need to iterate in case there's multiple levels
	#while(scalar @pops) {
	#	
	#	my $id_str;
	#	if(scalar @pops)  {
	#		$id_str = " IN (" . join(',', @pops). ")";
	#	} else {
	#		$id_str = " = ".$pops[0];
	#	}
	#	
	#	@pops = ();
	#	
	#	my $sth = $self->prepare(qq{
	#		SELECT sub_population_sample_id, super_population_sample_id
	#		FROM population_structure
	#		WHERE sub_population_sample_id $id_str
	#	});
	#	$sth->execute();
	#	
	#	my ($sub, $super);
	#	$sth->bind_columns(\$sub, \$super);
	#	while($sth->fetch) {
	#		push @{$new_pops{$sub}}, $super;
	#		push @pops, $super;
	#	}
	#	$sth->finish();
	#}
	#
	#foreach my $sub(keys %new_pops) {
	#	foreach my $super(@{$new_pops{$sub}}) {
	#		$ip_hash{$super}{$_} = 1 for keys %{$ip_hash{$sub}};
	#	}
	#}
	#
	return \%ip_hash;
}
#
# private method, creates population objects from an executed statement handle
# ordering of columns must be consistant
#
sub _objs_from_sth {
  my $self = shift;
  my $sth  = shift;

  my @pops;

  my ($pop_id, $name, $size, $desc, $collection, $freqs, $display, $display_name, $display_priority);

  $sth->bind_columns(\$pop_id, \$name, \$size, \$desc, \$collection, \$freqs, \$display, \$display_name, \$display_priority);

  while($sth->fetch()) {
    push @pops, Bio::EnsEMBL::Variation::Population->new(
        -dbID => $pop_id,
        -ADAPTOR => $self,
        -NAME => $name,
        -SIZE => $size,
        -DESCRIPTION => $desc,
        -COLLECTION => $collection,
        -FREQS => $freqs,
        -DISPLAY => $display,
        -DISPLAY_GROUP_NAME => $display_name,
        -DISPLAY_GROUP_PRIORITY => $display_priority,
    );
  }
  return \@pops;
}

sub _tables {
  return (['population','p'], ['display_group', 'dg']);
}

sub _columns {
    return qw(p.population_id p.name p.size p.description p.collection p.freqs_from_gts p.display dg.display_name dg.display_priority);
}

sub _left_join {
  my $self = shift;

  my @left_join = (
    ['display_group dg', 'p.display_group_id = dg.display_group_id'],
  );       

  return @left_join;
}

sub _default_where_clause {
  return '';
}

1;
