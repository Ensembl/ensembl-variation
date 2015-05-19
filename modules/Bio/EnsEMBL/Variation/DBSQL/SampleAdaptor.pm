=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::SampleAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::DBSQL::SampleAdaptor

=head1 SYNOPSIS

  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_registry_from_db(-host => 'ensembldb.ensembl.org', -user => 'anonymous');

  my $ia = $registry->get_adaptor('human', 'variation', 'individual');
  my $sa = $registry->get_adaptor('human', 'variation', 'sample');
  my $pa = $registry->get_adaptor('human', 'variation', 'population');


=head1 DESCRIPTION

This adaptor provides database connectivity for Individual objects.
Individuals may be retrieved from the ensembl variation database by
several means using this module.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::SampleAdaptor;

use Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Variation::Sample;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use DBI qw(:sql_types);
our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor');

sub store {
	my ($self, $sample) = @_;
	
	my $dbh = $self->dbc->db_handle;

  my ($individual, $individual_id); 
  $individual = $sample->individual;
  if ($individual) {
    $individual_id = $individual->dbID;   
  } else {
    $individual_id = $sample->{individual_id};
  }
  if (!$individual_id) {
    my $ia = $self->adaptor->db->get_IndividualAdaptor;
    $individual = Bio::EnsEMBL::Variation::Individual->new(
      name            => $sample->name,
      adaptor         => $ia,
      type_individual => 'outbred',
    );
    
    $individual = $individual->store($individual); 
    $individual_id = $individual->dbID();
    $sample->individual($individual);
  } 
	
	# add entry to sample table
	my $sth = $dbh->prepare(q{
		INSERT INTO sample (
      individual_id,
      name,
      description,
      study_id,
      display,
      has_coverage
		) VALUES (?,?,?,?,?,?)
	});
	$sth->execute(
    $individual_id,
		$sample->name,
    $sample->description,
    $sample->study ? $sample->study->dbID : undef,
    $sample->display,
    $sample->has_coverage
	);
	$sth->finish;
	my $dbID = $dbh->last_insert_id(undef, undef, 'sample', 'sample_id');
	$individual->{dbID}    = $dbID;
	$individual->{adaptor} = $self;

	# store individual/population relationships
	$sth = $dbh->prepare(q{
		INSERT INTO sample_population (
			sample_id,
			population_id
		) VALUES (?,?)
	});
	
	foreach my $population (@{$sample->{populations}}) {
		next unless defined($population->dbID);
		$sth->execute(
			$sample->dbID,
			$population->dbID
		);
	}
	
	$sth->finish;
}

=head2 fetch_by_synonym

    Arg [1]              : String $sample_synonym
    Arg [2]              : String $source (optional)
    Example              : my $ind = $sample_adaptor->fetch_by_synonym($sample_synonym, $source);
    Description          : Retrieves sample for the synonym given in the source. If no source is provided, retrieves all the synonyms
    Returntype           : listref of Bio::EnsEMBL::Variation::Sample objects
    Exceptions           : none
    Caller               : general
    Status               : At Risk

=cut

sub fetch_by_synonym {
  my $self = shift;
  my $synonym_name = shift;
  my $source = shift;
  my ($samples, $sample_ids, $sample_id);
  my $sql;
  if (defined $source) {
    $sql = qq{ SELECT syn.sample_id 
               FROM sample_synonym syn, source s
               WHERE syn.name = ? and syn.source_id = s.source_id AND s.name = "$source"};
  }
  else {
    $sql = qq{ SELECT sample_id FROM sample_synonym WHERE name = ?};
  }
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $synonym_name, SQL_VARCHAR);
  $sth->execute();    
  $sth->bind_columns(\$sample_id);
  while ($sth->fetch()) {
    push @{$sample_ids}, $sample_id;
  }

  foreach my $sample_id (@{$sample_ids}) {
    my $sample = $self->fetch_by_dbID($sample_id);
    push @{$samples}, $sample;
  }
  return $samples;
}

=head2 fetch_synonyms

    Arg [1]              : $sample_id
    Arg [2] (optional)   : $source
    Example              : my $dbSNP_synonyms = $sample_adaptor->fetch_synonyms($sample_id, $dbSNP);
                           my $all_synonyms = $sample_adaptor->fetch_synonyms($sample_id);
    Description: Retrieves synonyms for the source provided. Otherwise, return all the synonyms for the sample_id
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
  if (defined $source){
    $sql = qq{SELECT syn.name FROM sample_synonym syn, source s WHERE syn.sample_id = ? AND syn.source_id = s.source_id AND s.name = "$source"}
  } else{
    $sql = qq{SELECT name FROM sample_synonym WHERE sample_id = ?};
  }
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $dbID, SQL_INTEGER);
  $sth->execute();
  $sth->bind_columns(\$synonym);
  while ($sth->fetch){
    push @{$synonyms}, $synonym;
  }
  return $synonyms;
}

=head2 fetch_all_by_name

  Arg [1]    : string $name 
               The name of the samples to retrieve.
  Example    : my @samples = @{$sample_adaptor->fetch_all_by_name('CEPH1332.05')};
  Description: Retrieves all samples with a specified name. Sample
               names may be non-unique which is why this method returns a
               reference to a list.
  Returntype : listref of Bio::EnsEMBL::Variation::Sample objects
  Exceptions : throw if no argument passed
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_name {
  my $self = shift;
  my $name = shift;

  defined($name) || throw("Name argument expected.");

  my $sth = $self->prepare(q{
    SELECT sample_id, individual_id, name, description, study_id, display, has_coverage
    FROM   sample
    WHERE  name = ?;});

  $sth->bind_param(1, $name, SQL_VARCHAR);
  $sth->execute();
  my $samples =  $self->_objs_from_sth($sth);
  $sth->finish();
  return $samples;
}

=head2 fetch_all_by_name_list

  Arg [1]    : listref of samples names
  Example    : $samples = $sample_adaptor->fetch_all_by_name_list(["NA12347", "NA12348"]);
  Description: Retrieves a listref of Sample objects via a list of names
  Returntype : listref of Bio::EnsEMBL::Variation::Sample objects
  Exceptions : throw if list argument is not defined
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_name_list {
  my $self = shift;
  my $list = shift;

  if (!defined($list) || ref($list) ne 'ARRAY') {
    throw("list reference argument is required");
  }
  
  return [] unless scalar @$list >= 1;
  
  my $id_str = (@$list > 1)  ? " IN (". join(',', map {'"'.$_.'"'} @$list).")" : ' = \''.$list->[0].'\'';
  
  return $self->generic_fetch("name" . $id_str);
}

=head2 fetch_all_by_Population

  Arg [1]    : Bio::EnsEMBL::Variation::Population $population
  Example    : my $population = $population_adaptor->fetch_by_name('PACIFIC');
               foreach my $sample (@{$sample_adaptor->fetch_all_by_Population($population)}) {
                 print $sample->name(), "\n";
               }
  Description: Retrieves all samples from a specified population 
  Returntype : listref of Bio::EnsEMBL::Variation::Sample objects
  Exceptions : throw if incorrect argument is passed
               warning if provided Population does not have an dbID
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Population {
  my $self = shift;
  my $pop = shift;

  if (!ref($pop) || !$pop->isa('Bio::EnsEMBL::Variation::Population')) {
    throw("Bio::EnsEMBL::Variation::Population arg expected");
  }

  if (!$pop->dbID()) {
    warning("Population does not have dbID, cannot retrieve Samples");
    return [];
  }

  my $sth = $self->prepare(q{
    SELECT s.sample_id, s.individual_id, s.name, s.description, s.study_id, s.display, s.has_coverage 
    FROM   sample s, sample_population sp
    WHERE  s.sample_id = sp.sample_id
    AND    sp.population_id = ?});

  $sth->bind_param(1, $pop->dbID, SQL_INTEGER);
  $sth->execute();
  my $samples = $self->_objs_from_sth($sth);
  $sth->finish();
  return $samples;
}

sub _get_name_by_dbID {
  my $self = shift;
  my $dbID = shift;
  
  my $sth = $self->dbc->prepare(qq{
    SELECT name
    FROM sample
    WHERE sample_id = ?
  });
  $sth->bind_param(1, $dbID, SQL_INTEGER);
  $sth->execute();
  
  my $name;
  $sth->bind_columns(\$name);
  $sth->fetch;
  $sth->finish;
  
  return $name;
}

=head2 fetch_all_strains
    Args       : none
    Example    : my $strains = $sample_adaptor->fetch_all_strains();
    Description: Retrieves Samples that should be considered as strain (fully inbred).
    Returntype : listref of Bio::EnsEMBL::Variation::Sample objetcs
    Exceptions : none
    Caller     : Bio:EnsEMBL:Variation::Sample
    Status     : At Risk
=cut

sub fetch_all_strains {
    my $self = shift;
    return $self->generic_fetch("s.display in ('REFERENCE','DEFAULT','DISPLAYABLE')");
}


=head2 get_display_strains

    Args       : none
    Example    : my $strains = $sample_adaptor->get_display_strains();
    Description: Retrieves strain_names that are going to be displayed in the web (reference + default + others)
    Returntype : list of strings
    Exceptions : none
    Caller     : web
    Status     : At Risk

=cut

sub get_display_strains {
    my $self = shift;
    my @strain_names;
    my $name;
    #first, get the reference strain
    $name = $self->get_reference_strain_name();
    push @strain_names, $name;
    #then, get the default ones
    my $default_strains = $self->get_default_strains();
    push @strain_names, @{$default_strains};
    #and finally, get the others
    my $sth = $self->prepare(qq{SELECT name FROM sample WHERE display = ?});
    $sth->bind_param(1, 'DISPLAYABLE');
    $sth->execute;
    $sth->bind_columns(\$name);
    while ($sth->fetch()){
      push @strain_names, $name;
    }
    $sth->finish;
    return \@strain_names;
}

=head2 get_default_strains

    Args       : none
    Example    : my $strains = $sample_adaptor->get_default_strains();
    Description: Retrieves strain_names that are defined as default in the database(mainly, for web purposes)
    Returntype : list of strings
    Exceptions : none
    Caller     : web
    Status     : At Risk

=cut

sub get_default_strains {
    my $self = shift;
    my @strain_names;
    my $name;
    my $sth = $self->prepare(qq{SELECT name FROM sample WHERE display = ?});
    $sth->bind_param(1, 'DEFAULT');
    $sth->execute;
    $sth->bind_columns(\$name);
    while ($sth->fetch()){
        push @strain_names, $name;
    }
    $sth->finish;
    return \@strain_names;
}

=head2 get_reference_strain_name

    Args       : none
    Example    : my $reference_strain = $sample_adaptor->get_reference_strain_name();
    Description: Retrieves the reference strain_name that is defined as default in the database(mainly, for web purposes)
    Returntype : string
    Exceptions : none
    Caller     : web
    Status     : At Risk

=cut

sub get_reference_strain_name {
    my $self = shift;
    my $name;
    my $sth = $self->prepare(qq{SELECT name FROM sample WHERE display = ?});
    $sth->bind_param(1, 'REFERENCE');
    $sth->execute;
    $sth->bind_columns(\$name);
    $sth->fetch();
    $sth->finish;
    return $name;
}

#
# private method, constructs Samples from an executed statement handle
# ordering of columns must be consistant
#
sub _objs_from_sth {
  my $self = shift;
  my $sth = shift;

  my ($dbID, $individual_id, $name, $desc, $study_id, $display_flag, $has_coverage);
  $sth->bind_columns(\$dbID, \$individual_id, \$name, \$desc, \$study_id, \$display_flag, \$has_coverage);

  my %seen;
  my @samples;

  while ($sth->fetch()) {

    my $sample = $seen{$dbID} ||= Bio::EnsEMBL::Variation::Sample->new(
      -dbID          => $dbID,
      -adaptor       => $self,
      -individual_id => $individual_id,
      -name          => $name,
      -description   => $desc,
      -study_id      => $study_id,
      -display       => $display_flag,
      -has_coverage  => $has_coverage,
    );

    $seen{$dbID} = $sample;
    push @samples, $sample;
  }

  return \@samples;
}

sub _tables {
  return (['sample', 's']);
}

sub _columns {
  return qw(s.sample_id s.individual_id s.name s.description s.study_id s.display s.has_coverage);
}

1;
