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

  # Get all Samples with a particular name
  foreach my $sample (@{ $sa->fetch_all_by_name('CEPH1362.01') }) {
    print $sample->name(), "\n";
  }

=head1 DESCRIPTION

This adaptor provides database connectivity for Individual objects.
Individuals may be retrieved from the ensembl variation database by
several means using this module.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::IndividualAdaptor;

use Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Variation::Individual;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use DBI qw(:sql_types);
our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor');

sub store {
	my ($self, $sample) = @_;
	
	my $dbh = $self->dbc->db_handle;

  my $individiual  = $sample->individual;
  if (!$individual) {
    my $ia = $self->adaptor->db->get_IndividualAdaptor;
    $individual = Bio::EnsEMBL::Variation::Individual->new(
      name            => $sample->name,
      adaptor         => $ia,
      type_individual => 'outbred',
    );
    
    $individual->store($individual); 
    $sample->individual($individual);
  } 
	
	# add entry to sample table
	my $sth = $dbh->prepare(q{
		INSERT INTO sample (
      individual_id,
      name,
      description,
      display,
      has_coverage
		) VALUES (?,?,?,?,?,?,?,?)
	});
	$sth->execute(
    $sample->individual->dbID,
		$sample->name,
    $sample->description,
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
    $sql = qq{ SELECT is.individual_id 
               FROM sample_synonym syn, source s
               WHERE syn.name = ? and syn.source_id = s.source_id AND s.name = "$source"};
  }
  else {
    $sql = qq{ SELECT sample_id 
               FROM sample_synonym WHERE name = ?};
  }
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $synonym_name, SQL_VARCHAR);
  $sth->execute();    
  $sth->bind_columns(\$sample_id);
  while ($sth->fetch()) {
    push @{$sample_ids}, $sample_id;
  }

  foreach my $sample_id (@{$sample_ids}) {
    my $individual = $self->fetch_by_dbID($individual_id);
    push @{$individuals}, $individual;
  }
  return $individuals;
}

=head2 fetch_synonyms

    Arg [1]              : $individual_id
    Arg [2] (optional)   : $source
    Example              : my $dbSNP_synonyms = $ind_adaptor->fetch_synonyms($individual_id,$dbSNP);
                           my $all_synonyms = $ind_adaptor->fetch_synonyms($individual_id);
    Description: Retrieves synonyms for the source provided. Otherwise, return all the synonyms for the individual_id
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
	    $sql = qq{SELECT is.name FROM individual_synonym is, source s WHERE is.individual_id = ? AND is.source_id = s.source_id AND s.name = "$source"}
    } else{
	    $sql = qq{SELECT name FROM individual_synonym WHERE individual_id = ?};
    }
    my $sth = $self->prepare($sql);
    $sth->bind_param(1,$dbID,SQL_INTEGER);
    $sth->execute();
    $sth->bind_columns(\$synonym);
    while ($sth->fetch){
	    push @{$synonyms}, $synonym;
    }
    return $synonyms;
}

=head2 fetch_all_by_name

  Arg [1]    : string $name 
               The name of the individuals to retrieve.
  Example    : my @inds = @{$ind_adaptor->fetch_all_by_name('CEPH1332.05')};
  Description: Retrieves all individuals with a specified name.  Individual
               names may be non-unique which is why this method returns a
               reference to a list.
  Returntype : listref of Bio::EnsEMBL::Variation::Individual objects
  Exceptions : throw if no argument passed
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_name {
  my $self = shift;
  my $name = shift;

  defined($name) || throw("Name argument expected.");

  my $sth = $self->prepare(q{
    SELECT i.individual_id, i.name, i.description, i.gender, i.father_individual_id, i.mother_individual_id, it.name, it.description, i.display, i.has_coverage
    FROM   individual i, individual_type it
    WHERE  i.name = ?
    AND    it.individual_type_id = i.individual_type_id;});

  $sth->bind_param(1, $name, SQL_VARCHAR);
  $sth->execute();
  my $individuals =  $self->_objs_from_sth($sth);
  $sth->finish();
  return $individuals;
}

=head2 fetch_all_by_name_list

  Arg [1]    : listref of individual names
  Example    : $inds = $ind_adaptor->fetch_all_by_name_list(["NA12347", "NA12348"]);
  Description: Retrieves a listref of individual objects via a list of names
  Returntype : listref of Bio::EnsEMBL::Variation::Individual objects
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
  
  my $id_str = (@$list > 1)  ? " IN (".join(',', map {'"'.$_.'"'} @$list).")"   :   ' = \''.$list->[0].'\'';
  
  return $self->generic_fetch("i.name ".$id_str);
}

=head2 fetch_all_by_Population

  Arg [1]    : Bio::EnsEMBL::Variation::Population $pop
  Example    : my $pop = $pop_adaptor->fetch_by_name('PACIFIC');
               foreach my $ind (@{$ind_adaptor->fetch_all_by_Population($pop)}) {
                 print $ind->name(), "\n";
               }
  Description: Retrieves all individuals from a specified population 
  Returntype : listref of Bio::EnsEMBL::Variation::Individual objects
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
        SELECT i.individual_id, i.name, i.description, i.gender, i.father_individual_id, i.mother_individual_id, it.name, it.description, i.display, i.has_coverage 
        FROM   sample s, sample_population sp
        WHERE  s.sample_id = sp.individual_id
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
    FROM individual
    WHERE individual_id = ?
  });
  $sth->bind_param(1,$dbID,SQL_INTEGER);
  $sth->execute();
  
  my $name;
  $sth->bind_columns(\$name);
  $sth->fetch;
  $sth->finish;
  
  return $name;
}

sub _get_sample_name_by_dbID {
    my $self = shift;
    warn('The use od this method is deprecated. Use _get_name_by_dbID instead.')
}

#
# private method, constructs Samples from an executed statement handle
# ordering of columns must be consistant
#
sub _objs_from_sth {
  my $self = shift;
  my $sth = shift;

  my ($dbID, $name, $desc, $gender, $father_id, $mother_id, $it_name, $it_desc, $display_flag, $has_coverage);
  $sth->bind_columns(\$dbID, \$name, \$desc, \$gender, \$father_id, \$mother_id, \$it_name, \$it_desc, \$display_flag, \$has_coverage);

  my %seen;
  my @samples;

  while ($sth->fetch()) {

    my $sample = $seen{$dbID} ||= Bio::EnsEMBL::Variation::Sample->new(
      -dbID          => $dbID,
      -adaptor       => $self,
      -individual_id => $individual_id,
      -description   => $desc,
      -display       => $display_flag,
      -has_coverage  => $has_coverage,
      -name          => $name,
    );

    $seen{$dbID} = $sample;
    push @samples, $sample;
  }

  return \@samples;
}

sub _tables {
    return (['individual','i'],
		['individual_type','it'])}

sub _columns {
    return qw(i.individual_id i.name i.description i.gender i.father_individual_id i.mother_individual_id it.name it.description i.display i.has_coverage);
}

sub _default_where_clause {
    return 'i.individual_type_id = it.individual_type_id';
}

1;
