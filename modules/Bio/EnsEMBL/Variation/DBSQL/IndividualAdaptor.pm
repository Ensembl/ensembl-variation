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
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::IndividualAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::DBSQL::IndividualAdaptor

=head1 SYNOPSIS

  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_registry_from_db(-host => 'ensembldb.ensembl.org', -user => 'anonymous');

  my $ia = $registry->get_adaptor('human', 'variation', 'individual');
  my $pa = $registry->get_adaptor('human', 'variation', 'population');

  # Get all Individuals with a particular name
  foreach my $individual (@{ $ia->fetch_all_by_name('CEPH1362.01') }) {
    print $individual->name(), "\n";
  }

  # get all Individuals from a Population
  my $population = $pa->fetch_by_name('THOWARDEMORY:Coriell');
  foreach my $individual (@{ $ia->fetch_all_by_Population($population) }) {
    print $individual->name(), "\n";
  }

  # get all children of an Individual
  my $individuals = $ia->fetch_all_by_name('CEPH1362.01');
  my $individual  = $individuals->[0];
  foreach my $child (@{ $ia->fetch_all_by_parent_Individual($individual) }) {
    print $child->name(), " is a child of ", $individual->name(), "\n";
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
	my ($self, $individual) = @_;
	
	my $dbh = $self->dbc->db_handle;
	
	# retrieve individual type ID - default to 3 (outbred)
	my $individual_type_id = 3;
	
	if (defined($individual->type_individual)) {
		my $sth = $dbh->prepare(q{
			SELECT individual_type_id
			FROM individual_type
			WHERE name = ?
		});
		$sth->execute($individual->type_individual);
		$sth->bind_columns(\$individual_type_id);
		$sth->fetch();
		$sth->finish();
	}
	
	# add entry to individual table
	my $sth = $dbh->prepare(q{
		INSERT INTO individual (
      name,
      description,
      gender,
      father_individual_id,
      mother_individual_id,
      individual_type_id
		) VALUES (?,?,?,?,?,?)
	});
	$sth->execute(
		$individual->name,
    $individual->description,
    $individual->gender || 'Unknown',
		$individual->father_Individual ? $individual->father_Individual->dbID : undef,
		$individual->mother_Individual ? $individual->mother_Individual->dbID : undef,
		$individual_type_id,
	);
	$sth->finish;
	my $dbID = $dbh->last_insert_id(undef, undef, 'individual', 'individual_id');
	$individual->{dbID}    = $dbID;
	$individual->{adaptor} = $self;
  return $individual;
}

=head2 fetch_individual_by_synonym

    Arg [1]              : String $individual_synonym
    Arg [2]              : String $source (optional)
    Example              : my $ind = $ind_adaptor->fetch_individual_by_synonym($individual_synonym, $source);
    Description          : Retrieves individual for the synonym given in the source. If no source is provided, retrieves all the synonyms
    Returntype           : listref of Bio::EnsEMBL::Variation::Individual objects
    Exceptions           : none
    Caller               : general
    Status               : At Risk

=cut

sub fetch_individual_by_synonym {
    my $self = shift;
    my $synonym_name = shift;
    my $source = shift;
    my ($individuals, $individual_ids, $individual_id);
    my $sql;
    if (defined $source){
	    $sql = qq{
            SELECT isyn.individual_id 
            FROM individual_synonym isyn, source s
            WHERE isyn.name = ? and isyn.source_id = s.source_id AND s.name = "$source"};
    }
    else{
	    $sql = qq{
            SELECT individual_id 
            FROM individual_synonym WHERE name = ?};
    }
    my $sth = $self->prepare($sql);
    $sth->bind_param(1, $synonym_name, SQL_VARCHAR);
    $sth->execute();    
    $sth->bind_columns(\$individual_id);
    while ($sth->fetch()){
	    push @{$individual_ids}, $individual_id;
    }

    foreach my $individual_id (@{$individual_ids}){
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

sub fetch_synonyms{
    my $self = shift;
    my $dbID = shift;
    my $source = shift;
    my $synonym;
    my $synonyms;

    my $sql;
    if (defined $source){
	    $sql = qq{SELECT isyn.name FROM individual_synonym isyn, source s WHERE isyn.individual_id = ? AND isyn.source_id = s.source_id AND s.name = "$source"}
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
    SELECT i.individual_id, i.name, i.description, i.gender, i.father_individual_id, i.mother_individual_id, it.name, it.description
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
        warning("Population does not have dbID, cannot retrieve Individuals");
        return [];
    }

    my $sth = $self->prepare(q{
        SELECT i.individual_id, i.name, i.description, i.gender, i.father_individual_id, i.mother_individual_id, it.name, it.description 
        FROM   individual i, individual_type it, sample s, sample_population sp
        WHERE  sp.sample_id = s.sample_id
        AND    s.individual_id = i.individual_id   
        AND    i.individual_type_id = it.individual_type_id
        AND    sp.population_id = ?});

    $sth->bind_param(1, $pop->dbID,SQL_INTEGER);
    $sth->execute();
    my $individuals = $self->_objs_from_sth($sth);
    $sth->finish();
    return $individuals;
}


=head2 fetch_all_by_parent_Individual

  Arg [1]    : Bio::EnsEMBL::Variation::Individual
  Example    : my @children = @{$ind_adaptor->fetch_all_by_parent_Individual($ind)};
  Description: Retrieves all individuals which are children of a provided
               parent individual. This function operates under the assumptions
               that Male individuals can only be fathers, Female individuals
               can only be mothers and Unknown individuals can only be one
               or the other - not both.
  Returntype : listref of Bio::EnsEMBL::Variation::Individual objects
  Exceptions : throw if incorrect argument passed
               warning if provided individual has no dbID 
  Caller     : general, Individual::get_all_child_Individuals
  Status     : Stable

=cut

sub fetch_all_by_parent_Individual {
    my $self = shift;
    my $parent  = shift;

    if (!ref($parent) || !$parent->isa('Bio::EnsEMBL::Variation::Individual')) {
        throw("Bio::EnsEMBL::Variation::Individual argument expected");
    }

    if (!defined($parent->dbID())) {
        warning("Cannot fetch child Individuals for parent without dbID");
        return [];
    }
    my $gender = $parent->gender() || '';

    my $father_sql = q{
        SELECT i.individual_id, i.name, i.description, i.gender, i.father_individual_id, i.mother_individual_id, it.name, it.description 
        FROM   individual i, individual_type it
        WHERE  i.father_individual_id = ?
        AND    i.individual_type_id = it.individual_type_id;};
    my $mother_sql = q{
        SELECT i.individual_id, i.name, i.description, i.gender, i.father_individual_id, i.mother_individual_id, it.name, it.description
        FROM   individual i, individual_type it
        WHERE  i.mother_individual_id = ?
        AND    i.individual_type_id = it.individual_type_id;};

    if ($gender eq 'Male') {
        my $sth = $self->prepare($father_sql);
        $sth->bind_param(1,$parent->dbID,SQL_INTEGER);
        $sth->execute();
        my $individuals = $self->_objs_from_sth($sth);
        $sth->finish();
        return $individuals;
    } elsif ($gender eq 'Female') {
        my $sth = $self->prepare($mother_sql);
        $sth->bind_param(1,$parent->dbID,SQL_INTEGER);
        $sth->execute();
        my $individuals = $self->_objs_from_sth($sth);
        $sth->finish();
        return $individuals;
    } else { # unknown gender
        my $sth = $self->prepare($mother_sql);
        $sth->bind_param(1,$parent->dbID,SQL_INTEGER);
        $sth->execute();
        my $individuals = $self->_objs_from_sth($sth);
        $sth->finish();
        # if this parent was a mother, finish now and return results
        return $individuals if(@$individuals);
        # otherwise assume was a father (or nothing)
        $sth = $self->prepare($father_sql);
        $sth->bind_param(1,$parent->dbID,SQL_INTEGER);
        $sth->execute();
        $individuals = $self->_objs_from_sth($sth);
        $sth->finish();
        return $individuals;
    }
    return [];
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
# private method, constructs Individuals from an executed statement handle
# ordering of columns must be consistant
#
sub _objs_from_sth {
    my $self = shift;
    my $sth = shift;
    
    my ($dbID, $name, $desc, $gender, $father_id, $mother_id, $it_name, $it_desc);
    $sth->bind_columns(\$dbID, \$name, \$desc, \$gender, \$father_id, \$mother_id, \$it_name, \$it_desc);

    my %seen;
    my %wanted_fathers;
    my %wanted_mothers;
    my @inds;

    while ($sth->fetch()) {
        # get objects for mother and father if they were already constructed
        # otherwise may have to be lazy-loaded later
        my $father;
        if (defined($father_id)) {
            $father = $seen{$father_id};
            if (!$father) {
                $wanted_fathers{$dbID} ||= [];
                push @{$wanted_fathers{$father_id}}, $dbID;
            }
        }
        my $mother;
        if (defined($mother_id)) {
            $mother = $seen{$mother_id};
            if (!$mother) {
                $wanted_mothers{$mother_id} ||= [];
                push @{$wanted_mothers{$mother_id}}, $dbID;
            }
        }

        my $ind = $seen{$dbID} ||= Bio::EnsEMBL::Variation::Individual->new(
            -dbID        => $dbID,
            -adaptor     => $self,
            -description => $desc,
            -gender      => $gender,
            -name        => $name,
            -father_individual => $father,
            -mother_individual => $mother,
            -father_individual_id => $father_id,
            -mother_individual_id => $mother_id,
            -type_individual => $it_name,
            -type_description => $it_desc);

        $seen{$dbID} = $ind;
        push @inds, $ind;
    }

    # load any of the 'wanted' parent individuals that we did not have at the
    # of creation, but which we have now

    foreach my $wanted_id (keys %wanted_fathers) {
        if ($seen{$wanted_id}) {
            # add father to every child that wanted it
            foreach my $ind_id (@{$wanted_fathers{$wanted_id}}) {
                $seen{$ind_id}->father_Individual($seen{$wanted_id});
            }
        }
    }
    foreach my $wanted_id (keys %wanted_mothers) {
        if ($seen{$wanted_id}) {
            # add mother to every child that wanted it
            foreach my $ind_id (@{$wanted_mothers{$wanted_id}}) {
                $seen{$ind_id}->mother_Individual($seen{$wanted_id});
            }
        }
    }
    return \@inds;
}

sub _tables {
    return (['individual','i'],
		['individual_type','it'])}

sub _columns {
    return qw(i.individual_id i.name i.description i.gender i.father_individual_id i.mother_individual_id it.name it.description);
}

sub _default_where_clause {
    return 'i.individual_type_id = it.individual_type_id';
}

1;
