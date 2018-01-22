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
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::SampleAdaptor

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
    if (!$individual_id) {
      my $ia = $individual->adaptor;
      $individual = $ia->store($individual); 
      $individual_id = $individual->dbID;
    }
  } else {
    $individual_id = $sample->{individual_id};
    if (!$individual_id) {
      my $ia = $sample->adaptor->db->get_IndividualAdaptor;
      $individual = Bio::EnsEMBL::Variation::Individual->new(
        -name            => $sample->name,
        -adaptor         => $ia,
        -type_individual => 'outbred',
      );
      $individual = $ia->store($individual); 
      $individual_id = $individual->dbID();
      $sample->individual($individual);
    }
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
	$sample->{dbID}    = $dbID;
	$sample->{adaptor} = $self;

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

	# store synonyms
	$sth = $dbh->prepare(q{
		INSERT INTO sample_synonym (
			sample_id,
      source_id,
      name
		) VALUES (?, ?, ?)
	});

  my $source_adaptor = $self->db->get_SourceAdaptor();

  foreach my $source_name (keys %{$sample->{'synonyms'}}) {
    my $source = $self->{'sources'}->{$source_name};
    if (!$source) {
      $source = $source_adaptor->fetch_by_name($source_name);
      $self->{'sources'}->{$source_name} = $source;
    }   
    my $source_id = $source->dbID;
    foreach my $synonym_name (keys %{$sample->{'synonyms'}->{$source_name}}) {
		  $sth->execute(
			  $sample->dbID,
			  $source_id,
        $synonym_name
		  );
    }
  }

	$sth->finish;
}

sub update {
  my ($self, $sample) = @_;

  my $dbh = $self->dbc->db_handle;

	# store synonyms
  my $sth = $dbh->prepare(q{
		INSERT INTO sample_synonym (
			sample_id,
      source_id,
      name
		) VALUES (?, ?, ?)
	});

  # update synonyms
  my $source_adaptor = $self->db->get_SourceAdaptor();

  foreach my $source_name (keys %{$sample->{'synonyms'}}) {
    my $source = $self->{'sources'}->{$source_name};
    if (!$source) {
      $source = $source_adaptor->fetch_by_name($source_name);
      $self->{'sources'}->{$source_name} = $source;
    }   
    my $source_id = $source->dbID;
    foreach my $synonym_name (keys %{$sample->{'synonyms'}->{$source_name}}) {
      my $already_in_db = $self->fetch_synonyms($sample->dbID, $source_name);
      my ($synonym) = grep {$_ eq $synonym_name} @$already_in_db;
      if (!$synonym) {
		    $sth->execute(
			    $sample->dbID,
			    $source_id,
          $synonym_name
		    );
      }
    }
  } 
  $sth->finish;

  # update sample_population
	$sth = $dbh->prepare(q{
		INSERT INTO sample_population (
			sample_id,
			population_id
		) VALUES (?,?)
	});
	
	foreach my $population (@{$sample->{populations}}) {
    my $samples = $self->fetch_all_by_Population($population);
    my ($thisSample) = grep {$_->dbID eq $sample->dbID} @$samples;
    if (!$thisSample) {
		  $sth->execute(
			  $sample->dbID,
			  $population->dbID
		  );
    }
	}
  $sth->finish;

}

=head2 fetch_by_synonym

  Arg [1]     : String $sample_synonym
  Arg [2]     : String $source_name (optional)
  Example     : my $sample = $sample_adaptor->fetch_by_synonym($sample_synonym, $source_name);
  Description : Retrieves sample for the synonym and optional source name given. If no source is provided, retrieves all the synonyms
  Returntype  : Listref of Bio::EnsEMBL::Variation::Sample objects
  Exceptions  : Throws on missing synonym name and on wrong source name (Source name couldn't be found in database)
  Caller      : General
  Status      : Stable

=cut

sub fetch_by_synonym {
  my $self = shift;
  my $synonym_name = shift;
  my $source_name = shift;

  if (!defined($synonym_name)) {
    throw("Synonym name argument is required");
  }
  my $constraint = '';
  if ($synonym_name && $source_name) {
    my $source_adaptor = $self->db->get_SourceAdaptor();
    my $source = $source_adaptor->fetch_by_name($source_name);
    throw("Source name could not be found in the database") if (!$source);
    $constraint = qq{ss.name = ? AND ss.source_id = ?};
    $self->bind_param_generic_fetch($synonym_name, SQL_VARCHAR);
    $self->bind_param_generic_fetch($source->dbID, SQL_INTEGER);
  } else {
    $constraint = qq{ss.name = ?};
    $self->bind_param_generic_fetch($synonym_name, SQL_VARCHAR);
  }

  my $result = $self->generic_fetch($constraint);

  return $result;
}

=head2 fetch_synonyms

  Arg [1]            : $sample_id
  Arg [2] (optional) : $source
  Example            : my $dbSNP_synonyms = $sample_adaptor->fetch_synonyms($sample_id, $dbSNP);
                       my $all_synonyms = $sample_adaptor->fetch_synonyms($sample_id);
  Description        : Retrieves synonyms for the source provided. Otherwise, return all the synonyms for the sample_id
  Returntype         : Listref of strings
  Exceptions         : None
  Caller             : General
  Status             : Stable

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

  Arg [1]    : String $name 
               The name of the samples to retrieve.
  Example    : my @samples = @{$sample_adaptor->fetch_all_by_name('CEPH1332.05')};
  Description: Retrieves all samples with a specified name. Sample
               names may be non-unique which is why this method returns a
               reference to a list.
  Returntype : Listref of Bio::EnsEMBL::Variation::Sample objects
  Exceptions : Throw if no argument passed
  Caller     : General
  Status     : Stable

=cut

sub fetch_all_by_name {
  my $self = shift;
  my $name = shift;
  defined($name) || throw("Name argument expected.");

  my $constraint = qq{s.name = ?};
  $self->bind_param_generic_fetch($name, SQL_VARCHAR);
  my $result = $self->generic_fetch($constraint);
  if (scalar @$result == 0) {
    $result = $self->fetch_by_synonym($name);    
  }
  return $result;
}

=head2 fetch_all_by_name_list
  Arg [1]    : Listref of samples names
  Example    : $samples = $sample_adaptor->fetch_all_by_name_list(["NA12347", "NA12348"]);
  Description: Retrieves a listref of Sample objects via a list of names
  Returntype : Listref of Bio::EnsEMBL::Variation::Sample objects
  Exceptions : Throw if list argument is not defined
  Caller     : General
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
  
  return $self->generic_fetch("s.name" . $id_str . " OR ss.name" . $id_str);
}

=head2 fetch_all_by_Population

  Arg [1]    : Bio::EnsEMBL::Variation::Population $population
  Example    : my $population = $population_adaptor->fetch_by_name('PACIFIC');
               foreach my $sample (@{$sample_adaptor->fetch_all_by_Population($population)}) {
                 print $sample->name(), "\n";
               }
  Description: Retrieves all samples from a specified population 
  Returntype : Listref of Bio::EnsEMBL::Variation::Sample objects
  Exceptions : Throw if incorrect argument is passed
               Warning if provided Population does not have an dbID
  Caller     : General
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

  $self->{'_constrain_population'} = 1;

    # Add a constraint on the subsnp_id
  my $constraint = qq{sp.population_id = ?};

    # Bind the parameters
  $self->bind_param_generic_fetch($pop->dbID, SQL_INTEGER);

    # Get the results from generic fetch method
  my $result = $self->generic_fetch($constraint);

  delete($self->{'_constrain_population'});

  return $result;
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
    # temporary fix for release 82
    if ($name eq 'refGRCh37') {
      $name = 'Homo_sapiens';
    }
    return $name;
}

=head2 fetch_reference_strain

    Args       : none
    Example    : my $reference_strain = $sample_adaptor->fetch_reference_strain;
    Description: Retrieves the reference strain
    Returntype : Bio::EnsEMBL::Variation::Sample 
    Exceptions : none
    Caller     : general
    Status     : stable

=cut

sub fetch_reference_strain {
  my $self = shift;
  my $strains = $self->generic_fetch("s.display = 'REFERENCE'");
  return undef unless (scalar(@{$strains}));
  return $strains->[0];
}

sub _tables {
  my $self = shift;
  
  my @tables = (
    ['sample', 's'],
    ['sample_synonym', 'ss'],
    ['source', 's1']
  );

  push(@tables, ['sample_population', 'sp']) if ($self->{'_constrain_population'});

  return @tables;
}

sub _columns {
  my $self = shift;
  my @cols = (
    's.sample_id',
    's.individual_id',
    's.name',
    's.description',
    's.study_id',
    's.display',
    's.has_coverage',
    'ss.name AS synonym_name',
    's1.name AS synonym_source_name',
  );
  return @cols;
}

sub _left_join {
  my $self = shift;

  my @left_join = (
    ['sample_synonym', 's.sample_id = ss.sample_id'],
    ['source s1', 'ss.source_id = s1.source_id']
  );

  push (@left_join, ['sample_population', 's.sample_id = sp.sample_id']) if ($self->{'_constrain_population'});

  return @left_join;
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

  # Return the created objects
  return \@objs;
}

sub _obj_from_row {
  my ($self, $row) = @_;

  return undef unless $row->{sample_id};

  # If the sample for this sample_id hasn't already been created, do that
  my $obj = $self->{_temp_objs}{$row->{sample_id}};

  unless (defined($obj)) {
    # Create the variation object
    $obj = Bio::EnsEMBL::Variation::Sample->new(
      -adaptor => $self,
      -dbID => $row->{sample_id},
      -individual_id => $row->{individual_id},            
      -name => $row->{name},
      -description => $row->{description},
      -study_id => $row->{study_id},
      -display => $row->{display},
      -has_coverage => $row->{has_coverage}
    );
    $self->{_temp_objs}{$row->{sample_id}} = $obj;
  }

  # Add a synonym if available
  if (defined($row->{synonym_source_name}) && defined($row->{synonym_name})) {
    $obj->add_synonym($row->{synonym_source_name}, $row->{synonym_name});
  }

}

1;
