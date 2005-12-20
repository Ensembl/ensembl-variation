#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::PopulationAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::PopulationAdaptor

=head1 SYNOPSIS

  $db = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(...);

  $pa = $db->get_PopulationAdaptor();

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

=head1 AUTHOR - Graham McVicker

=head1 CONTACT

Post questions to the Ensembl development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::PopulationAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Variation::Population;

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');



=head2 fetch_by_dbID

  Arg [1]    : int $dbID
  Example    : $pop = $pop_adaptor->fetch_by_dbID(5543);
  Description: Retrieves a Population object via its internal identifier.
               If no such population exists undef is returned.
  Returntype : Bio::EnsEMBL::Variation::Population
  Exceptions : throw if dbID arg is not defined
  Caller     : general, IndividualAdaptor

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  throw('dbID argument expected') if(!defined($dbID));

  my $sth = $self->prepare(q{SELECT s.sample_id, s.name, s.size, s.description, p.is_strain
                             FROM   population p, sample s
                             WHERE  p.sample_id = s.sample_id
			     AND    p.sample_id = ?});
  $sth->execute($dbID);
  my $result = $self->_objs_from_sth($sth);
  $sth->finish();

  return undef if(!@$result);

  return $result->[0];
}



=head2 fetch_by_name

  Arg [1]    : string $name
  Example    : $pop = $pop_adaptor->fetch_by_name('NUSPAE:Singapore_HDL');
  Description: Retrieves a population object via its name
  Returntype : Bio::EnsEMBL::Variation::Population
  Exceptions : throw if name argument is not defined
  Caller     : general

=cut

sub fetch_by_name {
  my $self = shift;
  my $name = shift;

  throw('name argument expected') if(!defined($name));

  my $sth = $self->prepare(q{SELECT s.sample_id, s.name, s.size, s.description, p.is_strain
                             FROM   population p, sample s
                             WHERE  p.sample_id = s.sample_id 
			     AND    s.name = ?});

  $sth->execute($name);

  my $result = $self->_objs_from_sth($sth);

  $sth->finish();

  return undef if(!@$result);

  return $result->[0];
}


=head2 fetch_all_by_super_Population

  Arg [1]    : Bio::EnsEMBL::Variation::Population $pop
  Example    : foreach $sub_pop (@{$pa->fetch_all_by_super_Population($pop)}) {
                 print $sub_pop->name(), "\n";
               }
  Description: Retrieves all sub populations of a provided population.
  Returntype : Bio::EnsEMBL::Population
  Exceptions : throw on bad argument
  Caller     : general

=cut

sub fetch_all_by_super_Population {
  my $self = shift;
  my $pop  = shift;

  if(!ref($pop) || !$pop->isa('Bio::EnsEMBL::Variation::Population')) {
    throw('Bio::EnsEMBL::Variation::Population argument expected');
  }

  if(!$pop->dbID()) {
    warning("Cannot retrieve sub populations for population without dbID");
    return [];
  }

  my $sth = $self->prepare(q{SELECT s.sample_id, s.name, s.size,
                                    s.description, p.is_strain
                             FROM   population p, sample s, population_structure ps
                             WHERE  s.sample_id = ps.sub_population_sample_id
			     AND    s.sample_id = p.sample_id
                             AND    ps.super_population_sample_id = ?});

  $sth->execute($pop->dbID());

  my $result = $self->_objs_from_sth($sth);

  $sth->finish();

  return $result;
}



=head2 fetch_all_by_sub_Population

  Arg [1]    : Bio::EnsEMBL::Variation::Population $pop
  Example    : foreach $super_pop (@{$pa->fetch_all_by_sub_Population($pop)}) {
                 print $super_pop->name(), "\n";
               }
  Description: Retrieves all super populations for a provided population
  Returntype : Bio::EnsEMBL::Population
  Exceptions : throw on bad argument
  Caller     : general

=cut

sub fetch_all_by_sub_Population {
  my $self = shift;
  my $pop  = shift;

  if(!ref($pop) || !$pop->isa('Bio::EnsEMBL::Variation::Population')) {
    throw('Bio::EnsEMBL::Variation::Population argument expected');
  }

  if(!$pop->dbID()) {
    warning("Cannot retrieve super populations for population without dbID");
    return [];
  }

  my $sth = $self->prepare(q{SELECT s.sample_id, s.name, s.size,
                                    s.description, p.is_strain
                             FROM   population p, sample s, population_structure ps
                             WHERE  s.sample_id = ps.super_population_sample_id
			     AND    s.sample_id = p.sample_id
                             AND    ps.sub_population_sample_id = ?});

  $sth->execute($pop->dbID());

  my $result = $self->_objs_from_sth($sth);

  $sth->finish();

  return $result;
}

=head2 fetch_synonyms

    Arg [1]              : $pop_id
    Arg [2] (optional)   : $source
    Example              : my $dbSNP_synonyms = $pop_adaptor->fetch_synonyms($dbSNP);
                           my $all_synonyms = $pop_adaptor->fetch_synonyms();
    Description: Retrieves synonyms for the source provided. Otherwise, return all the synonyms for the population
    Returntype : list of strings
    Exceptions : none
    Caller     : Bio:EnsEMBL:Variation::Population

=cut

sub fetch_synonyms{
    my $self = shift;
    my $dbID = shift;
    my $source = shift;
    my $population_synonym;
    my $synonyms;

    my $sql;
    if (defined $source){
	$sql = qq{SELECT ss.name FROM sample_synonym ss, source s WHERE ss.sample_id = ? AND ss.source_id = s.source_id AND s.name = "$source"}
    }
    else{
	$sql = qq{SELECT name FROM sample_synonym WHERE sample_id = ?};
    }
    my $sth = $self->prepare($sql);
    $sth->execute($dbID);
    $sth->bind_columns(\$population_synonym);
    while ($sth->fetch){
	push @{$synonyms},$population_synonym;
    }
    return $synonyms;
}

=head2 fetch_population_by_synonym

    Arg [1]              : $population_synonym
    Example              : my $pop = $pop_adaptor->fetch_population_by_synonym($population_synonym,$source);
    Description          : Retrieves population for the synonym given in the source. If no source is provided, retrieves all the synonyms
    Returntype           : list of Bio::EnsEMBL::Variation::Population
    Exceptions           : none
    Caller               : general
d
=cut

sub fetch_population_by_synonym{
    my $self = shift;
    my $synonym_name = shift;
    my $source = shift;
    my $sql;
    my $population;
    my $population_array;

    if (defined $source){
	$sql = qq{SELECT ss.sample_id FROM sample_synonym ss, source s WHERE ss.name = ? and ss.source_id = s.source_id = s.name = "$source"};
    }
    else{
	$sql = qq{SELECT sample_id FROM sample_synonym WHERE name = ?};
    }
    my $population_id;
    my $sth = $self->prepare($sql);
    $sth->execute($synonym_name);
    $sth->bind_columns(\$population_id);
    while ($sth->fetch()){
	$population = $self->fetch_by_dbID($population_id);
	push @{$population_array}, $population;
    }
    return $population_array;
    
}

=head2 fetch_all_strains

    Args       : none
    Example    : my $strains = $pop_adaptor->fetch_all_strains();
    Description: Retrieves populations that should be considered as strain in the specie.
    Returntype : list of Bio::EnsEMBL::Variation::Population
    Exceptions : none
    Caller     : Bio:EnsEMBL:Variation::Population

=cut

sub fetch_all_strains{
    my $self = shift;
    
    return $self->generic_fetch("is_strain = 1");

}


=head2 fetch_all_by_Individual

  Arg [1]     : Bio::EnsEMBL::Variation::Individual $ind
  Example     : my $ind = $ind_adaptor->fetch_by_name('NA12004');
                foreach my $pop (@{$pop_adaptor->fetch_all_by_Individual($ind)}){
		    print $pop->name,"\n";
                }
  Description : Retrieves all populations from a specified individual
  ReturnType  : reference to list of Bio::EnsEMBL::Variation::Population objects
  Exceptions  : throw if incorrect argument is passed
                warning if provided individual does not have a dbID
  Caller      : general

=cut

sub fetch_all_by_Individual{
    my $self = shift;
    my $ind = shift;

    if(!ref($ind) || !$ind->isa('Bio::EnsEMBL::Variation::Individual')) {
	throw("Bio::EnsEMBL::Variation::Individual arg expected");
    }
    
    if(!$ind->dbID()) {
	warning("Individual does not have dbID, cannot retrieve Individuals");
	return [];
  } 

    my $sth = $self->prepare(qq{SELECT s.sample_id, s.name, s.size, s.description, p.is_strain
				FROM population p, sample s, individual_population ip
				WHERE s.sample_id = ip.population_sample_id
				AND   s.sample_id = p.sample_id
                                AND ip.individual_sample_id = ?
			    });
    $sth->execute($ind->dbID());

    my $results = $self->_objs_from_sth($sth);

    $sth->finish();

    return $results;
}


=head2 fetch_tagged_Population

  Arg [1]     : Bio::EnsEMBL::Variation::VariationFeature $vf
  Example     : my $vf = $vf_adaptor->fetch_by_name('rs205621');
                my $populations_tagged = $vf->is_tagged();
                foreach my $pop (@{$vf_adaptor->is_tagged}){
		    print $pop->name," has been tagged using a 0.99 r2 criteria\n";
                }
  Description : Retrieves all populations from a specified variation feature that have been tagged
  ReturnType  : reference to list of Bio::EnsEMBL::Variation::Population objects
  Exceptions  : throw if incorrect argument is passed
                warning if provided variation feature does not have a dbID
  Caller      : general

=cut

sub fetch_tagged_Population{
    my $self = shift;
    my $variation_feature = shift;

    if(!ref($variation_feature) || !$variation_feature->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
	throw("Bio::EnsEMBL::Variation::VariationFeature arg expected");
    }
    
    if(!$variation_feature->dbID()) {
	warning("Variation feature does not have dbID, cannot retrieve tagged populations");
	return [];
  } 

    my $sth = $self->prepare(qq{SELECT s.sample_id, s.name, s.size, s.description, p.is_strain
				FROM population p, sample s, tagged_variation_feature tvf
				WHERE s.sample_id = tvf.sample_id
				AND   s.sample_id = p.sample_id
                                AND tvf.variation_feature_id = ?
			    });
    $sth->execute($variation_feature->dbID());
    my $results = $self->_objs_from_sth($sth);

    $sth->finish();

    return $results;
}

#
# private method, creates population objects from an executed statement handle
# ordering of columns must be consistant
#
sub _objs_from_sth {
  my $self = shift;
  my $sth  = shift;

  my @pops;

  my ($pop_id, $name, $size, $desc, $is_strain);

  $sth->bind_columns(\$pop_id, \$name, \$size, \$desc, \$is_strain);

  while($sth->fetch()) {
    push @pops, Bio::EnsEMBL::Variation::Population->new
      (-dbID => $pop_id,
       -ADAPTOR => $self,
       -NAME => $name,
       -DESCRIPTION => $desc,
       -SIZE => $size,
       -IS_STRAIN => $is_strain);
  }

  return \@pops;
}

sub _tables{return['sample','s'],['population','p']}


sub _columns{
    return qw(s.sample_id s.name s.size s.description p.is_strain
	      );
}

sub _default_where_clause {
  my $self = shift;

  return 'p.sample_id = s.sample_id ';
}

1;
