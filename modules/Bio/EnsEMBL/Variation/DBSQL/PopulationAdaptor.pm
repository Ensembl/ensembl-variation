=head1 LICENSE

 Copyright (c) 1999-2012 The European Bioinformatics Institute and
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

use Bio::EnsEMBL::Variation::DBSQL::SampleAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Scalar qw(wrap_array);

use Bio::EnsEMBL::Variation::Population;

our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::SampleAdaptor');

sub store {
	my ($self, $pop) = @_;
	
	my $dbh = $self->dbc->db_handle;
    
    my $sth = $dbh->prepare(q{
        INSERT INTO sample (
            name,
			size,
			description,
			freqs_from_gts
        ) VALUES (?,?,?,?)
    });
	
	$sth->execute(
		$pop->name,
		$pop->size,
		$pop->description,
		$pop->_freqs_from_gts || 0,
	);
	$sth->finish;
	
	# get the sample_id inserted
	my $dbID = $dbh->last_insert_id(undef, undef, 'sample', 'sample_id');
	
	$pop->{dbID}    = $dbID;
	$pop->{adaptor} = $self;
	
	# add entry to population table also
	$sth = $dbh->prepare(q{
		INSERT INTO population (sample_id) VALUES (?)
	});
	$sth->execute($dbID);
	$sth->finish;
}

=head2 fetch_population_by_synonym

    Arg [1]              : $population_synonym
    Example              : my $pop = $pop_adaptor->fetch_population_by_synonym($population_synonym,$source);
    Description          : Retrieves populations for the synonym given in the source. If no source is provided, retrieves all the synonyms
    Returntype           : list of Bio::EnsEMBL::Variation::Population
    Exceptions           : none
    Caller               : general
    Status               : Stable

=cut

sub fetch_population_by_synonym{
    my $self = shift;
    my $synonym_name = shift;
    my $source = shift;
    my $pops;
    my $pop;
    #return all sample_id from the database
    my $samples = $self->SUPER::fetch_sample_by_synonym($synonym_name, $source);
    foreach my $sample_id (@{$samples}){
	#get the ones that are individuals
	$pop = $self->fetch_by_dbID($sample_id);
	push @{$pops}, $pop if (defined $pop);
    }
    return $pops;
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

  throw('name argument expected') if(!defined($name));

  my $sth = $self->prepare(q{SELECT p.sample_id, s.name, s.size, s.description, s.freqs_from_gts
                             FROM   population p, sample s
                             WHERE  s.name = ?
			     AND    s.sample_id = p.sample_id});

  $sth->bind_param(1,$name,SQL_VARCHAR);
  $sth->execute();

  my $result = $self->_objs_from_sth($sth);

  $sth->finish();

  return undef if(!@$result);

  return $result->[0];
}

=head2 fetch_all_by_dbID_list

  Arg [1]    : listref $list
  Example    : $pops = $pop_adaptor->fetch_all_by_dbID_list([907,1132]);
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

  if(!defined($list) || ref($list) ne 'ARRAY') {
    throw("list reference argument is required");
  }
  
  return [] unless scalar @$list >= 1;
  
  my $id_str = (@$list > 1)  ? " IN (".join(',',@$list).")"   :   ' = \''.$list->[0].'\'';

  my $sth = $self->prepare(qq{SELECT p.sample_id, s.name, s.size, s.description, s.freqs_from_gts
                              FROM   population p, sample s
                              WHERE  s.sample_id $id_str AND s.sample_id = p.sample_id});
  $sth->execute();

  my $result = $self->_objs_from_sth($sth);

  $sth->finish();

  return undef if(!@$result);

  return $result;
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

  throw('name argument expected') if(!defined($name));

  my $sth = $self->prepare(q{SELECT p.sample_id, s.name, s.size, s.description, s.freqs_from_gts
                             FROM   population p, sample s
                             WHERE  s.name like concat('%', ?, '%')
			     AND    s.sample_id = p.sample_id});

  $sth->bind_param(1,$name,SQL_VARCHAR);
  $sth->execute();

  my $result = $self->_objs_from_sth($sth);

  $sth->finish();

  return $result;
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
  Status     : At Risk

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

  my $sth = $self->prepare(q{SELECT p.sample_id, s.name, s.size,
                                    s.description, s.freqs_from_gts
                             FROM   population p, population_structure ps, sample s
                             WHERE  p.sample_id = ps.sub_population_sample_id
                             AND    ps.super_population_sample_id = ?
			     AND    p.sample_id = s.sample_id});

  $sth->bind_param(1,$pop->dbID,SQL_INTEGER);
  $sth->execute();

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
  Status     : At Risk

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

  my $sth = $self->prepare(q{SELECT p.sample_id, s.name, s.size,
                                    s.description, s.freqs_from_gts
                             FROM   population p, population_structure ps, sample s
                             WHERE  p.sample_id = ps.super_population_sample_id
                             AND    ps.sub_population_sample_id = ?
			     AND    p.sample_id = s.sample_id});

  $sth->bind_param(1,$pop->dbID,SQL_INTEGER);
  $sth->execute();

  my $result = $self->_objs_from_sth($sth);

  $sth->finish();

  return $result;
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

sub fetch_default_LDPopulation{
    my $self = shift;
    my $population_id;
    
    my $sth = $self->prepare(qq{SELECT meta_value from meta where meta_key = ?});

    $sth->bind_param(1,'pairwise_ld.default_population',SQL_VARCHAR);
    $sth->execute();
    $sth->bind_columns(\$population_id);
    $sth->fetch();
    $sth->finish;

    if (defined $population_id){
	return $self->fetch_by_dbID($population_id);
    }
    else{
	return undef;
    }
}


=head2 fetch_all_LD_Populations

    Example     : @populations = @{$pop_adaptor->fetch_all_LD Populations();
    Description : Gets all populations that can be used in the LD display
    ReturnType  : listref of Bio::EnsEMBL::Variation::Population objects
    Exceptions  : none
    Caller      : general
    Status      : At Risk

=cut

sub fetch_all_LD_Populations{
    my $self = shift;
	
	return [grep {$_->name !~ /ALL|AFR|AMR|ASN|EUR/} @{$self->generic_fetch(qq{ s.display = 'LD' })}];
}


=head2 fetch_all_HapMap_Populations

    Example     : @populations = @{$pop_adaptor->fetch_all_HapMap_populations();
    Description : Gets all populations from the HapMap project (human only!)
    ReturnType  : listref of Bio::EnsEMBL::Variation::Population objects
    Exceptions  : none
    Caller      : general
    Status      : At Risk

=cut

sub fetch_all_HapMap_Populations {
    my $self = shift;
	
	return $self->generic_fetch(qq{ s.name like 'cshl-hapmap%' });
}


=head2 fetch_all_1KG_Populations

    Example     : @populations = @{$pop_adaptor->fetch_all_1KG_populations();
    Description : Gets all populations from the 1000 genomes project (human only!)
    ReturnType  : listref of Bio::EnsEMBL::Variation::Population objects
    Exceptions  : none
    Caller      : general
    Status      : At Risk

=cut

sub fetch_all_1KG_Populations{
    my $self = shift;
	
	return $self->generic_fetch(qq{ s.name like '1000GENOMES%' });
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
  Status      : Stable

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

    my $sth = $self->prepare(qq{SELECT p.sample_id, s.name, s.size, s.description, s.freqs_from_gts
				FROM population p, individual_population ip, sample s
				WHERE s.sample_id = ip.population_sample_id
				AND s.sample_id = p.sample_id
                                AND ip.individual_sample_id = ?
			    });
    $sth->bind_param(1,$ind->dbID,SQL_INTEGER);
    $sth->execute();

    my $results = $self->_objs_from_sth($sth);

    $sth->finish();

    return $results;
}


=head2 fetch_all_by_Individual_list

  Arg [1]     : reference to list of of Bio::EnsEMBL::Variation::Individual objects
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

sub fetch_all_by_Individual_list{
	my $self = shift;
	my $list = shift;
	
	if(!ref($list) || !$list->[0]->isa('Bio::EnsEMBL::Variation::Individual')) {
		throw("Listref of Bio::EnsEMBL::Variation::Individual arg expected");
	}
	
	if(!$list->[0]->dbID()) {
		warning("First Individual does not have dbID, cannot retrieve Populations");
		return [];
	}
	
	my $id_str = " IN (" . join(',', map {$_->dbID} @$list). ")";	
	
	my $sth = $self->prepare(qq{
		SELECT p.sample_id, s.name, s.size, s.description, s.freqs_from_gts
		FROM population p, individual_population ip, sample s
		WHERE s.sample_id = ip.population_sample_id
		AND s.sample_id = p.sample_id
		AND ip.individual_sample_id $id_str
	});
	$sth->execute();
	
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
  Status      : At Risk

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

    my $sth = $self->prepare(qq{
		SELECT p.sample_id, s.name, s.size, s.description, s.freqs_from_gts
		FROM population p, tagged_variation_feature tvf, sample s
		WHERE p.sample_id = tvf.sample_id
		AND   s.sample_id = p.sample_id
		AND   tvf.tagged_variation_feature_id = ?
	});
    $sth->bind_param(1,$variation_feature->dbID,SQL_INTEGER);
    $sth->execute();
    my $results = $self->_objs_from_sth($sth);

    $sth->finish();

    return $results;
}

=head2 fetch_tag_Population

  Arg [1]     : Bio::EnsEMBL::Variation::VariationFeature $vf
  Example     : my $vf = $vf_adaptor->fetch_by_name('rs205621');
                my $populations_is_tag = $vf->is_tag();
                foreach my $pop (@{$vf_adaptor->is_tag}){
		    print $pop->name," has been tagged using a 0.99 r2 criteria\n";
                }
  Description : Retrieves all populations in which the specified variation feature is a tag
  ReturnType  : reference to list of Bio::EnsEMBL::Variation::Population objects
  Exceptions  : throw if incorrect argument is passed
                warning if provided variation feature does not have a dbID
  Caller      : general
  Status      : At Risk

=cut

sub fetch_tag_Population{
    my $self = shift;
    my $variation_feature = shift;

    if(!ref($variation_feature) || !$variation_feature->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
	throw("Bio::EnsEMBL::Variation::VariationFeature arg expected");
    }
    
    if(!$variation_feature->dbID()) {
	warning("Variation feature does not have dbID, cannot retrieve tag populations");
	return [];
  } 

    my $sth = $self->prepare(qq{
		SELECT p.sample_id, s.name, s.size, s.description, s.freqs_from_gts
		FROM population p, tagged_variation_feature tvf, sample s
		WHERE p.sample_id = tvf.sample_id
		AND   s.sample_id = p.sample_id
		AND   tvf.variation_feature_id = ?
	});
    $sth->bind_param(1,$variation_feature->dbID,SQL_INTEGER);
    $sth->execute();
    my $results = $self->_objs_from_sth($sth);

    $sth->finish();

    return $results;
}


=head2 get_sample_id_for_population_names

  Arg [1]     : $population_names reference to list of population names
  Example     : my $ids = $pop_adaptor->get_sample_id_for_population_names(['CSHL-HAPMAP:HAPMAP-MEX','1000GENOMES:pilot_1_CHB+JPT_low_coverage_panel']);
                map {printf("Population: \%s has sample_id \%d\n",$ids->{$_},$_)} keys(%{$ids});
  Description : Retrieve the sample_ids for a list of population names
  ReturnType  : reference to hash with sample_ids as keys and population names as values
  Caller      : web
  Status      : At Risk

=cut

sub get_sample_id_for_population_names {
    my $self = shift;
    my $population_names = shift;

    # Wrap the argument into an arrayref
    $population_names = wrap_array($population_names);
    
    # Define a statement handle for the lookup query
    my $stmt = qq{
        SELECT
            s.sample_id
            s.name
        FROM
            sample s
        WHERE
            s.name = ?
        LIMIT 1
    };
    my $sth = $self->prepare($stmt);
    
    # Loop over the population names and query the db
    my %sample_ids;
    foreach my $name (@{$population_names}) {
        $sth->execute($name);
        
        my ($sid,$sname);
        $sth->bind_columns(\$sid,\$sname);
        $sth->execute();
        
        $sample_ids{$sid} = $sname if (defined($sid));
    }
    
    return \%sample_ids;
}

#
# private method, creates population objects from an executed statement handle
# ordering of columns must be consistant
#
sub _objs_from_sth {
  my $self = shift;
  my $sth  = shift;

  my @pops;

  my ($pop_id, $name, $size, $desc, $freqs);
  
  $DB::single = 1;

  $sth->bind_columns(\$pop_id, \$name, \$size, \$desc, \$freqs);

  while($sth->fetch()) {
	
    push @pops, Bio::EnsEMBL::Variation::Population->new
      (-dbID => $pop_id,
       -ADAPTOR => $self,
       -NAME => $name,
       -DESCRIPTION => $desc,
       -SIZE => $size,
       -FREQS => $freqs);
  }

  return \@pops;
}

sub _tables{return (['population','p'],
		    ['sample','s']);}

sub _columns{
    return qw(s.sample_id s.name s.size s.description s.freqs_from_gts);
}

sub _default_where_clause{
    return 's.sample_id = p.sample_id';
}

1;
