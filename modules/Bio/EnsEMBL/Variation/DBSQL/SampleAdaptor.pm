=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/legal/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::SampleAdaptor
#
# Copyright (c) 2005 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::SampleAdaptor

=head1 SYNOPSIS

Abstract class - should not be instantiated.  Implementation of
abstract methods must be performed by subclasses.

Base adaptors provides:

#using the adaptor of the subclass, given the id of the sample, returns all synonyms in database
$synonyms = $sample_adaptor->fetch_synonyms($sample_id);

#using the adaptor of the subclass and given the name of the synonym, returns the sample
$sample = $sample_adaptor->fetch_sample_by_synonym($sample_synonym_id);

=head1 DESCRIPTION

This is a base class implementing common methods in population, individual and strain. This base
class is simply a way of merging similar concepts that should have the same ID

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::SampleAdaptor;
use vars qw(@ISA @EXPORT);
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);


our @ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');
@EXPORT = (@{$DBI::EXPORT_TAGS{'sql_types'}});

=head2 fetch_synonyms

    Arg [1]              : $sample_id
    Arg [2] (optional)   : $source
    Example              : my $dbSNP_synonyms = $pop_adaptor->fetch_synonyms($sample_id,$dbSNP);
                           my $all_synonyms = $pop_adaptor->fetch_synonyms($sample_id);
    Description: Retrieves synonyms for the source provided. Otherwise, return all the synonyms for the sample
    Returntype : list of strings
    Exceptions : none
    Caller     : Bio:EnsEMBL:Variation::Sample
    Status     : Stable

=cut

sub fetch_synonyms{
    my $self = shift;
    my $dbID = shift;
    my $source = shift;
    my $sample_synonym;
    my $synonyms;

    my $sql;
    if (defined $source){
	$sql = qq{SELECT ss.name FROM sample_synonym ss, source s WHERE ss.sample_id = ? AND ss.source_id = s.source_id AND s.name = "$source"}
    }
    else{
	$sql = qq{SELECT name FROM sample_synonym WHERE sample_id = ?};
    }
    my $sth = $self->prepare($sql);
    $sth->bind_param(1,$dbID,SQL_INTEGER);
    $sth->execute();
    $sth->bind_columns(\$sample_synonym);
    while ($sth->fetch){
	push @{$synonyms},$sample_synonym;
    }
    return $synonyms;
}

=head2 fetch_sample_by_synonym

    Arg [1]              : $sample_synonym
    Example              : my $pop = $pop_adaptor->fetch_sample_by_synonym($sample_synonym,$source);
    Description          : Retrieves sample for the synonym given in the source. If no source is provided, retrieves all the synonyms
    Returntype           : list of Bio::EnsEMBL::Variation::Sample
    Exceptions           : none
    Caller               : general
    Status               : At Risk

=cut

sub fetch_sample_by_synonym{
    my $self = shift;
    my $synonym_name = shift;
    my $source = shift;
    my $sql;
    my $sample;
    my $sample_array;
    if (defined $source){
	$sql = qq{SELECT sample_id FROM sample_synonym ss, source s WHERE ss.name = ? and ss.source_id = s.source_id AND s.name = "$source"};
    }
    else{
	$sql = qq{SELECT sample_id FROM sample_synonym WHERE name = ?};
    }
    my $sample_id;
    my $sth = $self->prepare($sql);
    $sth->bind_param(1,$synonym_name,SQL_VARCHAR);
    $sth->execute();    
    $sth->bind_columns(\$sample_id);
    while ($sth->fetch()){
	push @{$sample_array}, $sample_id;
    }
    return $sample_array;
    
}

=head2 fetch_by_dbID

  Arg [1]    : int $id
               The unique sample identifier for the sample to be obtained
  Example    : $population = $population_adaptor->fetch_by_dbID(1234);
  Description: Returns the feature sample from the database defined by the
               the id $id.  
  Returntype : Bio::EnsEMBL::Variation::Sample
  Exceptions : thrown if $id arg is not provided
               does not exist
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_dbID{
  my ($self,$id) = @_;

  throw("id argument is required") if(!defined $id);

  #construct a constraint like 't1.table1_id = 123'
  my @tabs = $self->_tables;
  my ($name, $syn) = @{$tabs[0]};
  #the constraint must contain the sample_id, that it is used either for individuals or populations
  my $constraint = "${syn}.sample_id = ?";
  $self->bind_param_generic_fetch($id,SQL_INTEGER);

  #Should only be one
  my ($feat) = @{$self->generic_fetch($constraint)};

  return undef if(!$feat);

  return $feat;
}


=head2 fetch_all_by_dbID_list

  Arg [1]    : listref of ints $id_list
               The unique database identifiers for the samples to be obtained
  Example    : @individuals = @{$individual_adaptor->fetch_by_dbID_list([1234, 2131, 982]))};
  Description: Returns the samples created from the database defined by the
               the ids in contained in the id list $id_list.  If none of the
               samples are found in the database a reference to an empty 
               list is returned.
  Returntype : listref of Bio::EnsEMBL::Variation::Sample
  Exceptions : thrown if $id arg is not provided
               does not exist
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_dbID_list {
  my ($self,$id_list_ref) = @_;

  if(!defined($id_list_ref) || ref($id_list_ref) ne 'ARRAY') {
    throw("id_list list reference argument is required");
  }

  return [] if(!@$id_list_ref);

  my @out;
  #construct a constraint like 't1.table1_id = 123'
  my @tabs = $self->_tables;
  my ($name, $syn) = @{$tabs[0]};

  # mysql is faster and we ensure that we do not exceed the max query size by
  # splitting large queries into smaller queries of 200 ids
  my $max_size = 200;
  my @id_list = @$id_list_ref;

  while(@id_list) {
    my @ids;
    if(@id_list > $max_size) {
      @ids = splice(@id_list, 0, $max_size);
    } else {
      @ids = splice(@id_list, 0);
    }

    my $id_str;
    if(@ids > 1)  {
      $id_str = " IN (" . join(',', @ids). ")";
    } else {
      $id_str = " = ?";
      $self->bind_param_generic_fetch($ids[0],SQL_INTEGER);
    }

    my $constraint = "${syn}.sample_id $id_str";

    push @out, @{$self->generic_fetch($constraint)};
  }

  return \@out;
}


sub _get_individual_population_hash {
	my $self = shift;
	my $id_list_ref = shift;

	if(!defined($id_list_ref) || ref($id_list_ref) ne 'ARRAY') {
		throw("id_list list reference argument is required");
	}
	
	return [] if (!@$id_list_ref);
	
	my %ip_hash;
	my $max_size = 200;
	my @id_list = @$id_list_ref;
	
	while(@id_list) {
		my @ids;
		if(@id_list > $max_size) {
			@ids = splice(@id_list, 0, $max_size);
		} else {
			@ids = splice(@id_list, 0);
		}
		
		my $id_str;
		if(@ids > 1)  {
			$id_str = " IN (" . join(',', @ids). ")";
		} else {
			$id_str = " = ".$ids[0];
		}
		
		my $sth = $self->prepare(qq/
			SELECT individual_sample_id, population_sample_id
			FROM individual_population
			WHERE individual_sample_id $id_str
		/);
		
		$sth->execute();
		
		my ($ind, $pop);
		$sth->bind_columns(\$ind, \$pop);
		$ip_hash{$pop}{$ind} = 1 while $sth->fetch;
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

sub _get_sample_name_by_dbID {
  my $self = shift;
  my $dbID = shift;
  
  my $sth = $self->dbc->prepare(qq{
    SELECT name
    FROM sample
    WHERE sample_id = ?
  });
  $sth->bind_param(1,$dbID,SQL_INTEGER);
  $sth->execute();
  
  my $sample_name;
  $sth->bind_columns(\$sample_name);
  $sth->fetch;
  $sth->finish;
  
  return $sample_name;
}

1;
