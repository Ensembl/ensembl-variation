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
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeAdaptor
#
# Copyright (c) 2005 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeAdaptor

=head1 SYNOPSIS

Adaptor for IndividualGenotype objects.

=head1 DESCRIPTION

This adaptor provides database connectivity for IndividualGenotype objects.
IndividualGenotypes may be retrieved from the Ensembl variation database by
several means using this module.

=head1 METHODS

=cut
package Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeAdaptor;

use strict;
use warnings;

use vars qw(@ISA);

use Bio::EnsEMBL::Variation::DBSQL::BaseGenotypeAdaptor;
use Bio::EnsEMBL::Variation::IndividualGenotype;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

@ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseGenotypeAdaptor');

# stores a listref of individual genotype objects
sub store {
	my ($self, $gts, $merge) = @_;
	
	throw("First argument to store is not a listref") unless ref($gts) eq 'ARRAY';
	
	# sort genotypes into rows by variation
	my %by_var;
	push @{$by_var{$_->variation->dbID.'_'.($_->{subsnp} ? $_->{subsnp} : '')}}, $_ for @$gts;
	
	# get unique genotypes and codes
	my %unique_gts = map {$_->genotype_string() => 1} @$gts;
	$unique_gts{$_} = $self->_genotype_code([split /\|/, $_]) for keys %unique_gts;
	
	# get variation objects
	my %var_objs = map {$_->variation->dbID => $_->variation} @$gts;
	
	my $dbh = $self->dbc->db_handle;
	
	my $sth = $dbh->prepare_cached(qq{
		INSERT INTO compressed_genotype_var(
			variation_id,
			subsnp_id,
			genotypes
		) VALUES (?,?,?)
	});
	
	my $update_sth = $dbh->prepare(qq{
		UPDATE compressed_genotype_var
		SET genotypes = ?
		WHERE variation_id = ?
		AND subsnp_id = ?
	});
	
	my $rows_added = 0;
	
	foreach my $combo_id(keys %by_var) {
		my $genotype_string = '';
		
		my ($var_id, $subsnp_id) = split /\_/, $combo_id;
		$subsnp_id = undef if $subsnp_id eq '';
		
		# check for existing genotypes
		my @existing_gts = ($merge ?
			grep {
				(defined($_->{subsnp}) && defined($subsnp_id) && $_->{subsnp} eq $subsnp_id) ||
				(!defined($_->{subsnp}) && !defined($subsnp_id))
			} @{$self->fetch_all_by_Variation($var_objs{$var_id})} : ());
		
		# update if existing
		if(@existing_gts) {
			
			# refresh unique_gts
			%unique_gts = map {$_->genotype_string() => 1} (@existing_gts, @$gts);
			$unique_gts{$_} = $self->_genotype_code([split /\|/, $_]) for keys %unique_gts;
			
			# make sure we don't put in duplicates
			my %by_ind = map {$_->individual->dbID => $_} (@existing_gts, @$gts);
			
			$genotype_string .= pack("ww", $_->individual->dbID, $unique_gts{$_->genotype_string}) for values %by_ind;
			
			$update_sth->execute(
				$genotype_string,
				$var_id,
				$subsnp_id
			);
		}
		
		else {
			$genotype_string .= pack("ww", $_->individual->dbID, $unique_gts{$_->genotype_string}) for @$gts;
			
			$sth->execute(
				$var_id,
				$subsnp_id,
				$genotype_string
			);
			
			$rows_added++;
		}
		
	}
	
	$sth->finish;
	$update_sth->finish;
	
	return $rows_added;
}


=head2 fetch_all_by_Variation

  Arg [1]    : Bio::EnsEMBL::Variation $variation
  Example    : my $var = $variation_adaptor->fetch_by_name( "rs1121" )
               $igtypes = $igtype_adaptor->fetch_all_by_Variation( $var )
  Description: Retrieves a list of individual genotypes for the given Variation.
               If none are available an empty listref is returned.
  Returntype : listref Bio::EnsEMBL::Variation::IndividualGenotype 
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut


sub fetch_all_by_Variation {
    my $self = shift;
    my $variation = shift;
	my $individual = shift;

    if(!ref($variation) || !$variation->isa('Bio::EnsEMBL::Variation::Variation')) {
		throw('Bio::EnsEMBL::Variation::Variation argument expected');
    }

    if(!defined($variation->dbID())) {
		warning("Cannot retrieve genotypes for variation without dbID");
		return [];
    }
	
	my $results;
	
	# check cache
	if(!defined($self->{_cache}->{$variation->dbID})) {
		$self->{_cache}->{$variation->dbID} = $self->generic_fetch("g.variation_id = " . $variation->dbID());
		$_->variation($variation) for @{$self->{_cache}->{$variation->dbID}};
	}
	
	if (defined $individual && defined $individual->dbID){
		if($individual->isa('Bio::EnsEMBL::Variation::Individual')) {
			@$results = grep {$_->individual->dbID == $individual->dbID} @{$self->{_cache}->{$variation->dbID}};
		}
		elsif($individual->isa('Bio::EnsEMBL::Variation::Population')) {
			my %include = map {$_->dbID => 1} @{$individual->get_all_Individuals};			
			@$results = grep {$include{$_->individual->dbID}} @{$self->{_cache}->{$variation->dbID}};
		}
		else {
			throw("Argument supplied is not of type Bio::EnsEMBL::Variation::Sample");
		}
	}
	else {
		$results = $self->{_cache}->{$variation->dbID};
	}
	
	# flip genotypes for flipped variations
	#if(defined $variation->flipped && $variation->flipped == 1) {
	#	foreach my $gt(@$results) {
	#		my @new_gt;
	#		
	#		foreach my $allele(@{$gt->{genotype}}) {
	#			reverse_comp(\$allele);
	#			push @new_gt, $allele;
	#		}
	#		$gt->{genotype} = \@new_gt;
	#	}
	#}
	
	return $results;
}

sub fetch_all_by_Slice {
	my $self = shift;
	
	my $cga = $self->db->get_IndividualGenotypeFeatureAdaptor();
	
	return $cga->fetch_all_by_Slice(@_);
}

sub _tables{
    my $self = shift;

	return (['compressed_genotype_var','g'],['failed_variation','fv']);
}

#ÊAdd a left join to the failed_variation table
sub _left_join { return ([ 'failed_variation', 'fv.variation_id = g.variation_id']); }

sub _columns{
    return qw(g.variation_id g.subsnp_id g.genotypes);
}

sub _objs_from_sth{
	my $self = shift;
	my $sth = shift;
	
	my ($variation_id, $subsnp_id, $genotypes);
	
	$sth->bind_columns(\$variation_id, \$subsnp_id, \$genotypes);
	
	my (%individual_hash, %gt_code_hash, @results);
	
	while($sth->fetch) {
		my @genotypes = unpack("(ww)*", $genotypes);
		
		while(@genotypes) {
			my $sample_id = shift @genotypes;
			my $gt_code = shift @genotypes;
			
			my $igtype  = Bio::EnsEMBL::Variation::IndividualGenotype->new_fast({
				_variation_id => $variation_id,
				subsnp        => $subsnp_id,
				adaptor       => $self,
			});
			
			$individual_hash{$sample_id} ||= [];
			push @{$individual_hash{$sample_id}}, $igtype;
			
			$gt_code_hash{$gt_code} ||= [];
			push @{$gt_code_hash{$gt_code}}, $igtype;
			
			push @results, $igtype;
		}
	}
	
	# fetch individuals
	my $ia = $self->db()->get_IndividualAdaptor();
	my $inds = $ia->fetch_all_by_dbID_list([keys %individual_hash]);
	
	foreach my $i (@$inds) {
		foreach my $igty (@{$individual_hash{$i->dbID()}}) {
			$igty->{individual} = $i;
		}
	}
	
	# get all genotypes from codes
	my $gtca = $self->db->get_GenotypeCodeAdaptor();
	my $gtcs = $gtca->fetch_all_by_dbID_list([keys %gt_code_hash]);
	
	foreach my $gtc(@$gtcs) {
		foreach my $igty(@{$gt_code_hash{$gtc->dbID}}) {
			$igty->{genotype} = $gtc->genotype;
		}
	}
	
	# check for defined genotype and individual
	@results = grep {defined $_->{genotype} && defined $_->{individual}} @results;
	
	return \@results;
}

1;
