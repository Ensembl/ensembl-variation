=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
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

use Scalar::Util qw(weaken);

our $CACHE_SIZE = 5;

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

# similar to store above but writes to old-style non-compressed table
# defaults to individual_genotype_multiple_bp since
# tmp_individual_genotype_single_bp will only accept single nucleotide genotypes
sub store_uncompressed {
	my ($self, $gts, $table) = @_;
	
	$table ||= 'individual_genotype_multiple_bp';
	
	throw("First argument to store is not a listref") unless ref($gts) eq 'ARRAY';
	
	my $dbh = $self->dbc->db_handle;
	
	my $q_string = join ",", map {'(?,?,?,?,?)'} @$gts;
	
	my @args = map {
		$_->{_variation_id} || $_->variation->dbID,
		$_->{subsnp},
		$_->allele1,
		$_->allele2,
		$_->individual ? $_->individual->dbID : undef,
	} @$gts;
	
	my $sth = $dbh->prepare_cached(qq{
		INSERT INTO $table(
			variation_id,
			subsnp_id,
			allele_1,
			allele_2,
			sample_id
		) VALUES $q_string
	});
	
	$sth->execute(@args);
	
	$sth->finish;
	
	return scalar @$gts;
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
  Status     : Stable

=cut


sub fetch_all_by_Variation {
    my $self = shift;
    my $variation = shift;
	my $sample = shift;

    if(!ref($variation) || !$variation->isa('Bio::EnsEMBL::Variation::Variation')) {
		throw('Bio::EnsEMBL::Variation::Variation argument expected');
    }

    if(!defined($variation->dbID())) {
		warning("Cannot retrieve genotypes for variation without dbID");
		return [];
    }
	
	my $results;
	
	# check cache
	# All ind gt objects for a variation are stored on the adaptor.
	# They are stored in an array of limited size so that the memory usage
	# doesn't go out of control for when the API requests genotypes for many
	# variations sequentially.
	my $cached;
	
	if(defined($self->{_cache})) {
		foreach my $stored(@{$self->{_cache}}) {
			my @keys = keys %{$stored};
			$cached = $stored->{$keys[0]} if $keys[0] eq $variation->dbID;
			last if defined($cached);
		}
	}

	if(!defined($cached)) {
		$cached = $self->generic_fetch("g.variation_id = " . $variation->dbID());
		for(@$cached) {
			$_->variation($variation);
			weaken($_->{'variation'});
		}
		
		# add genotypes for this variant to the cache
		push @{$self->{_cache}}, {$variation->dbID => $cached};
		
		# shift off first element to keep cache within size limit
		shift @{$self->{_cache}} if scalar @{$self->{_cache}} > $CACHE_SIZE;
	}
	
	if (defined $sample && defined $sample->dbID){
		if($sample->isa('Bio::EnsEMBL::Variation::Individual')) {
			@$results = grep {$_->individual->dbID == $sample->dbID} @{$cached};
		}
		elsif($sample->isa('Bio::EnsEMBL::Variation::Population')) {
			my %include = map {$_->dbID => 1} @{$sample->get_all_Individuals};			
			@$results = grep {$include{$_->individual->dbID}} @{$cached};
		}
		else {
			throw("Argument supplied is not of type Bio::EnsEMBL::Variation::Sample");
		}
	}
	
	return $cached;
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
