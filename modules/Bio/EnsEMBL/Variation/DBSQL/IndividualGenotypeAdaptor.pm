#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeAdaptor

=head1 SYNOPSIS

  $db = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(...);

  $iga = $db->get_IndividualGenotypeAdaptor();
  $ia  = $db->get_IndividualAdaptor();

  # Get an IndividualGenotype by its internal identifier
  $igtype = $ia->fetch_by_dbID(145);

  # Get all individual genotypes for an individual
  $ind = $ia->fetch_all_by_Individual(1219);

  foreach $igtype (@{$iga->fetch_all_by_Individual($ind)}) {
    print $igtype->variation()->name(),  ' ',
          $igtype->allele1(), '/', $igtype->allele2(), "\n";
  }



=head1 DESCRIPTION

This adaptor provides database connectivity for IndividualGenotype objects.
IndividualGenotypes may be retrieved from the Ensembl variation database by
several means using this module.

=head1 AUTHOR - Graham McVicker

=head1 CONTACT

Post questions to the Ensembl development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Variation::IndividualGenotype;

use Data::Dumper;
our @ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');


=head2 fetch_all_by_Individual

  Arg [1]    : Bio::EnsEMBL::Variation::Individual
  Example    : $ind = $ind_adaptor->fetch_by_dbID(1345);
               @gtys = $igty_adaptor->fetch_all_by_Individual($ind);
  Description: Retrieves all genotypes which are stored for a specified
               individual.
  Returntype : Bio::EnsEMBL::Variation::IndividualGenotype
  Exceptions : throw on incorrect argument
  Caller     : general

=cut

sub fetch_all_by_Individual {
  my $self = shift;
  my $ind = shift;

  if(!ref($ind) || !$ind->isa('Bio::EnsEMBL::Variation::Individual')) {
    throw('Bio::EnsEMBL::Variation::Individual argument expected');
  }

  if(!defined($ind->dbID())) {
    warning("Cannot retrieve genotypes for individual without set dbID");
    return [];
}

  $self->_tables(['individual_genotype_single_bp','igs']);
  my $res = $self->generic_fetch("individual_id = " . $ind->dbID()); #to select data from individual_genotype_single_bp
  $self->_tables(['individual_genotype_multiple_bp','igm']);
  push @{$res},@{$self->generic_fetch("individual_id = " . $ind->dbID())}; #to select data from individual_genotype_multiple_bp  
  return $res
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

=cut


sub fetch_all_by_Variation {
    my $self = shift;
    my $variation = shift;

    if(!ref($variation) || !$variation->isa('Bio::EnsEMBL::Variation::Variation')) {
	throw('Bio::EnsEMBL::Variation::Variation argument expected');
    }

    if(!defined($variation->dbID())) {
	warning("Cannot retrieve genotypes for variation without set dbID");
	return [];
    }	
    $self->_tables(['individual_genotype_single_bp','igs']);
    my $res = $self->generic_fetch("variation_id = " . $variation->dbID()); #to select data from individual_genotype_single_bp
    $self->_tables(['individual_genotype_multiple_bp','igm']);
    push @{$res},@{$self->generic_fetch("variation_id = " . $variation->dbID())}; #to select data from individual_genotype_multiple_bp
    return $res;
}

#Getter/Setter of the table where to select the data from

sub _tables{
    my $self = shift;
    return $self->{'_table'} = shift if (@_);
    return $self->{'_table'};
}

sub _columns{
    return qw(individual_id variation_id allele_1 allele_2);
}

sub _objs_from_sth{
    my $self = shift;
    my $sth = shift;
    
    my @results;
    my ($individual_id, $variation_id, $allele_1, $allele_2);
    $sth->bind_columns(\$individual_id, \$variation_id, \$allele_1, \$allele_2);
    
    my %individual_hash;
    my %variation_hash;
    while($sth->fetch()){
	my $igtype = Bio::EnsEMBL::Variation::IndividualGenotype->new
	    (-adaptor => $self,
	     -allele1 => $allele_1,
	     -allele2 => $allele_2);
	$individual_hash{$individual_id} ||= [];
	$variation_hash{$variation_id} ||= [];
	push @{$variation_hash{$variation_id}}, $igtype; #store the variations to get the objects once
	push @{$individual_hash{$individual_id}}, $igtype; #store the individuals to get the objects once
	push @results, $igtype;
    }

    # get all variations in one query (faster)
    # and add to already created genotypes
    my @var_ids = keys %variation_hash;
    my $va = $self->db()->get_VariationAdaptor();
    my $vars = $va->fetch_all_by_dbID_list(\@var_ids);
    
    foreach my $v (@$vars) {
	foreach my $igty (@{$variation_hash{$v->dbID()}}) {
	    $igty->variation($v);
	}
    }

    # get all individual in one query (faster)
    # and add to already created genotypes
    my @ind_ids = keys %individual_hash;
    my $ia = $self->db()->get_IndividualAdaptor();
    my $inds = $ia->fetch_all_by_dbID_list(\@ind_ids);
    
    foreach my $i (@$inds) {
	foreach my $igty (@{$individual_hash{$i->dbID()}}) {
	    $igty->individual($i);
	}
    }
    return \@results;   

}

1;
