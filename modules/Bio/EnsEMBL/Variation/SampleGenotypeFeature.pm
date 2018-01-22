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

# Ensembl module for Bio::EnsEMBL::Variation::SampleGenotypeFeature
#
#


=head1 NAME

Bio::EnsEMBL::Variation::SampleGenotypeFeature-Module representing the genotype
of a single sample at a single position

=head1 SYNOPSIS

    print $genotype->variation()->name(), "\n";
    print $genotype->allele1(), '/', $genotype->allele2(), "\n";
    print $genotype->sample()->name(), "\n";

=head1 DESCRIPTION

This is a class representing the genotype of a single diploid sample at
a specific position

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::SampleGenotypeFeature;

use Bio::EnsEMBL::Variation::SampleGenotype;
use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Variation::SampleGenotype Bio::EnsEMBL::Feature);



=head2 new

  Arg [-adaptor] :
    Bio::EnsEMBL::Variation::DBSQL::IndividualAdaptor
  Arg [-START] :
    see superclass constructor
  Arg [-END] :
    see superclass constructor
  Arg [-STRAND] :
    see superclass constructor
  Arg [-SLICE] :
    see superclass constructor
  Arg [-genotype] :
    arrayref - arrayref of alleles making up this genotype (in haplotype order)
  Arg [-variation] :
    Bio::EnsEMBL::Variation::Variation - The variation associated with this
    genotype
  Arg [-sample] :
    Bio::EnsEMBL::Sample - The sample this genotype is for.
  Example    : $ind_genotype = Bio:EnsEMBL::Variation::SampleGenotype->new(
					-start      => 100,
					-end        => 100,
					-strand     => 1,
					-slice      => $slice,
					-genotype   => ['A','T'],
					-variation  => $variation,
					-sample => $sample
				);
  Description: Constructor.  Instantiates SampleGenotype object.
  Returntype : Bio::EnsEMBL::Variation::SampleGenotype
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub new {
	my $caller = shift;
	my $class = ref($caller) || $caller;
	
	my $self = $class->SUPER::new(@_);
	
	my ($adaptor, $genotype, $var, $var_id, $sample) =
	rearrange([qw(adaptor genotype variation _variation_id sample)],@_);
	
	if(defined($var) && (!ref($var) || !$var->isa('Bio::EnsEMBL::Variation::Variation'))) {
		throw("Bio::EnsEMBL::Variation::Variation argument expected");
	}
	
	if(defined($sample) && (!ref($sample) || !$sample->isa('Bio::EnsEMBL::Variation::Sample'))) {
		throw("Bio::EnsEMBL::Variation::Sample argument expected");
	}
	
	$self->{'adaptor'} = $adaptor;
	$self->{'genotype'} = $genotype;
	$self->{'sample'} = $sample;
	$self->{'variation'} = $var;
	$self->{'_variation_id'} = $var_id;
	
	return $self;
}

sub variation_feature {
  my $self = shift;

  $self->{variation_feature} = shift if @_;

  if(!defined($self->{variation_feature})) {
    $self->{variation_feature} = (grep {$_->{start} == $self->{start} && $_->seq_region_name cmp $self->{slice}->seq_region_name} @{$self->variation->get_all_VariationFeatures})[0];
  }
 
  return $self->{variation_feature};
}

=head2 differences
  Arg [1]    : Hashref $differences (optional)
  Example    : $differences = $sample_genotype_feature->differences();
               while (my ($sample_name, $genotype_string) = each %$differences) {
                 print "$sample_name $genotype_string\n";
               }  
  Description: Getter/Setter for Hashref of differences as computed in
               Bio::EnsEMBL::Variation::DBSQL::SampleGenotypeFeatureAdaptor::fetch_all_differences_by_Slice.
               Keys are the sample names and corresponding values are the genotype strings.
  Returntype : Hashref
  Exceptions : None
  Caller     : General
  Status     : Stable
=cut
sub differences {
  my $self = shift;
  $self->{differences} = shift if @_;
  return $self->{differences};
}

1;
