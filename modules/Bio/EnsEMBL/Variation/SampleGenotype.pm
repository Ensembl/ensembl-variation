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

# Ensembl module for Bio::EnsEMBL::Variation::SampleGenotype
#
#


=head1 NAME

Bio::EnsEMBL::Variation::SampleGenotype-Module representing the genotype
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

package Bio::EnsEMBL::Variation::SampleGenotype;

use Bio::EnsEMBL::Variation::Genotype;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Variation::Genotype);



=head2 new

  Arg [-adaptor] :
    Bio::EnsEMBL::Variation::DBSQL::SampleGenotypeAdaptor
  Arg [-genotype] :
    arrayref - arrayref of alleles making up this genotype (in haplotype order)
  Arg [-allele2] :
    string - One of the two alleles defining this genotype
  Arg [-variation] :
    Bio::EnsEMBL::Variation::Variation - The variation associated with this
    genotype
  Arg [-sample] :
    Bio::EnsEMBL::Variation::Sample - The sample this genotype is for.
  Example    : $sample_genotype = Bio:EnsEMBL::Variation::SampleGenotype->new(
				-start      => 100,
				-end        => 100,
				-strand     => 1,
				-slice      => $slice,
				-genotype   => ['A','T'],
				-variation  => $variation,
				-sample => $sample
			   );
  Description: Constructor.  Instantiates a SampleGenotype object.
  Returntype : Bio::EnsEMBL::Variation::SampleGenotype
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my ($adaptor, $genotype, $var, $varid, $ssid, $sample, $phased) =
    rearrange([qw(adaptor genotype variation _variation_id subsnp sample phased)],@_);

  if(defined($var) && (!ref($var) || !$var->isa('Bio::EnsEMBL::Variation::Variation'))) {
    throw("Bio::EnsEMBL::Variation::Variation argument expected");
  }

  if(defined($sample) && (!ref($sample) || !$sample->isa('Bio::EnsEMBL::Variation::Sample'))) {
    throw("Bio::EnsEMBL::Variation::Sample argument expected");
  }

  $self->{'adaptor'}       = $adaptor;
  $self->{'genotype'}      = $genotype;
  $self->{'sample'}        = $sample;
  $self->{'variation'}     = $var;
  $self->{'subsnp'}        = $ssid;
  $self->{'phased'}        = $phased;
  $self->{'_variation_id'} = $varid unless defined $var;

  return $self;
}


=head2 sample

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::Sample $ind
  Example    : $sample = $sample_genotype->sample();
  Description: Getter/Setter for the sample associated with this genotype
  Returntype : Bio::EnsEMBL::Variation::Sample
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut


sub sample {
  my $self = shift;
  if (@_) {
    my $sample = shift;
    if (defined($sample) && (!ref($sample) || !$sample->isa('Bio::EnsEMBL::Variation::Sample'))) {
      throw('Bio::EnsEMBL::Variation::Sample argument expected');
    }
    return $self->{'sample'} = $sample;
  }

  if (!defined($self->{sample}) && defined($self->{sample_id})) {
    my $sa = $self->adaptor->db->get_SampleAdaptor;

    if (defined($sa)) {
      my $sample = $sa->fetch_by_dbID($self->{sample_id});

      if (defined($sample)) {
        $self->{sample} = $sample;
        delete $self->{sample_id};
      }
    }
  }
  return $self->{'sample'};
}

1;
