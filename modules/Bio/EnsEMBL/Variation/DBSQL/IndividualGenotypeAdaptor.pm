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
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeAdaptor

=head1 SYNOPSIS

  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_registry_from_db(-host => 'ensembldb.ensembl.org', -user => 'anonymous');

  my $sa = $registry->get_adaptor('human', 'core', 'slice');
  my $igta = $registry->get_adaptor('human', 'variation', 'individualgenotype');

  # Fetch region for which we want to get all individual genotypes
  my $slice = $sa->fetch_by_region('chromosome', '3', 52_786_960, 52_786_970);
  my $individual_genotypes = $igta->fetch_all_by_Slice($slice);

  foreach my $igt (@$individual_genotypes) {
    my $variation_name = $igt->variation()->name;
    my $genotype = $igt->genotype_string;
    my $individual_name = $igt->individual()->name;
    print "$variation_name\t$genotype\t$individual_name\n";
  }

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
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);
use Scalar::Util qw(weaken);

our $CACHE_SIZE = 5;

@ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseGenotypeAdaptor');

# stores a listref of individual genotype objects
sub store {
  deprecate("Please use Bio::EnsEMBL::Variation::SampleGenotypeAdaptor::store.\n");
}

# similar to store above but writes to old-style non-compressed table
# defaults to individual_genotype_multiple_bp since
# tmp_individual_genotype_single_bp will only accept single nucleotide genotypes
sub store_uncompressed {
  deprecate("Please use Bio::EnsEMBL::Variation::SampleGenotypeAdaptor::store_uncompressed.\n");
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
  deprecate("Please use Bio::EnsEMBL::Variation::SampleGenotypeAdaptor::fetch_all_by_Variation.\n");
}

sub fetch_all_by_Variation_dbID {
  deprecate("Please use Bio::EnsEMBL::Variation::SampleGenotypeAdaptor::fetch_all_by_Variation_dbID.\n");
}

sub fetch_all_by_Slice {
  deprecate("Please use Bio::EnsEMBL::Variation::SampleGenotypeAdaptor::fetch_all_by_Slice.\n");
}

1;
