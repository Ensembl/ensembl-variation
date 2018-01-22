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
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeFeatureAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeFeatureAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $iga = $reg->get_adaptor("human","variation","individualgenotype");

  #returns all genotypes in a certain Slice

  $genotypes = $iga->fetch_by_Slice($slice);



=head1 DESCRIPTION

This adaptor provides database connectivity for IndividualGenotypeFeature objects.
IndividualGenotypeFeatures may be retrieved from the Ensembl variation database by
several means using this module.

=head1 METHODS

=cut

package Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeFeatureAdaptor;

use strict;
use warnings;

use vars qw(@ISA);

use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::BaseGenotypeAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);
use Bio::EnsEMBL::Variation::IndividualGenotypeFeature;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor Bio::EnsEMBL::Variation::DBSQL::BaseGenotypeAdaptor);


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
  deprecate("Please use Bio::EnsEMBL::Variation::SampleGenotypeFeatureAdaptor::fetch_all_by_Variation.\n");
}

=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL:Slice $slice
  Arg [2]    : (optional) Bio::EnsEMBL::Variation::Individual $individual
  Example    : my @IndividualGenotypesFeatures = @{$ca->fetch_all_by_Slice($slice)};
  Description: Retrieves all IndividualGenotypeFeature features for a given slice for
               a certain individual (if provided). 
  Returntype : reference to list Bio::EnsEMBL::Variation::IndividualGenotype
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Slice{
  deprecate("Please use Bio::EnsEMBL::Variation::SampleGenotypeFeatureAdaptor::fetch_all_by_Slice.\n");
}

1;
