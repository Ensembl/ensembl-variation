=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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
package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::RunCADD;

use strict;
use Bio::EnsEMBL::Variation::Utils::CADDProteinFunctionAnnotation;
use File::Path qw(make_path);

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

sub run {
  my $self = shift;
  my $working_dir = $self->required_param('cadd_working');
  die "Working directory ($working_dir) doesn't exist" unless (-d $working_dir);

  my $assembly = $self->get_assembly();
  my $cadd_annotation = $self->required_param('cadd_annotation');
  die "Assembly ($assembly) is not supported" if (!defined $cadd_annotation->{$assembly});
  my $annotation_file = $cadd_annotation->{$assembly}->{file};
  my $annotation_file_version = $cadd_annotation->{$assembly}->{version};

  my $cadd = Bio::EnsEMBL::Variation::Utils::CADDProteinFunctionAnnotation->new(
    -registry_file => $self->required_param('ensembl_registry'),
    -species => $self->required_param('species'),
    -working_dir => $working_dir,
    -annotation_file =>  $annotation_file,
    -assembly => $assembly,
    -annotation_file_version => $annotation_file_version,
  );

  my $translation_md5 = $self->required_param('translation_md5');
  $cadd->run($translation_md5);

}

1;
