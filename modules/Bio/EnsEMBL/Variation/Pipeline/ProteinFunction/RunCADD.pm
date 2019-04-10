=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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
use Data::Dumper;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

sub run {
  my $self = shift;
  my $working_dir = $self->param('cadd_working');
  unless (-d $working_dir) {
    my $err;
    make_path($working_dir, {error => \$err});
    die "make_path failed: ".Dumper($err) if $err && @$err;
  }

  my $assembly = $self->get_assembly();
  my $cadd_annotation = $self->param('cadd_annotation');
  my $annotation_file = $cadd_annotation->{$assembly}->{file};
  my $annotation_file_version = $cadd_annotation->{$assembly}->{version};

  my $cadd = Bio::EnsEMBL::Variation::Utils::CADDProteinFunctionAnnotation->new(
    -registry_file => $self->param('ensembl_registry'),
    -species => $self->param('species'),
    -working_dir => $working_dir,
    -annotation_file =>  $annotation_file,
    -assembly => $assembly,
    -annotation_file_version => $annotation_file_version,
  );

  my $translation_md5 = $self->param('translation_md5');
  $cadd->run($translation_md5);

}

1;
