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
 developers list at <https://lists.ensembl.org/mailman/listinfo/dev>.
 Questions may also be sent to the Ensembl help desk at
 <https://www.ensembl.org/Help/Contact>.
=cut

package Bio::EnsEMBL::Variation::Pipeline::RemappingVCF::JoinRemapping;

use base ('Bio::EnsEMBL::Hive::Process');

use strict;
use warnings;
use FileHandle;

sub run {
  my $self = shift;
  my $population = $self->param('population');
  my $pipeline_dir = $self->param('pipeline_dir');

  my $vcf_file = "$pipeline_dir/$population.vcf";

  if (-e $vcf_file) {
    $self->run_system_command("rm $vcf_file");
  }

  my $header = "$pipeline_dir/h";

 $self->run_system_command(
    sprintf(
      'cat %s >> %s',
      $header,
      $vcf_file
    )
  );
  my @files = ();
  opendir(DIR, $pipeline_dir) or die $!;
  while (my $file = readdir(DIR)) {
    if ($file =~ /new_assembly$/) {
      push @files, $file;
    }
  }
  closedir(DIR);

  foreach my $file (@files) {
    $self->run_system_command(
      sprintf(
        'cat %s >> %s',
        "$pipeline_dir/$file",
        $vcf_file
      )
    );
  }


 $self->run_system_command("vcf-sort < $vcf_file | bgzip > $vcf_file.gz");

 $self->run_system_command("tabix $vcf_file.gz");

}

1;
