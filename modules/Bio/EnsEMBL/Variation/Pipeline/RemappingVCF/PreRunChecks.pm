=head1 LICENSE
Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute
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
package Bio::EnsEMBL::Variation::Pipeline::RemappingVCF::PreRunChecks;

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use IPC::Cmd qw(can_run);

use base qw(Bio::EnsEMBL::Hive::Process);

sub fetch_input {
  my $self = shift;
}

sub run {
  my $self = shift;

  my $pipeline_dir = $self->param('pipeline_dir');
  die "$pipeline_dir doesn't exist" unless (-d $pipeline_dir);		

  foreach (qw//) {
    my $file = $self->param($_);
    die "File $_ doesn't exist: $file" unless (-e $file);		
  }

  can_run('tabix') or die "Cannot run tabix";

}

1;
