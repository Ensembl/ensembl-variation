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

=head1 NAME

Bio::EnsEMBL::Variation::Pipeline::DBSNPImport::ReportDBSNPImport

=head1 DESCRIPTION

Reports on the dbSNP import

=cut

package Bio::EnsEMBL::Variation::Pipeline::DBSNPImport::ReportDBSNPImport;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

use FileHandle;
use Bio::EnsEMBL::Registry;
use File::Path qw(make_path);
use Bio::EnsEMBL::Variation::Utils::Date;
use POSIX;

sub fetch_input {
  my $self = shift;
}

sub run {
  my $self = shift;
  my $dir = $self->param_required('pipeline_dir');
  my $file = 'ReportDBSNPImport_report.txt';
  $self->warning("writing to $dir/$file");
  open(my $report, ">", "$dir/$file") or die("Failed to open $dir/$file.txt : $!");
  my $run_time = strftime("%Y-%m-%d_%H%M%S", localtime);
  print $report "Completed run ", $run_time, "\n";
  close($report);
}

1;
