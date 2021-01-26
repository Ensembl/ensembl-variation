=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     https://www.apache.org/licenses/LICENSE-2.0

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


=head1 FinishPhenotypeAnnotation

This module runs at the end of the phenotype annotation import pipeline and produces
a summary report of the run times based on the hive db.

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::FinishPhenotypeAnnotationPipeline;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation);

sub run {
  my $self = shift;

  my $hive_dba = $self->dbc;

  my $runTime_sth = $hive_dba->prepare(qq[ select timediff(max(when_died), min(when_born)) from worker ]);
  $runTime_sth->execute() || die ("Failed to fetch timediff from hive db: $!\n");
  my $time = $runTime_sth->fetchall_arrayref();

  my $runTime_imports_sth = $hive_dba->prepare(qq[
            SELECT ab.logic_name, timediff(max(when_finished), min(when_started))
            FROM analysis_base ab, role r
            WHERE ab.logic_name like 'import_%' AND ab.analysis_id = r.analysis_id
            GROUP BY ab.logic_name ]);
  $runTime_imports_sth->execute() || die ("Failed to fetch runtime stats from hive db: $!\n");

  my $dir =$self->required_param('pipeline_dir');
  open(my $report, ">$dir/REPORT_hive_pipe.txt") || die ("Failed to open report file for summary info: $!\n");

  print $report "PhenotypeAnnotation pipeline finished! \n";
  print $report "running time: $time->[0]->[0] \n";
  print $report "pipeline_dir: ", $self->required_param('pipeline_dir'), "\n";

  while(my @row = $runTime_imports_sth->fetchrow_array()) {
    print $report join("\t", @row)."\n";
  }

  close $report;
}

1;

