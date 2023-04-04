=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2023] EMBL-European Bioinformatics Institute

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

=head1 NAME

Bio::EnsEMBL::Variation::Utils::Reports

=head1 DESCRIPTION

Create files with reports

=cut

package Bio::EnsEMBL::Variation::Utils::Reports;

use strict;
use warnings;
use base qw(Exporter);
use Bio::EnsEMBL::Variation::Utils::QCUtils qw(count_rows);
use Bio::EnsEMBL::Variation::Utils::Date qw(log_time);


our @EXPORT_OK = qw(report_counts);


=head2 report_counts

  Example     : report_counts($dba, $flag, $tables);
  Description : writes to file the number of rows for each table before/after an import
  ReturnType  : None
  Exceptions  : None
  Caller      : General


=cut

sub report_counts {
  my $dba    = shift;
  my $flag   = shift;
  my $tables = shift;
  my $output_file = shift;

  die ("ERROR: cannot report the number of variation counts\n") unless ($dba && $flag && $tables);

  my $file_counts = $output_file ? $output_file : "report_variation_counts.txt";

  # Open file to write
  my $type;
  my $message;
  if($flag eq 'before') {
    $type = ">";
    $message = "(" . log_time() . ") Counts before import:";
  }
  else {
    $type = ">>";
    $message = "\n" . "(" . log_time() . ") Counts after import:";
  }

  open(my $fh, $type, $file_counts) or die "Cannot write to file: $!\n";
  print $fh $message . "\n";

  foreach my $table (@{$tables}) {
    my $count = count_rows($dba, $table);
    print $fh $table . "\t" . $count . "\n";
  }

  close($fh);
}

1;
