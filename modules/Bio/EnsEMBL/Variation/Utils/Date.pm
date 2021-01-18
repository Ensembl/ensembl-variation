


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

=head1 NAME

Bio::EnsEMBL::Variation::Utils::Date

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Variation::Utils::Date;

use strict;
use warnings;
use base qw(Exporter);
use POSIX;

our @EXPORT_OK = qw(run_date log_time);


=head2 run_date


  Example     : my $run_date = run_date();
  Description : returns the current date in a standard format to use in the meta
                table as a process run date
  ReturnType  : String
  Exceptions  : None
  Caller      : General


=cut

sub run_date{

    return strftime("%Y-%m-%d", localtime);
}

=head2 log_time


  Example     : my $filename = $process_name . log_time() . 'txt';
  Description : returns the current date/time in a standard format to use in
                output files
  ReturnType  : String
  Exceptions  : None
  Caller      : General


=cut

sub log_time{

    return strftime("%Y-%m-%d_%H%M%S", localtime);
}


1;
