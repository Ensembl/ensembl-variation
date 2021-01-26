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

Bio::EnsEMBL::Variation::Pipeline::DBSNPImport::InitDBSNPImport

=head1 DESCRIPTION

Initialises the dbSNP import

=cut

package Bio::EnsEMBL::Variation::Pipeline::DBSNPImport::InitDBSNPImport;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

use FileHandle;
use Bio::EnsEMBL::Registry;
use File::Path qw(make_path);
use Bio::EnsEMBL::Variation::Utils::Date;
use POSIX;

my @chrs = (1..22, 'X', 'Y', 'MT');

sub fetch_input {
  my $self = shift;

  my $data_dir = $self->param_required('data_dir');

}

sub run {
  my $self = shift;
  
  my $data_dir = $self->param_required('data_dir');
  my $rpt_dir = $self->param_required('rpt_dir');

  if (! -d $data_dir) {
    die("No data directory ($data_dir)");
  }

  if (! -d $rpt_dir) {
    die("No rpt directory ($rpt_dir)");
  }

  my @sub_dirs = map('chr' . $_, @chrs);
  for my $sub_dir (@sub_dirs) {
    if (! -d "$data_dir/$sub_dir") {
      die("No data directory for ${data_dir}/${sub_dir}");
    }
    if (! -d "$rpt_dir/$sub_dir") {
      die("No rpt directory for ${rpt_dir}/${sub_dir}");
    }
  }
  # set up the list of sub_dir
  $self->param('sub_dirs', [ map { {sub_dir => $_} } @sub_dirs]);

  # Check the assembly is valid
  my $assembly = $self->param_required('assembly');
  if ($assembly !~ /^(GRCh37|GRCh38)/) {
      die("Assembly ($assembly) is invalid. Please specify GRCh37 or GRCh38");
  }
}

sub write_output {
  my $self = shift @_;

  $self->dataflow_output_id($self->param('sub_dirs'), 2);
}


1;
