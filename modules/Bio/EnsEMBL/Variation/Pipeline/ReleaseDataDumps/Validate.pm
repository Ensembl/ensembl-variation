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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::Validate;

use strict;
use base ('Bio::EnsEMBL::Hive::Process');
use File::Basename;

sub run {
  my $self = shift;
  my $file_type = $self->param('file_type');
  if ($file_type eq 'gvf') {
    $self->validate_gvf();
  } elsif ($file_type eq 'vcf') {
    $self->validate_vcf();
  } else {
    die "File type ($file_type) not correct. It should be gvf of vcf";
  }
}

sub validate_gvf {
  my $self = shift;
  my $gvf_validator = $self->param('gvf_validator');
  my $so_file       = $self->param('so_file');

  my $working_dir   = $self->param('working_dir');
  my $file_name     = $self->param('file_name');	

  my $err = "$working_dir/Validate\_$file_name.err";
  my $out = "$working_dir/Validate\_$file_name.out";

  my $file = "$working_dir/$file_name.gvf";
  my $file_for_validation = "$working_dir/$file_name\_validate.gvf";
  $self->run_cmd("head -2500 $file > $file_for_validation");
  my $cmd = "$gvf_validator --so_file $so_file $file_for_validation";
  $self->run_cmd("$cmd 1>$out 2>$err");	
  $self->run_cmd("rm $file_for_validation");
}

sub validate_vcf {
  my $self = shift;
  my $vcf_validator = $self->param('vcf_validator');
  my $vcf_sort = $self->param('vcf_sort');
  my $vcf_file = $self->param('vcf_file');
  $vcf_file =~ s/--vcf_file //;
  my ($file_name, $working_dir, $suffix) = fileparse($vcf_file, qr/\.[^.]*/);

  my $err = "$working_dir/Validate\_vcf\_$file_name.err";
  my $out = "$working_dir/Validate\_vcf\_$file_name.out";

  # create short version of file for validation   
  my $file_for_validation = "$working_dir/$file_name\_validate.vcf";
  $self->run_cmd("head -2500 $vcf_file > $file_for_validation");

  # sort and bgzip
  my $cmd = "$vcf_sort < $vcf_file | bgzip > $vcf_file.gz";
  $self->run_cmd($cmd);
  $cmd = "$vcf_sort < $file_for_validation | bgzip > $file_for_validation.gz";
  $self->run_cmd($cmd);

  # validate
  $cmd = "$vcf_validator $file_for_validation.gz";
  $self->run_cmd("$cmd 1>$out 2>$err");
  $self->run_cmd("rm $file_for_validation");
}

sub run_cmd {
  my $self = shift;
  my $cmd = shift;
  if (my $return_value = system($cmd)) {
    $return_value >>=8;
    die "system($cmd) failed: $return_value";
  }
}

1;
