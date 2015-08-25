=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
package Bio::EnsEMBL::Variation::Pipeline::DumpVEP::QC;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::DumpVEP::BaseVEP);

use File::Path qw(make_path remove_tree);

sub run {
  my $self = shift;

  my $ensembl_cvs_root_dir = $self->param('ensembl_cvs_root_dir'); 
  my $ensembl_variation = "$ensembl_cvs_root_dir/ensembl-variation/";
  my $ensembl_tools = "$ensembl_cvs_root_dir/ensembl-tools/";
  my $dir = $self->param('pipeline_dir');
  my $cache_dir = "$dir/test_cache_dir";

  remove_tree($cache_dir) if (-d $cache_dir);
  make_path($cache_dir);

  my $version = $self->param('ensembl_release');
  
  foreach my $assembly (qw/GRCh37 GRCh38/) {
    foreach my $type (qw/vep refseq_vep merged_vep/) {
      my $file = "homo_sapiens_$type\_$version\_$assembly.tar.gz";
      if (-f "$dir/$file") { 
        $self->run_cmd("cp $dir/$file $cache_dir");
        $self->run_cmd("tar -C $cache_dir -xzf $cache_dir/$file");
      }
    } 
  }                 
  my $cmd = "perl $ensembl_variation/scripts/misc/test_chrom_coverage_in_cache_files.pl -version $version -cache_dir $cache_dir -script_dir $ensembl_tools/scripts/variant_effect_predictor/";
  my $out_file = "$dir/test_chrom_coverage.out";
  my $err_file = "$dir/test_chrom_coverage.err";
  $self->run_cmd("$cmd 1>$out_file 2>$err_file");

}

1;
