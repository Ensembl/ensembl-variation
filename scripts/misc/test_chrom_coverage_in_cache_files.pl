#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


## script to create a table of HGVS stings and variation_ids for search index building


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

use strict;
use warnings;

use FileHandle;
use Getopt::Long;

my $config = {};
GetOptions(
  $config,
  'version=i',
  'script_dir=s',
  'vep_input=s',
  'cache_dir=s',
  'no_cleanup',
  'help!',
) or die "Error: Failed to parse command line arguments\n";

if ($config->{help}) {
  usage();
  exit(0);
}

if (!$config->{version}) {
  die "Version is missing. Set version with --version.";
}

my $version = $config->{version};

my @assemblies = ('GRCh37', 'GRCh38');

my $cache_dirs = {
  'homo_sapiens' => '',
  'homo_sapiens_merged' => '--merged',
  'homo_sapiens_refseq' => '--refseq',
};

foreach my $dir_name (qw/script_dir cache_dir/) {
  if (! $config->{$dir_name}) {
    die "Parameter ($dir_name) is not set.";
  }
  if (!(-d $config->{$dir_name})) {
    die "$dir_name: ", $config->{$dir_name}, " is not a directory.";
  }
}

my $script_dir = $config->{script_dir};
die "script_dir is missing file variant_effect_predictor.pl" unless (-f "$script_dir/variant_effect_predictor.pl");

$config->{vep_input} ||= "$script_dir/t/testdata/test_vep_input.txt.gz";
die "vep_input: ", $config->{vep_input}, " is not a file." unless (-f $config->{vep_input});

$config->{cache_dir} ||= $ENV{HOME} . '/.vep';

#all_cache_files_are_installed($config);
run_vep($config);
tests($config);
cleanup($config) unless($config->{no_cleanup});

sub all_cache_files_are_installed {
  my $config = shift;
  my $root_cache_dir = $config->{cache_dir};

  foreach my $cache_dir (keys %$cache_dirs) {
    foreach my $assembly (@assemblies) {
      my $dir = "$root_cache_dir/$cache_dir/$version\_$assembly";
      if (!(-d $dir)) {
        die "$dir is not a directory. Cache files for $cache_dir/$version\_$assembly are missing.";
      }
    }
  }
  return 1;
}

sub run_vep {
  my $config = shift;
  my $script_dir = $config->{script_dir};
  my $input = $config->{vep_input};
  my $root_cache_dir = $config->{cache_dir};
  my $version = $config->{version};

  my $output_files = {};

  foreach my $cache_dir (keys %$cache_dirs) {
    foreach my $assembly (@assemblies) {
      my $dir = "$root_cache_dir/$cache_dir/$version\_$assembly";
      if (-d $dir) {
        my $vep_run_name = "$cache_dir\_$version\_$assembly";
        my $params = $cache_dirs->{$cache_dir} . " --assembly $assembly";
        my $output = "$script_dir/test_vep_output_$vep_run_name"; 
        my $err_file = "$script_dir/err_$vep_run_name"; 
        my $out_file = "$script_dir/out_$vep_run_name"; 
        
        $output_files->{$vep_run_name}->{vep_output} = $output;
        $output_files->{$vep_run_name}->{err} = $err_file;
        $output_files->{$vep_run_name}->{out} = $out_file;
        # -cache_version
        my $cmd = "perl $script_dir/variant_effect_predictor.pl --cache_version $version --cache --offline --dir $root_cache_dir -i $input -o $output --force_overwrite --no_stats --regulatory --sift b --polyphen b $params";

        run_cmd("$cmd 1>$out_file 2>$err_file");
      }
    }
  } 

  $config->{output_files} = $output_files;
}

sub tests {
  my $config = shift;
  my $output_files = $config->{output_files};

  my @annotations = qw/SIFT PolyPhen MotifFeature RegulatoryFeature/;

  foreach my $vep_run_name (keys %$output_files) {
    my $vep_out_file = $output_files->{$vep_run_name}->{vep_output};
    my $fh = FileHandle->new($vep_out_file, 'r');
    my $covered_chroms = {};
    my $has_annotation = {};
    while (<$fh>) {
      chomp;
      next if /^#/;
      my ($name, $location, $rest) = split("\t", $_, 3);
      my ($chrom, $position) = split(':', $location);
      $covered_chroms->{$chrom} = 1;
      foreach my $annotation (@annotations) {
        if ($rest =~ /$annotation/) {
          $has_annotation->{$annotation} = 1;
        } 
      }
    }
    $fh->close();
    foreach my $chrom (1..22, 'X', 'Y', 'MT') {
      if (!$covered_chroms->{$chrom}) {
        die "Chrom $chrom is missing from VEP output $vep_out_file. Need to check cache files are dumped correctly.";
      }
    }
    print STDOUT "All chromosomes are covered in $vep_out_file\n";
    foreach my $annotation (@annotations) {
      if (!$has_annotation->{$annotation}) {
        die "Annotation: $annotation is missing from VEP output $vep_out_file. Need to check cache files are dumped correctly.";
      }
    }
    print STDOUT "Annotations (", join(', ', @annotations), ") are contained in $vep_out_file\n";
  }
  return 1;
}

sub run_cmd {
  my $cmd = shift;
  if (my $return_value = system($cmd)) {
    $return_value >>= 8;
    die "system($cmd) failed: $return_value";
  }
}

sub cleanup {
  my $config = shift;
  my $output_files = $config->{output_files};

  foreach my $vep_run_name (keys %$output_files) {
    foreach my $file_type (qw/vep_output err out/) {
      my $file = $output_files->{$vep_run_name}->{$file_type};
      run_cmd("rm $file");
    }
  }
}

sub usage {
    my $usage =<<END;

Usage:
perl test_chrom_coverage_in_cache_files.pl [arguments]

The script runs for human only. It checks that all human chromosomes have been dumped to the
cache files (default, refseq, merged) for GRCh37 and GRCh38.

The script should be run after the cache file generation. Copy and unpack all the human cache files
to the cache file directory (--cache_dir).

bsub -J test_chrom_coverage -o out -e err -R"select[mem>2500] rusage[mem=2500]" -M2500 perl test_chrom_coverage_in_cache_files.pl -version 78 -cache_dir /lustre/scratch110/ensembl/at7/vep/ -script_dir ~/DEV/ensembl-tools/scripts/variant_effect_predictor/

Options
=======
--help       Display this message and quit
--version    Set the version for the new release
--script_dir Location of variant_effect_predictor.pl script
--cache_dir  Cache file directory
--no_cleanup Don't clean up err, out, vep_output files

END

    print $usage;
}
