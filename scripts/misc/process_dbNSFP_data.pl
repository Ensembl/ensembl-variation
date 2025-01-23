# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2025] EMBL-European Bioinformatics Institute
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


=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.
Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

process_dbNSFP_data.pl 

=head1 DESCRIPTION

After a new dbNSFP release, this script can be used to merge, sort and compress
the chromosome level dbNSFP files into a single file. The file is then tabixed.
The script generates files for GRCh37 and GRCh38.

=head1 SYNOPSIS

process_dbNSFP_data.pl [arguments]

=head1 OPTIONS

=over 4

=item B<--help>

Displays this documentation

=item B<--dbNSFP_dir DIR>

Directory which contains the new dbNSFP files.

=item B<--dbNSFP_version>

The new dbNSFP version.

=item B<--tmp_dir DIR >

A directory which can be used by the linux sort
program for storing temporary files that a created
in the sorting process.

=back

=cut

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage qw(pod2usage);

my $args = scalar @ARGV;
my $config = {};
GetOptions(
  $config,
  'tmp_dir=s',
  'dbNSFP_version=s',
  'dbNSFP_dir=s',
);

pod2usage(1) if ($config->{'help'} || !$args);

foreach my $arg (qw/tmp_dir dbNSFP_version dbNSFP_dir/) {
  if (!$config->{$arg}) {
    die("Argument --$arg is required.");
  }
}

my $tmp_dir = $config->{tmp_dir};
my $dir = $config->{dbNSFP_dir};
my $version = $config->{dbNSFP_version};

my @assemblies = qw/grch37 grch38/;

# check chromosome file exists
if (!-e "$dir/dbNSFP$version\_variant.chr1.gz") {
  die("File $dir/dbNSFP$version\_variant.chr1.gz doesn't exist");
}

# create header file
run_system_cmd("zcat $dir/dbNSFP$version\_variant.chr1.gz | head -n1 > $dir/h");

my $cmd_args = {
  'extra' => {
    'grch37' => '| awk \'$8 != "." \'',
    'grch38' => '',
  },
  'sort' => {
    'grch37' => '-k8,8 -k9,9n',
    'grch38' => '-k1,1 -k2,2n',
  },
  tabix => {
    'grch37' => '-s 8 -b 9 -e 9',
    'grch38' => '-s 1 -b 2 -e 2',
  }
};

foreach my $assembly (@assemblies) {
  my $tabix_arg = $cmd_args->{tabix}->{$assembly};
  my $sort_arg = $cmd_args->{sort}->{$assembly};
  my $extra_arg = $cmd_args->{extra}->{$assembly};
  run_system_cmd("zgrep -h -v ^#chr $dir/dbNSFP$version\_variant.chr* $extra_arg | sort -T $tmp_dir $sort_arg - | cat $dir/h - | bgzip -c > $dir/dbNSFP$version\_$assembly.gz");
  run_system_cmd("tabix $tabix_arg $dir/dbNSFP$version\_$assembly.gz");
}

sub run_system_cmd {
  my $cmd = shift;
  print STDERR "Run CMD: $cmd\n";
  my $rc = system($cmd);
  if ($rc != 0) {
    die "Command failed: $cmd\n";
  }
}
