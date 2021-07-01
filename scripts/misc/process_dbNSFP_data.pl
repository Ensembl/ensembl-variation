# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2021] EMBL-European Bioinformatics Institute
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
=cut

use strict;
use warnings;

use Getopt::Long;

my $config = {};
GetOptions(
  $config,
  'tmp_dir|t',
  'dbNSFP_version|v',
  'dbNSFP_dir|d',
);

print $ENV{USER}, "\n";
#foreach my $arg (qw/tmp_dir dbNSFP_version dbNSFP_dir/) {
#  if (!$config->{$arg}) {
#    die("Argument --$arg is required.");
#  }
#}

my $tmp_dir = $config->{tmp_dir} || '/hps/nobackup2/production/ensembl/anja';
my $dir = $config->{dbNSFP_dir} || '/nfs/production/panda/ensembl/variation/data/dbNSFP/4.2a';
my $version = $config->{dbNSFP_version} || '4.2a';

my @assemblies = qw/grch37 grch38/;

# check chromosome file exists
if (!-e "$dir/dbNSFP$version\_variant.chr1.gz") {
  die("File $dir/dbNSFP$version\_variant.chr1.gz doesn't exist");
}

# create header file
run_system_cmd("zcat $dir/dbNSFP$version\_variant.chr1.gz | head -n1 > $dir/h");

my $args = {
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
  my $tabix_args = $args->{tabix}->{$assembly};
  my $sort_args = $args->{sort}->{$assembly};
  run_system_cmd("zgrep -h -v ^#chr $dir/dbNSFP$version\_variant.chr* | sort -T $tmp_dir $sort_args - | cat $dir/h - | bgzip -c > $dir/dbNSFP$version\_$assembly.gz");
  run_system_cmd("tabix $tabix_args $dir/dbNSFP$version\_$assembly.gz");
}

sub run_system_cmd {
  my $cmd = shift;
  print "$cmd\n"; 

}
