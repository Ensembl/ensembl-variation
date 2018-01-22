# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

use strict;
use warnings;
#use diagnostics;

use SubGVFDumper;
use IndividualGVFDumper;
use PopulationGVFDumper;
use Bio::EnsEMBL::Registry;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Variation::Utils::EnsEMBL2GFF3;

my $species;
my $output_file;
my $registry;
my $host;
my $user;
my $port;
my $compress;
my $chunk_size;
my @seq_regions;

my $somatic;
my $use_iterator;
my $include_failed;
my $just_failed;
my $set_name;
my $include_consequences;
my $include_coding_details;
my $include_global_maf;
my $include_svs;
my $just_svs;

my $individual;
my $population;

my $load_balance;

my $help;
my $log_file;

GetOptions(
    "species|s=s"                   => \$species,
    "output|o=s"                    => \$output_file,
    "registry|r=s"                  => \$registry,
    "host=s"                        => \$host,
    "user=s"                        => \$user,
    "individual|i=s"                => \$individual,
    "population|p=s"                => \$population,
    "set=s"                         => \$set_name,
    "somatic"                       => \$somatic,
    "compress|c"                    => \$compress,
    "use_iterator"                  => \$use_iterator,
    "include_failed"                => \$include_failed,
    "just_failed"                   => \$just_failed,
    "chunk_size|cs=i"               => \$chunk_size,
    "include_consequences"          => \$include_consequences,
    "include_coding_details"        => \$include_coding_details,
    "include_global_maf"            => \$include_global_maf,
    "seq_regions|sr=s{,}"           => \@seq_regions,
    "include_structural_variations" => \$include_svs,
    "just_structural_variations"    => \$just_svs,
    "load_balance"                => \$load_balance, 
    "help|h"                        => \$help,
    "log_file=s"                    => \$log_file,
) or die pod2usage(1);

pod2usage(1) if $help;

die "species argument required, try --help for usage instructions\n" unless $species;

die "Can't fetch for a population and an individual at once" 
    if $population && $individual;

# chunk size is in kilobases
$chunk_size *= 100 if defined $chunk_size;

# if the user wants just structural variants that implies that SVs should be included
$include_svs = 1 if $just_svs;

# default to a sensible file name
$output_file ||= "$species.gvf";

my $reg = 'Bio::EnsEMBL::Registry';

if ($host) {
    # if we are supplied with a host, try to load the registry from there
    $user ||= 'anonymous';
    $reg->load_registry_from_db(-host => $host, -user => $user, -port => $port);
}
else {
    # otherwise use the registry file supplied (or default to the 
    # ENSEMBL_REGISTRY environment variable)
    $reg->load_all($registry);
}

my $translate_to_seq_region = {
    1  => [1, 2],
    2  => [3, 4],
    3  => [5, 6],
    4  => [7, 8],
    5  => [9, 10],
    6  => [11, 12],
    7  => [13, 14],
    8  => [15, 16],
    9  => [17, 18],
    10 => [19, 20],
    11 => [21, 22],
    12 => ['X', 'Y', 'MT'],
};

if ($load_balance) {
    $load_balance = $ENV{'LSB_JOBINDEX'};
    $output_file = $output_file . '.' . $ENV{'LSB_JOBINDEX'} . '.gvf';
    push @seq_regions, @{$translate_to_seq_region->{$load_balance}};
} else {
    $output_file = $output_file . '.gvf';
}

my $cdba = $reg->get_DBAdaptor($species, 'core') 
    or die "Failed to get core DBAdaptor";

my $vdba = $reg->get_DBAdaptor($species, 'variation') 
    or die "Failed to get variation DBAdaptor";

my %arguments = (
    'cdba'                   => $cdba,
    'vdba'                   => $vdba,
    'use_iterator'           => $use_iterator,
    'somatic'                => $somatic,
    'set_name'               => $set_name,

    'include_failed'            => $include_failed,
    'just_failed'               => $just_failed,
    'include_consequences'      => $include_consequences,
    'include_coding_details'    => $include_coding_details,
    'include_global_maf'        => 1,
    'include_validation_states' => 1,
    'include_clinical_significance' => 1,
    'include_svs'               => $include_svs,
    'just_svs'                  => $just_svs,

    'chunk_size'             => $chunk_size,
    'seq_region_names'       => \@seq_regions,
    
    'individual'             => $individual,
    'population'             => $population,

    'load_balance'           => $load_balance,

    'cache_dir'              => '/Path_to/cache_files_70/',
    'output'                 => $output_file,
);

my $dumper;
if ($individual) {
    $dumper = IndividualGVFDumper->new(%arguments);
} elsif ($population) {
    $dumper = PopulationGVFDumper->new(%arguments);
} else {
    $dumper = GVFDumper->new(%arguments);
}

$dumper->dump();

