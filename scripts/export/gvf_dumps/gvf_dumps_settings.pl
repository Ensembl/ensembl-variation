#!/usr/bin/env perl
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


# release_gvf.pl - Create GVF dumps for an ensembl release
#
# see the end of the file for documentation

use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Path qw(make_path);
use FindBin qw($Bin);

my $registry;
my $toplevel_dir;
my $script;
my $test;
my @species;
my @ignore;
my $help;
my $opts;

GetOptions(
    "registry|r=s"      => \$registry,
    "output_dir|o=s"    => \$toplevel_dir,
    "script|s=s"        => \$script, 
    "test|t"            => \$test,
    "species|s=s"       => \@species,
    "ignore|i=s"        => \@ignore,
    "help|h"            => \$help,
) or die pod2usage(1);

pod2usage(1) if $help;

die "--output_dir argument is required, try --help for usage" 
    unless $toplevel_dir;

die "--registry argument is required, try --help for usage" 
    unless $registry;

$script ||= "$Bin/dump_gvf.pl";

my $default_rc = '-q normal -R"select[mem>1500] rusage[mem=1500]" -M1500000';
my %to_ignore = map { $_ => 1 } @ignore;

# differenct categories of data in variation database
my $data = {
    all => {
        postfix => '',
        argument => '',
        all => 1,
    },
    include_consequence => {
        postfix  => '_incl_consequences',
        argument => '--include_consequences',
        all => 1,
    },
    failed => {
        postfix  => '_failed',
        argument => '--just_failed',
        all => 1,
    },
    structural_variations => {
        postfix  => '_structural_variations',
        argument => '--just_structural_variations',
        all => 0,
    },
    somatic => {
        postfix => '_somatic',
        argument => '--somatic',
        all => 0,
    },
    somatic_include_consequence => {
        postfix => '_somatic_incl_consequences',
        argument => '--somatic --include_consequences',
        all => 0,
    },
    individuals => {
        postfix => '_',
        argument => '',
        all => 0, 
    },
    populations => {
        postfix => '_',
        argument => '',
        all => 0,
    }
};

my $species_config = {
    'Homo_sapiens' => {
        settings  => {
            load_balance => {
                populations => 1,
                include_consequence => 1,
            },
            rc => '-q normal -R"select[mem>4000] rusage[mem=4000]" -M4000000',
        },
        #individuals     => {
        #    Watson  => 'ENSEMBL:Watson',
        #    Venter  => 'ENSEMBL:Venter',
        #},
        populations     => {
            '1000GENOMES:phase_1_AFR' => '',
            '1000GENOMES:phase_1_AMR' => '',
            '1000GENOMES:phase_1_ASN' => '',
            '1000GENOMES:phase_1_EUR' => '',
            #'1000GENOMES' => '',

            #'CSHL-HAPMAP:HapMap-CEU'  => '',
            #'CSHL-HAPMAP:HapMap-HCB'  => '',
            #'CSHL-HAPMAP:HapMap-JPT'  => '',
            #'CSHL-HAPMAP:HapMap-YRI'  => '',
            #'CSHL-HAPMAP:HAPMAP-ASW'  => '',
            #'CSHL-HAPMAP:HAPMAP-CHB'  => '',
            #'CSHL-HAPMAP:HAPMAP-CHD'  => '',
            #'CSHL-HAPMAP:HAPMAP-GIH'  => '',
            #'CSHL-HAPMAP:HAPMAP-LWK'  => '',
            #'CSHL-HAPMAP:HAPMAP-MEX'  => '',
            #'CSHL-HAPMAP:HAPMAP-MKK'  => '',
            #'CSHL-HAPMAP:HAPMAP-TSI'  => '',
        },
        #somatic                 => 1,
        #structural_variations   => 1,
        #somatic_include_consequence => 1,
    },
    'Mus_musculus'	            => {structural_variations => 1},
    'Macaca_mulatta'	        => {structural_variations => 1},
    'Danio_rerio'	            => {structural_variations => 1},
    'Felis_catus'	            => {},
    'Bos_taurus'	            => {structural_variations => 1},
    'Canis_familiaris'	        => {},
    'Equus_caballus'	        => {structural_variations => 1},
    'Taeniopygia_guttata'	    => {},
    'Tetraodon_nigroviridis'    => {},
    'Monodelphis_domestica'	    => {},
    'Sus_scrofa'	            => {},
    'Ornithorhynchus_anatinus'	=> {},
    'Pan_troglodytes'	        => {},
    'Pongo_abelii'	            => {},
    'Rattus_norvegicus'	        => {},
    'Gallus_gallus'	            => {},
    'Drosophila_melanogaster'	=> {},
    'Saccharomyces_cerevisiae'	=> {},
};

unless (-e $toplevel_dir || $test) {
    make_path($toplevel_dir) or die "Failed to create toplevel output directory '$toplevel_dir': $!";
}

my $normal_opts = "--species %s --registry %s --compress --chunk_size 1000 ";

my @cmds;

@species = keys %$species_config unless @species;

for my $species_name (@species) {
    next if $to_ignore{$species_name};
    $opts = sprintf $normal_opts, $species_name, $registry;
    my $config = $species_config->{$species_name};
    
    my $dir = $toplevel_dir.'/'.lc($species_name);
    #unless (-e $dir || $test) {
    unless (-e $dir) {
        make_path($dir) or die "Failed to create output dir '$dir': $!";
    }

    foreach my $type (keys %$data) {
        if  ($data->{$type}->{all} || $config->{$type}) {
            my $load_balance = $config->{settings}->{load_balance}->{$type} || 0;
            my $rc = $config->{settings}->{rc} || $default_rc;

            my $name = $species_name . $data->{$type}->{postfix};

            if ($type eq 'individuals') {
                for my $individual (keys %{ $config->{individuals} || {} }) {
                    make_cmd($individual, "--individual '$individual'", $dir, $rc, $load_balance);
                }
            } elsif ($type eq 'populations') {
                for my $population (keys %{ $config->{populations} || {} }) {
                    my $name = $population;
                    $name =~ s/:/-/g;
                    make_cmd($name, "--population '$population'", $dir, $rc, $load_balance);
                }
            } else {
                make_cmd($name, $data->{$type}->{argument}, $dir, $rc, $load_balance);
            }
        }
    }
}

sub make_cmd {
    my ($name, $extra_opts, $dir, $rc, $load_balance) = @_;
    my ($cmd, $gvf, $out, $err);

    if ($load_balance) {
        my $job_arg = ' --load_balance ';
        $extra_opts .= $job_arg;
        $gvf = "$dir/$name";
        $out = "$dir/$name.%I.out";
        $err = "$dir/$name.%I.err";
        $name = '"' . $name . '[1-12]%4"'; 
    } else { 
        $gvf = "$dir/$name";
        $out = "$dir/$name.out";
        $err = "$dir/$name.err";
    }

    $cmd = "bsub -o $out -e $err -J $name $rc perl $script $opts $extra_opts --output $gvf";
    if ($test) {
        print "DEBUG ", $cmd, "\n\n";
    } else {
        print $cmd, "\n\n";
    }
    push @cmds, $cmd;
}

__END__
perl release_gvf.pl --output_dir /lustre/scratch110/ensembl/at7/Dumps/release_69_gvf/ --registry /lustre/scratch110/ensembl/at7/Dumps/ensembl.registry.69 --species Danio_rerio
=head1 NAME

release_gvf.pl

=head1 DESCRIPTION

Create GVF dumps for an ensembl release, submitting a farm job for each
file to be created (possibly several per species, especially for human).

You can control which individuals and populations etc. to create files for
by editing the $species_config hash at the top of this script.

=head1 SYNOPSIS

release_gvf.pl --output_dir DIR [options]

=head1 EXAMPLE COMMAND LINES

  release_gvf.pl --output_dir /lustre/scratch103/ensembl/gr5/release_61_GVFs \
    --registry ensembl.registry
    
    Create GVF dumps for all species in the specified directory, using the 
    given registry file

  release_gvf.pl --output_dir human_and_mouse_GVFs --species Mus_musculus --species Homo_sapiens

    Just create GVF dumps for human and mouse

  release_gvf.pl --output_dir non_human_GVFs --ignore Homo_sapiens

    Create dumps for all species except human
  
  release_gvf.pl --output_dir human_GVFs --species Homo_sapiens --load_file load_balance.txt

  	Create dumps for subset of regions

=head1 OPTIONS

=over 4

=item B<--registry FILE>

Read database connection details from this registry file

=item B<--output_dir DIR>

Create the files under this directory, a sub-directory will be created for each 
species. The directory will be created if it doesn't exist. This option is required.

=item B<--script NAME>

Use this script to actually create the GVF files, defaults to dump_gvf.pl in the same
directory as this script

=item B<--test>

Don't actually submit the farm jobs, just show the bsub command generated for each file

=item B<--species NAME>

Only create files for the specified species, defaults to all in the $species_config
hash at the top of this script. You can specify multiple species by supplying multiple
options for each.

=item B<--ignore NAME>

Don't create files for the specified species, you can specify multiple species to ignore.

=item B<--help>

Display this documentation

=head1

For help with this script address questions to http://lists.ensembl.org/mailman/listinfo/dev
