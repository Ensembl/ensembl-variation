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



=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut


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

my $default_rc = '-q long';

my %to_ignore = map { $_ => 1 } @ignore;

my $species_config = {
    'Homo_sapiens'              => {
        individuals     => {
            Watson  => 'ENSEMBL:Watson',
            Venter  => 'ENSEMBL:Venter',
        },
        populations     => {
            '1000GENOMES:pilot_1_CEU_low_coverage_panel'    => '1000 genomes - Low coverage - CEU',
            '1000GENOMES:pilot_1_CHB+JPT_low_coverage_panel'=> '1000 genomes - Low coverage - CHB+JPT',
            '1000GENOMES:pilot_1_YRI_low_coverage_panel'    => '1000 genomes - Low coverage - YRI',
            '1000GENOMES:pilot_3_CEU_exon_capture_panel'    => '1000 genomes - High coverage exons - CEU',
            '1000GENOMES:pilot_3_CHB_exon_capture_panel'    => '1000 genomes - High coverage exons - CHB',
            '1000GENOMES:pilot_3_CHD_exon_capture_panel'    => '1000 genomes - High coverage exons - CHD',
            '1000GENOMES:pilot_3_JPT_exon_capture_panel'    => '1000 genomes - High coverage exons - JPT',
            '1000GENOMES:pilot_3_LWK_exon_capture_panel'    => '1000 genomes - High coverage exons - LWK',
            '1000GENOMES:pilot_3_TSI_exon_capture_panel'    => '1000 genomes - High coverage exons - TSI',
            '1000GENOMES:pilot_3_YRI_exon_capture_panel'    => '1000 genomes - High coverage exons - YRI',
            'CSHL-HAPMAP:HapMap-CEU'            => '',
            'CSHL-HAPMAP:HapMap-HCB'            => '',
            'CSHL-HAPMAP:HapMap-JPT'            => '',
            'CSHL-HAPMAP:HapMap-YRI'            => '',
            'CSHL-HAPMAP:HAPMAP-ASW'            => '',
            'CSHL-HAPMAP:HAPMAP-CHB'            => '',
            'CSHL-HAPMAP:HAPMAP-CHD'            => '',
            'CSHL-HAPMAP:HAPMAP-GIH'            => '',
            'CSHL-HAPMAP:HAPMAP-LWK'            => '',
            'CSHL-HAPMAP:HAPMAP-MEX'            => '',
            'CSHL-HAPMAP:HAPMAP-MKK'            => '',
            'CSHL-HAPMAP:HAPMAP-TSI'            => '',
        },
        rc                      => '-q long -R"select[mem>15000] rusage[mem=15000]" -M15000000',
        somatic                 => 1,
        structural_variations   => 1,
    },
    'Mus_musculus'	            => {structural_variations => 1},
    'Danio_rerio'	            => {},
    'Felis_catus'	            => {},
    'Bos_taurus'	            => {},
    'Canis_familiaris'	        => {structural_variations => 1},
    'Equus_caballus'	        => {},
    'Taeniopygia_guttata'	    => {},
    'Tetraodon_nigroviridis'    => {},
    'Monodelphis_domestica'	    => {},
    'Sus_scrofa'	            => {structural_variations => 1},
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

my $normal_opts = "--species %s --registry %s --compress --chunk_size 100";

my @cmds;

@species = keys %$species_config unless @species;

for my $species_name (@species) {

    next if $to_ignore{$species_name};
    
    my $config = $species_config->{$species_name};
    my $rc = $config->{rc} || $default_rc;

    # the webteam want the toplevel species directory to be lower case
    my $dir = $toplevel_dir.'/'.lc($species_name);
    
    unless (-e $dir || $test) {
        make_path($dir) or die "Failed to create output dir '$dir': $!";
    }

    my $opts = sprintf $normal_opts, $species_name, $registry;    

    my $cmd_root =  "bsub $rc perl $script $opts --output";

    my $make_cmd = sub {
        my ($name, $extra_opts) = @_;
        
        $name       ||= $species_name;
        $extra_opts ||= '';
        
        if ($name =~ /^_/) {
            $name = $species_name.$name;
        }
        
        my $gvf = "$dir/$name.gvf";
        my $out = "$dir/$name.out";
        my $err = "$dir/$name.err";

        my $cmd = "bsub -o $out -e $err -J $name $rc perl $script $opts $extra_opts --output $gvf";

        return $cmd;
    };

    # first the normal complete variation feature dump

    push @cmds, $make_cmd->();

    # then the complete dump with consequences
    
    push @cmds, $make_cmd->('_incl_consequences', '--include_consequences');

    # then the failed variations

    push @cmds, $make_cmd->('_failed', '--just_failed');

    # then any more specific dumps

    if ($config->{somatic}) {
        push @cmds, $make_cmd->('_somatic', '--somatic');
        push @cmds, $make_cmd->('_somatic_incl_consequences', '--somatic --include_consequences');
    }
    
    if ($config->{structural_variations}) {
        push @cmds, $make_cmd->('_structural_variations', '--just_structural_variations');
    }

    for my $ind (keys %{ $config->{individuals} || {} }) {
        push @cmds, $make_cmd->($ind, "--individual '$ind'");
        
        #my $set = $config->{individuals}->{$ind};
        #push @cmds, $make_cmd->($ind, "--individual '$ind' --set '$set'");
    }
    
    for my $pop (keys %{ $config->{populations} || {} }) {
        my $name = $pop;
        $name =~ s/:/-/g;
        push @cmds, $make_cmd->($name, "--population '$pop'");

        #my $set = $config->{populations}->{$pop};
        #push @cmds, $make_cmd->($name, "--population '$pop' --set '$set'");
    }
}

for my $cmd (@cmds) {
    warn "Submitting:\n$cmd\n\n";
    unless ($test) {
        system($cmd) == 0 or die "Failed to submit last command";
    }
}

__END__

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

