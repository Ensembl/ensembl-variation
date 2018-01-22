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

# see end of file for documentation

use strict;
use warnings;

use FileHandle;
use Getopt::Long;
use Pod::Usage;
use Bio::DB::Fasta;
use Bio::EnsEMBL::Registry;
use Compress::Zlib;

my $species;
my ($gvf_file, $validate_gvf, $vcf_file);
my ($ref_fasta_file, $fasta_file, $fasta_db, $use_fasta_file);
my ($registry_file, $host, $user, $port, $slice_adaptor);
my ($structural_variation);
my ($help, $debug, $compress_output);

$| = 1;
GetOptions(
    "species=s"         => \$species,
    "registry_file|r=s" => \$registry_file,
    "host=s"            => \$host,
    "user=s"            => \$user,
    "port=i"            => \$port,
    "ref_fasta_file=s"  => \$ref_fasta_file,
    "fasta_file=s"      => \$fasta_file,
    "gvf_file=s"        => \$gvf_file,
    "validate_gvf"      => \$validate_gvf,
    "vcf_file=s"        => \$vcf_file,
    "structural_variation" => \$structural_variation,
    "compress"          => \$compress_output,
    "help|h"            => \$help,
    "debug|d"           => \$debug,
) or die pod2usage(1);

pod2usage(1) if $help;

die "Species argument is needed" unless ($species);
die "Input file argument needed" unless ($gvf_file);
die "Input file is not GVF?" unless ($gvf_file =~ m/\.gvf/ || $gvf_file =~ m/\.gvf.gz/);
warn "'Reference fasta file' argument is recommended for the VCF file" unless ($ref_fasta_file);
die "Cannot read from fasta file at the moment" if ($fasta_file);
die "Either a connection to the Ensembl core database or a fasta file is needed" unless ($registry_file || ($host && $user));

my $registry = 'Bio::EnsEMBL::Registry';
if ($host || $registry_file) {
    if ($host) {
        $user ||= 'anonymous';
        if ($port) {
            $registry->load_registry_from_db(-host => $host, -user => $user, -port => $port);
        } else {
            $registry->load_registry_from_db(-host => $host, -user => $user);
        }
    } else {
        $registry->load_all($registry_file);
    }
    my $cdba = $registry->get_DBAdaptor($species, 'core');
    $slice_adaptor = $cdba->get_SliceAdaptor;
} else {
    $fasta_db = Bio::DB::Fasta->new($fasta_file);
    $use_fasta_file = 1;
}

my ($gvf_fh, $vcf_body_fh, $vcf_header_fh);

if ($vcf_file) {
    $vcf_file =~ s/\.vcf//;
} else {
    $vcf_file = $gvf_file;
    if ($gvf_file =~ m/\.gvf/) {
        $vcf_file =~ s/\.gvf//;
    } else {
        $vcf_file =~ s/\.gvf.gz//;
    }
} 

$gvf_fh = gzopen($gvf_file, "rb") or die "Error reading $gvf_file: $gzerrno\n";
$vcf_body_fh = FileHandle->new($vcf_file . '_body.vcf', 'w');
$vcf_header_fh = FileHandle->new($vcf_file . '_header.vcf', 'w');

if ($debug) {
    print 'DEBUG GVF file ', $gvf_file, "\n";
    print 'DEBUG VCF file ', $vcf_file, "\n";
}

my ($validation_states, $evidence, $minor_allele_frequency, $genotype, @sample_ids);
my $map_validation_states = {
    'cluster'     => 'VS_C',
    'freq'        => 'VS_Freq',
    'submitter'   => 'VS_S',
    'doublehit'   => 'VS_DH',
    'hapmap'      => 'VS_HM',
    '1000Genome'  => 'VS_1000G',
    'failed'      => 'VS_Failed',
    'precious'    => 'VS_P',
    '-'           => '',
};
my $map_evidence = {
    'Multiple_observation' => 'E_MO',
    'Frequency'            => 'E_Freq',
    'HapMap'               => 'E_HM',
    '1000Genomes'          => 'E_1000G',
    'Cited'                => 'E_C',
};
my %evidence_levels;
# GVF
my ($seq_id, $source, $type, $start, $end, $score, $strand, $phase, $attrib);
my %gvf_header = ();
my %dbxrefs = ();
my $attributes = {};
my $genotypes;
# VCF
my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $db, $look_up);

parse_gvf2vcf();

sub parse_gvf2vcf {
    parse_gvf_body();  
    if ($gvf_header{End_range} || $gvf_header{Start_range} || $gvf_header{study_accession} || $gvf_header{Parent}) {
        print "Parsing structural variation data is not supported at the moment.\n";
    } 
    print_vcf_header();
    $vcf_header_fh->close();
    $vcf_body_fh->close();
    combine_vcf_header_and_body();
    if ($compress_output) {
        system("gzip -f $vcf_file.vcf") == 0 or die "Failed to gzip $vcf_file.vcf";
    }
    # remove tmp files
    system("rm $vcf_file\_body.vcf") == 0 or die "Failed to rm $vcf_file\_body.vcf";
    system("rm $vcf_file\_header.vcf") == 0 or die "Failed to rm $vcf_file\_header.vcf";
}

sub parse_gvf_header {
    my $count = 0;
    while ($gvf_fh->gzreadline($_) > 0) {
        chomp;
        # parse gvf header
        if ($_ =~ /^##/) {
            unless ($_ =~ /^##sequence-region/) {
                my @values = split(' ', $_);
                $gvf_header{$values[0]} = $values[1];
            }
            next;
        } else {
            # which attributes are set?
            # user defined attribute values should appear in gvf header
            # check attribute columns just in case -> make it an option?
            ($seq_id, $source, $type, $start, $end, $score, $strand, $phase, $attrib) = split/\t/;
            $attributes = get_attributes(\$attrib); 
            map { $gvf_header{$_} = 1 } keys %{$attributes};
            if ($attributes->{Dbxref}) {
                my @dbxref = split(':', $attributes->{Dbxref});
                $dbxrefs{$dbxref[0]} = 1;
            }
        }
    }
    if ($debug) {
        print "DEBUG Parse GVF header\n";
        print join("\n", map { 'DEBUG ' . $_ . ' => ' . $gvf_header{$_} } keys %gvf_header), "\n\n";
    }
    die "Error reading $gvf_file: $gzerrno\n" if $gzerrno != Z_STREAM_END;
    $gvf_fh->gzclose();
    $gvf_fh = gzopen($gvf_file, "rb") or die "Error reading $gvf_file: $gzerrno\n";
}

sub print_vcf_header {
    my $vcf_version = '4.1';

    my ( $sec, $min, $hr, $mday, $mon, $year ) = localtime;
    $year += 1900;    # correct the year
    $mon++;           # correct the month
    my $vcf_file_date = sprintf "%4d%02d%02d", $year, $mon, $mday;

    # use ensembl version used for gvf dump as source
    my $gvf_source_field = $gvf_header{'##data-source'};
    $gvf_source_field =~ s/Source=//;
    my $vcf_source_field = $gvf_source_field . ';This file was parsed from a GVF file.';

    # Meta-information
    print $vcf_header_fh "##fileformat=VCFv$vcf_version", "\n";
    print $vcf_header_fh "##fileDate=$vcf_file_date", "\n";
    print $vcf_header_fh "##source=$vcf_source_field", "\n";
    print $vcf_header_fh "##reference=$ref_fasta_file", "\n"
        if $ref_fasta_file;
    
    #INFO=<ID=ID,Number=number,Type=type,Description=”description”>
    #Types for INFO fields are: Integer, Float, Flag, Character, and String.
    #
    print $vcf_header_fh '##INFO=<ID=TSA,Number=1,Type=String,Description="A SO term describing the type of sequence_alteration (child term of SO sequence_alteration).">', "\n";

    if ($gvf_header{validation_states}) {
        # Where possible dbSNP provides the types of evidence which were used to confirm a variation. We import that information and present it as validation status of a variation.
        print $vcf_header_fh '##INFO=<ID=VS_C,Number=0,Type=Flag,Description="Cluster. Validated by multiple, independent submissions to the refSNP cluster. Source=dbSNP">', "\n";
        print $vcf_header_fh '##INFO=<ID=VS_Freq,Number=0,Type=Flag,Description="Freq. Validated by frequency or genotype data: minor alleles observed in at least two chromosomes. Source=dbSNP">', "\n";
        print $vcf_header_fh '##INFO=<ID=VS_S,Number=0,Type=Flag,Description="Submitter. Validated by submitter confirmation. Source=dbSNP">', "\n";
        print $vcf_header_fh '##INFO=<ID=VS_DH,Number=0,Type=Flag,Description="Doublehit. All alleles have been observed in at least two chromosomes apiece. Source=dbSNP">', "\n";
        print $vcf_header_fh '##INFO=<ID=VS_HM,Number=0,Type=Flag,Description="Hapmap. Genotyped by HapMap project. Source=dbSNP">', "\n";
        print $vcf_header_fh '##INFO=<ID=VS_1000G,Number=0,Type=Flag,Description="1000Genomes. SNP has been sequenced in 1000 Genomes project. Source=dbSNP">', "\n";
        print $vcf_header_fh '##INFO=<ID=VS_P,Number=0,Type=Flag,Description="Precious. Variants for which dbSNP holds citations from PubMed. Source=ensembl">', "\n";
    }
    if ($gvf_header{evidence}) {
        print $vcf_header_fh '##INFO=<ID=E_MO,Number=0,Type=Flag,Description="">', "\n";
        print $vcf_header_fh '##INFO=<ID=E_Freq,Number=0,Type=Flag,Description="">', "\n";
        print $vcf_header_fh '##INFO=<ID=E_HM,Number=0,Type=Flag,Description="">', "\n";
        print $vcf_header_fh '##INFO=<ID=E_1000G,Number=0,Type=Flag,Description="">', "\n";
        print $vcf_header_fh '##INFO=<ID=E_C,Number=0,Type=Flag,Description="">', "\n";
    }
    if ($gvf_header{Variant_effect}) {
        print $vcf_header_fh '##INFO=<ID=SV,Number=1,Type=String,Description="Sequence variant">', "\n";
        print $vcf_header_fh '##INFO=<ID=IDX,Number=1,Type=Integer,Description="0-based index value that identifies which variant sequence the effect is being described for.">', "\n";
        print $vcf_header_fh '##INFO=<ID=FT,Number=1,Type=String,Description="Feature type that is being affected. This term must be the SO term sequence_feature or one of its children.">', "\n";
        print $vcf_header_fh '##INFO=<ID=FID,Number=.,Type=ListOfString,Description="Feature IDs correspond to ID attributes in a GFF3 file that describe the sequence features (for example genes or mRNAs).">', "\n";
        print $vcf_header_fh '##INFO=<ID=VE,Number=.,Type=ListOfString,Description="Effect that a sequence alteration has on a sequence feature that overlaps it.",Format=SV|IDX|FT|FID">', "\n";
    }

    if ($gvf_header{ensembl_failure_reason}) {
        print $vcf_header_fh '##INFO=<ID=EFR,Number=1,Type=String,Description="Ensembl failure reason">', "\n";
    }

    if ($gvf_header{Variant_freq}) {
        print $vcf_header_fh '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">', "\n";
    } else {
        if ($gvf_header{global_minor_allele_frequency}) {
            print $vcf_header_fh '##INFO=<ID=MA,Number=1,Type=String,Description="Minor Alelele">', "\n";
            print $vcf_header_fh '##INFO=<ID=MAF,Number=1,Type=Float,Description="Minor Alelele Frequency">', "\n";
            print $vcf_header_fh '##INFO=<ID=MAC,Number=1,Type=Integer,Description="Minor Alelele Count">', "\n";
        }
    }
    #if (Structural_variation) {
        ##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=“Imprecise structural variation”>
        ##INFO=<ID=NOVEL,Number=0,Type=Flag,Description=“Indicates a novel structural variation”>
        ##INFO=<ID=END,Number=1,Type=Integer,Description=“End position of the variant described in this record”>
        ##INFO=<ID=SVTYPE,Number=1,Type=String,Description=“Type of structural variant”> (DEL, INS, DUP, INV, CNV, BND)

        ##INFO=<ID=DGVID,Number=1,Type=String,Description=“ID of this element in Database of Genomic Variation”>
        ##INFO=<ID=DBVARID,Number=1,Type=String,Description=“ID of this element in DBVAR”>
        ##INFO=<ID=DBRIPID,Number=1,Type=String,Description=“ID of this element in DBRIP”>
    #}

    if ($gvf_header{Genotype}) {
        print $vcf_header_fh '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', "\n";
    }

    foreach my $db (keys %dbxrefs) {
        print $vcf_header_fh "##INFO=<ID=$db,Number=0,Type=Flag,Description=\"$db membership\">\n"; 
    }

    my @header_line = ('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO');
    if ($gvf_header{Genotype}) {
        push @header_line, 'FORMAT';
        #push @header_line, @sample_ids;
    };

    print $vcf_header_fh join("\t", @header_line), "\n";
}

sub parse_gvf_body {
    while ($gvf_fh->gzreadline($_) > 0) {
        chomp;
        if ($_ =~ /^##/) {
            unless ($_ =~ /^##sequence-region/) {
                my @values = split(' ', $_);
                $gvf_header{$values[0]} = $values[1];
            }
        } else {
            parse_row($_);
        }
    }
    die "Error reading $gvf_file: $gzerrno\n" if $gzerrno != Z_STREAM_END;
    $gvf_fh->gzclose();
}

sub parse_row {
    # parse gvf header
    my $line = shift;
   
    # seqid source type start end score strand phase attributes
    ($seq_id, $source, $type, $start, $end, $score, $strand, $phase, $attrib) = split(/\t/, $line);
    $attributes = get_attributes(\$attrib); 
    map { $gvf_header{$_} = 1 } keys %{$attributes};
    if ($attributes->{Dbxref}) {
        my @dbxref = split(':', $attributes->{Dbxref});
        $dbxrefs{$dbxref[0]} = 1;
    }

    ##CHROM
    $chrom = $seq_id;

    # POS, REF, ALT
    $ref = $attributes->{Reference_seq};
    $alt = $attributes->{Variant_seq};
    ($pos, $ref, $alt, $look_up) = @{get_position_and_alleles([$ref, $alt, $start, $end, $chrom])};

    # ID
    my @dbxref = split(':', $attributes->{Dbxref});
    ($id, $db) = ($dbxref[1], $dbxref[0]);
    
    # QUAL assign missing value because quality score is unknown
    $qual = '.';

    # FILTER assign missing value because filters have not been applied/not stored in ensembl variation db 
    $filter = '.'; 
    # INFO ------------------------------------------------------------------------------------------------------------------------
    my @info_data = ();
    push @info_data, "TSA=$type";

    if ($attributes->{ensembl_failure_reason}) {
        push @info_data, "EFA=" . $attributes->{ensembl_failure_reason}; 
    }
    if ($attributes->{Variant_freq}) {
        push @info_data, "AF=" . $attributes->{Variant_freq};    
    } else {
        if ($attributes->{global_minor_allele_frequency} && !$attributes->{Genotype}) {
            # global_minor_allele_frequency=@ 0.435 950
            # global_minor_allele_frequency=0 0.4666 1019
            my $ma;
            my ($index, $maf, $mac) = split(' ', $attributes->{global_minor_allele_frequency});
            if ($maf && $mac) {
                if ($index eq '@') {
                    $ma = $ref;
                } else {
                    my @alleles = split(',', $attributes->{Variant_seq});
                    $ma = $alleles[$index];
                }
                push @info_data, ("MA=$ma", "MAF=$maf", "MAC=$mac");
            }
        }
    }
    if ($attributes->{Variant_effect}) {
        # Variant_effect=sequence_variant index feature_type feature_ID feature_ID
        # consider index -> increase by one to take ref base into account?
        my $effects = join(',', map { join('|', split(' ', $_))} split(',', $attributes->{Variant_effect}) );
        push @info_data, "VE=$effects";
    }
    if ($attributes->{validation_states}) {
        unless ($attributes->{validation_states} eq '-') {
            push @info_data, map { $map_validation_states->{$_} } split(',', $attributes->{validation_states}); 
        }
    } 
    if ($attributes->{evidence}) {
        push @info_data, map { $map_evidence->{$_} } split(',', $attributes->{evidence}); 
    }
    if ($attributes->{clinical_significance}) {
        unless ($attributes->{clinical_significance} eq '-') {
            push @info_data, "CS=$attributes->{clinical_significance}";
        }
    }

    # DB membership
    push @info_data, $db;
    
    # ancestral allele (AA), at the moment not stored in GVF files
    # my $ancestral_alleles = ancestral_allele();
    $info = join(';', @info_data) || '.';

    # END INFO -------------------------------------------------------------------------------------------------------------------------

    # allele count in genotypes (AC)
    # allele frequency for each ALT allele in the same order as listed 

    # Genotypes. Zygosity
    $genotypes = '';
    if ($attributes->{Genotype}) {
        my @gvf_alleles = split(',', $attributes->{Variant_seq});
        my @gvf_genotype_idxs = split(':', $attributes->{Genotype});
        my @vcf_alleles = ();
        push @vcf_alleles, $ref;
        push @vcf_alleles, split(',', $alt);
        my @vcf_gtypes_idxs = ();

        my %look_up_pos = ();
        my $i = 0;
        for my $a (@vcf_alleles) {
            if ($look_up->{$a}) {
                $look_up_pos{$look_up->{$a}} = $i;
            } else {
                $look_up_pos{$a} = $i;
            }
            $i++;
        }
        for my $gt_idx (@gvf_genotype_idxs) {
            my $allele = $gvf_alleles[$gt_idx];
            my $new_idx = $look_up_pos{$allele};
            #my ($new_idx) = grep { $vcf_alleles[$_] eq $allele } 0 .. $#vcf_alleles;
            push @vcf_gtypes_idxs, $new_idx;
        }
        $genotypes = join('/', @vcf_gtypes_idxs);
    }
    $format = '';
    my $sample = '';
    # FORMAT
    if ($genotypes) {
        $format = 'GT';
        $sample = $genotypes;
        print $vcf_body_fh join("\t", $chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample), "\n";
    } else {
        print $vcf_body_fh join("\t", $chrom, $pos, $id, $ref, $alt, $qual, $filter, $info), "\n";
    }
}

sub get_attributes {
    # my %attributes = map { split('=', $_) } split(';', $attrib);
    # need a work around, currently for variant_effect there are more than one key with the same name: 'Variant_effect'
    my $attrib = shift;
    my @pairs = split /;/, $$attrib;
    my @variant_effects = ();
    $attributes = {}; 
    for my $pair (@pairs) {
        my ($key, $value) = split /=/, $pair;
        if ($key eq 'Variant_effect') {
            push @variant_effects, $value;
        } else {
            $attributes->{$key} = $value; 
        } 
    }
    if (@variant_effects) {
        $attributes->{Variant_effect} = join(',', @variant_effects); 
    }
    return $attributes;
}

sub get_position_and_alleles {
    my $input = shift;
    my ($ref, $alt, $start, $end, $seq_region_name) = @{$input};
    my @alts = split(',', $alt);

    my @alleles = ();
    my %allele_lengths = ();

    push @alleles, $ref;
    # if gvf contains genotypes for heterozygous site variant_seq contains all alleles
    push @alleles, grep {$_ ne $ref} @alts;
   
    foreach my $allele (@alleles) {
        $allele =~ s/\-//g;
        $allele_lengths{length($allele)} = 1;
    }
    # in/del/unbalanced
    my $look_up = ();
    if (scalar keys %allele_lengths > 1) {
        # we need ref base before the variation; default to N
        my $prev_base = get_padding_base($start - 1, $start - 1, $seq_region_name);
        for my $i (0..$#alleles) {
            $alleles[$i] =~ s/\-//g;
            $look_up->{$prev_base.$alleles[$i]} = $alleles[$i];
            $alleles[$i] = $prev_base.$alleles[$i];
        }
        $pos = $start - 1;
        $ref = shift @alleles;
    } else {
        $pos = $start;
        shift @alleles;
    }
    $alt = join(',', @alleles); 
    return [$pos, $ref, $alt, $look_up];
}

sub get_padding_base {
    my ($start, $end, $seq_region_name) = @_;
    my $base;
    if ($use_fasta_file) {
        $base = $fasta_db->seq("$seq_region_name:$start,$end"); 
    } else {
        my $slice = $slice_adaptor->fetch_by_toplevel_location("$seq_region_name:$start-$end");
        $base = $slice->seq();
    }
    return $base || 'N';
}

sub combine_vcf_header_and_body {
    my $vcf_fh        = FileHandle->new($vcf_file . '.vcf', 'w');
    my $vcf_header_fh = FileHandle->new($vcf_file . '_header.vcf', 'r');
    while (<$vcf_header_fh>) {
        chomp;
        print $vcf_fh $_, "\n";
    }
    my $vcf_body_fh   = FileHandle->new($vcf_file . '_body.vcf', 'r');

    while (<$vcf_body_fh>) {
        chomp;
        print $vcf_fh $_, "\n";
    }
    $vcf_fh->close();
    return 1;
}


__END__

=head1 NAME 

gvf2vcf.pl

=head1 DESCRIPTION

Covert GVF files to VCF files. GVF files are dumped from variation database.
Only a part of the GVF specification is parsed to VCF. 
VCF files contain customised data information. Please, read VCF header of converted for more information.

=head1 SYNOPSIS

gvf2vcf.pl --species NAME [options]

=head1 EXAMPLE COMMAND LINES
    
  gvf2vcf.pl --gvf_file Pongo_abelii.gvf.gz --ref_fasta_file ftp://ftp.ensembl.org/pub/release-70/fasta/pongo_abelii/dna/Pongo_abelii.PPYG2.70.dna.toplevel.fa.gz --host mysql-ensembl-mirror.ebi.ac.uk --port 4240 --user anonymous --species orangutan --debug

=head1 OPTIONS

=over 4

=item B<--species NAME>

Fetch data for this species, accepts latin name, or any alias set up in your 
registry file. This option is required.

=item B<--gvf_file FILE>

GVF file while which is to be parsed to VCF.

=item B<--validate_gvf>

Check if GVF file as well formated for parsing. Currently not implemented.

=item B<--registry_file FILE>

Use database connection details from this registry file, will default to the
registry specified in the $ENSEMBL_REGISTRY environment variable if not
specified.

=item B<--host>

Alternatively to a registry file host, user and port can be defined to connect to a core database.

=item B<--user>

Alternatively to a registry file host, user and port can be defined to connect to a core database.

=item B<--port>

Alternatively to a registry file host, user and port can be defined to connect to a core database.

=item B<--ref_fasta_file>

A location to a fasta file which is used as reference for variations in the VCF file

=item B<--vcf_file>

Dump VCF to the specified file. Defaults to same location as the given GVF file but file extension will be changed to .vcf 

=item B<--compress>

gzip the file when the dump has finished

=item B<--debug>

Display progress on intermediate results of parsing 

=item B<--help>

Display this documentation

=head1

For help with this script address questions to http://lists.ensembl.org/mailman/listinfo/dev

=cut
