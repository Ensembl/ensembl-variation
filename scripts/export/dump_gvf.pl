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


# dump_gvf.pl - Create GVF dumps from Ensembl variation databases
#
# see end of file for documentation

use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Utils::EnsEMBL2GFF3;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

my $species;
my $output_file;
my $registry;
my $host;
my $user;
my $somatic;
my $compress;
my $use_iterator;
my $include_failed;
my $just_failed;
my $individual_name;
my $population_name;
my $set_name;
my $chunk_size;
my $include_consequences;
my $include_coding_details;
my $include_global_maf;
my @seq_region_names;
my $include_svs;
my $just_svs;
my $help;

$| = 1;

GetOptions(
    "species|s=s"                   => \$species,
    "output|o=s"                    => \$output_file,
    "registry|r=s"                  => \$registry,
    "host=s"                        => \$host,
    "user=s"                        => \$user,
    "individual|i=s"                => \$individual_name,
    "population|p=s"                => \$population_name,
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
    "seq_region|sr=s"               => \@seq_region_names,
    "include_structural_variations" => \$include_svs,
    "just_structural_variations"    => \$just_svs,
    "help|h"                        => \$help,
) or die pod2usage(1);

pod2usage(1) if $help;

# the only parameter we require is species
die "species argument required, try --help for usage instructions\n" unless $species;

die "Can't fetch for a population and an individual at once" 
    if $population_name && $individual_name;

# chunk size is in kilobases
$chunk_size *= 1000 if defined $chunk_size;

# if the user wants just structural variants that implies that SVs should be included
$include_svs = 1 if $just_svs;

# default to a sensible file name
$output_file ||= "$species.gvf";

my $reg = 'Bio::EnsEMBL::Registry';

if ($host) {
    # if we are supplied with a host, try to load the registry from there
    $user ||= 'anonymous';
    $reg->load_registry_from_db(-host => $host, -user => $user);
}
else {
    # otherwise use the registry file supplied (or default to the 
    # ENSEMBL_REGISTRY environment variable)
    $reg->load_all($registry);
}

my $cdba = $reg->get_DBAdaptor($species, 'core') 
    or die "Failed to get core DBAdaptor";

my $vdba = $reg->get_DBAdaptor($species, 'variation') 
    or die "Failed to get variation DBAdaptor";

my $failed_descs;
my $failed_ids;

if ($include_failed || $just_failed) {
    
    # include failed variations, but flag them with an attribute
    $vdba->include_failed_variations(1);

    # joining to variation and failed variation is too slow, so
    # we cache the entire failed_variation table in memory, just 
    # storing the failed_description_id to save memory. We also 
    # store the failed_description table in memory to map these 
    # ids to the description

    my $dbh = $vdba->dbc->db_handle;

    my $sth = $dbh->prepare(qq{
        SELECT  failed_description_id, description
        FROM    failed_description
    });

    $sth->execute;

    while (my ($id, $desc) = $sth->fetchrow_array) {
        $failed_descs->{$id} = $desc;
    }

    $sth = $dbh->prepare(qq{
        SELECT  variation_id, failed_description_id
        FROM    failed_variation
    });

    $sth->execute;
    
    my $v_id;
    my $desc_id;

    $sth->bind_columns(\$v_id, \$desc_id);

    while ($sth->fetch) {
        $failed_ids->{$v_id} = $desc_id;
    }
}


my $individual;

if ($individual_name) {
    my $ia = $vdba->get_IndividualAdaptor;
    $individual = $ia->fetch_all_by_name($individual_name)->[0] 
        or die "Didn't find individual '$individual_name'";
}

my $population;

my $alleles_by_variation_id;

if ($population_name) {
    my $pa = $vdba->get_PopulationAdaptor;
    $population = $pa->fetch_by_name($population_name)
        or die "Didn't find population '$population_name'";

    # joining to variation and allele etc. is too slow so we cache
    # all alleles and frequencies for the population in memory

    my $dbh = $vdba->dbc->db_handle;
    my $sth = $dbh->prepare(qq{
        SELECT  a.variation_id, ac.allele, a.frequency
        FROM    allele a, allele_code ac
        WHERE   a.allele_code_id = ac.allele_code_id
        AND     population_id = ?
    });

    $sth->execute($population->dbID);

    my $v_id;
    my $allele;
    my $freq;

    $sth->bind_columns(\$v_id, \$allele, \$freq);

    while ($sth->fetch) {
        $alleles_by_variation_id->{$v_id}->{$allele} = defined $freq ? $freq : '';
    }
}

my $set;

if ($set_name) {
    my $vsa = $vdba->get_VariationSetAdaptor;
    $set = $vsa->fetch_by_name($set_name)
        or die "Didn't find set '$set_name'";
}

# open the output file

open GVF, ">$output_file" or die $!;

# get our slices

my $slices;

my $sa = $cdba->get_SliceAdaptor;

if (@seq_region_names) {
    for my $name (@seq_region_names) {
        my $slice = $sa->fetch_by_region('toplevel', $name) 
            or die "Failed to get a slice for seq region name: $name";
        push @$slices, $slice;
    }
}
else {

#    # find all the unique seq_region_ids in the database
#    # and fetch a slice for each of them
#
#    my $dbh = $vdba->dbc->db_handle;
#   
#    my %sr_ids;
#
#    if ($include_svs) {
#        my $sth = $dbh->prepare(qq{
#            SELECT  DISTINCT(seq_region_id)
#            FROM    structural_variation
#        });
#
#        $sth->execute;
#
#        while (my ($sr_id) = $sth->fetchrow_array) {
#            $sr_ids{$sr_id} = 1;
#        }
#    }
#
#    unless ($just_svs) {
#        my $sth = $dbh->prepare(qq{
#            SELECT  DISTINCT(seq_region_id)
#            FROM    variation_feature
#        });
#
#        $sth->execute;
#        
#        while (my ($sr_id) = $sth->fetchrow_array) {
#            $sr_ids{$sr_id} = 1;
#        }
#    }
#
#    for my $sr_id (keys %sr_ids) {
#        my $slice = $sa->fetch_by_seq_region_id($sr_id);
#        die "seq_region_id $sr_id found in variation database but not in core?" unless $slice;
#        push @$slices, $slice;
#    }
    
    #$slices = $sa->fetch_all('toplevel', undef, 1, 0);
    $slices = $sa->fetch_all('toplevel');
    die "Didn't find any toplevel slices for $species?" unless @$slices;
}

# print the GVF header for the first slice, without sequence_region information
# because we add that for all slices later, but include details of the individual
# or population as necessary

print GVF $slices->[0]->gvf_header(
    no_sequence_region  => 1, 
    individual          => $individual,
    population          => $population,
);

# print a sequence-region line in the GVF file for each slice

for my $slice (@$slices) {
    print GVF '##sequence-region ', $slice->seq_region_name, ' ',$slice->start, ' ', $slice->end, "\n"; 
}

# just use a file-wide count for the ID attribute
my $id_count = 0;

my $prev_svs;

while (my $slice = shift @$slices) {

    # only operate on chunk_size sub slices to avoid running out of memory

    my $sub_slice_start = 1;

    $chunk_size ||= $slice->end;

    #warn "Fetching by slice: ".$slice->name." chunk size: ".$chunk_size;
    
    while ($sub_slice_start < $slice->end) {

        my $sub_slice_end = $sub_slice_start + $chunk_size - 1;
        
        $sub_slice_end = $slice->end if $sub_slice_end > $slice->end;

        my $sub_slice = $slice->sub_Slice($sub_slice_start, $sub_slice_end)
            or die "No sub slice for $sub_slice_start - $sub_slice_end?";
        
        $sub_slice_start = $sub_slice_end + 1;

        unless ($just_svs) {

            my $vfa  = $vdba->get_VariationFeatureAdaptor;

            my $vfs;

            if ($use_iterator) {
                $vfs = $vfa->fetch_Iterator_by_Slice_constraint($sub_slice, '');
            }
            elsif ($somatic) {
                $vfs = $vfa->fetch_all_somatic_by_Slice_constraint_with_TranscriptVariations($sub_slice);
            }
            elsif ($set) {
                $vfs = $vfa->fetch_all_by_Slice_VariationSet($sub_slice, $set);
            }
            elsif ($include_consequences) {
                $vfs = $vfa->fetch_all_by_Slice_constraint_with_TranscriptVariations($sub_slice);
            }
            else {
                $vfs = $vfa->fetch_all_by_Slice($sub_slice);
            }

            # we pre-fetch all the genotypes for this slice and individual to 
            # speed things up below

            my $gt_hash;

            if ($individual) {

                my $igta = $vdba->get_IndividualGenotypeAdaptor;

                my $igts = $igta->fetch_all_by_Slice($sub_slice, $individual);

                for my $igt (@$igts) {

                    my $key = join('_', 
                        $igt->slice->seq_region_name,
                        $igt->seq_region_start,
                        $igt->seq_region_end,
                    );

                    my $gts = $gt_hash->{$key} ||= [];

                    my $seen = 0;

                    if (@$gts) {
                        for my $gt (@$gts) {
                            if (($gt->allele1 eq $igt->allele1 && $gt->allele2 eq $igt->allele2) ||
                                ($gt->allele1 eq $igt->allele2 && $gt->allele2 eq $igt->allele1)) {
                                $seen++;
                                last;
                            }
                        }
                    }

                    push @$gts, $igt unless $seen;
                }
            } # end individual

            # iterate over our variation_features

            VF : while (my $vf = ($use_iterator ? $vfs->next : shift @$vfs)) {

                my $attrs;

                my $comment = '';

                if ($individual) { 

                    # we're making a file for an individual, so fetch the genotypes

                    my $key = join('_',
                        $vf->slice->seq_region_name,
                        $vf->seq_region_start,
                        $vf->seq_region_end,
                    );

                    my @gts = @{ $gt_hash->{$key} || [] };

                    GT : for my $gt (@gts) { 

                        my $gt1 = $gt->allele1;
                        my $gt2 = $gt->allele2;

                        unless ($vf->strand == $gt->strand) {
                            reverse_comp(\$gt1);
                            reverse_comp(\$gt2);
                        }

                        my @alleles = split /\//, $vf->allele_string;

                        my $ref_allele = shift @alleles;

                        my $hom = $gt1 eq $gt2;

                        my $gt1_is_ref = $gt1 eq $ref_allele;
                        my $gt2_is_ref = $gt2 eq $ref_allele;

                        my $gt1_in_alt;
                        my $gt2_in_alt;

                        if ($hom && $gt1_is_ref) {
                            # homozygous for the reference, don't include
                            next VF;
                        }

                        my $variant_seq;

                        for my $allele (@alleles) {
                            $gt1_in_alt = 1 if $allele eq $gt1;
                            $gt2_in_alt = 1 if $allele eq $gt2;
                            last if ($gt1_in_alt && $gt2_in_alt);
                        }

                        if ($hom && $gt1_in_alt) {
                            # homozygous for alt
                            $variant_seq = "$gt1";
                        }
                        elsif ($gt1_is_ref && $gt2_in_alt) {
                            # het for ref and alt
                            $variant_seq = "$gt1,$gt2";
                        }
                        elsif ($gt2_is_ref && $gt1_in_alt) {
                            # het for ref and alt
                            $variant_seq = "$gt2,$gt1";
                        }
                        elsif ($gt1_in_alt && $gt2_in_alt) {
                            # het for 2 alts
                            $variant_seq = "$gt1,$gt2";
                        }
                        else {

                            # the genotype sequence isn't in the allele string

                            # first check if the classes match

                            my $gt_incl_indel = ($gt1 eq '-' || $gt eq '-');
                            my $vf_is_indel = ($vf->var_class eq 'deletion' || $vf->var_class eq 'insertion');

                            next GT unless $gt_incl_indel == $vf_is_indel;

                            my $vf_is_sub = $vf->var_class eq 'substitution';
                            my $gt_is_sub = (length $gt1 > 1) || (length $gt2 > 1);

                            next GT unless $vf_is_sub == $gt_is_sub;

                            # sometimes there are multiple genotypes at the same location,
                            # but this only seems to happen when the associated variation
                            # differs

                            next GT if ($gt->variation->name ne $vf->variation_name);

                            # otherwise sometimes the variations match, but the alleles
                            # don't match the genotype e.g. rs55844409 has allele string 
                            # T/C but Watson has a G/G genotype. We suspect these are dbSNP 
                            # bugs, but we still include them

                            $variant_seq = $hom ? $gt1 : "$gt1,$gt2";
                        }

                        $attrs->{Variant_seq} = $variant_seq;

                        $attrs->{Zygosity} = $hom ? 'homozygous' : 'heterozygous';	

						my @variant_seqs = split ',', $variant_seq;
						my $index = 0;
						my %allele_index = map {$_ => $index++} @variant_seqs;
						my $genotype = $allele_index{$gt1} . ":" . $allele_index{$gt2};
						$attrs->{Genotype} = $genotype;	
                        #$comment .= "real genotype: $gt1/$gt2"
                    }

                    unless ($attrs->{Zygosity}) {
                        # this vf isn't in this individual, so don't output it
                        next VF;
                    }
                }

                if ($population) {

                    # we're making a file for a population, so look up the alleles and their frequencies

                    my @vf_alleles = split /\//, $vf->allele_string;

                    my $found = 0;

                    my @alleles;
                    my @freqs;

                    if (my $allele_hash = $alleles_by_variation_id->{$vf->{_variation_id}}) {
                        for my $allele (keys %$allele_hash) {
                            my $freq = $allele_hash->{$allele};
                            if ($allele eq $vf_alleles[0]) {
                                $attrs->{Reference_seq} = $allele;
                            }
                            else {
                                push @alleles, $allele;
                                push @freqs, $freq;
                                $found = 1;
                            }
                        }
                    }

                    if ($found) {
                        $attrs->{Variant_seq}   = join ',', @alleles;

                        if (grep {$_ =~ /\d/} @freqs) {
                            $attrs->{Variant_freq}  = join ',', @freqs;
                        }

                        #warn $vf->to_gvf(extra_attrs => $attrs)."\n";
                    }
                    else {
                        next VF;
                    }
                }

                if ($include_failed || $just_failed) {

                    # we're including failed variations, so look up the description and add to the output

                    if (my $desc_id = $failed_ids->{$vf->{_variation_id}}) {
                        my $desc = $failed_descs->{$desc_id};
                        $attrs->{ensembl_failure_reason} = $desc;
                    }
                    else {
                        next VF if $just_failed;
                    }
                }

                # if we get here, we want this vf included in the dump, so convert
                # it to GVF including any extra attributes defined above

                $attrs->{ID} = ++$id_count;

                my $gvf_line = $vf->to_gvf(
                    extra_attrs             => $attrs, 
                    include_consequences    => $include_consequences,
                    include_coding_details  => $include_coding_details,
                    include_global_maf      => $include_global_maf,
                );

                if ($id_count % 1000 == 1) {
                    #check_mem();
                }

                print GVF $gvf_line, ($comment ? " # $comment" : ''), "\n" if $gvf_line;
            } # end variation feature
        }

        if ($include_svs) {

            my $svfa = $vdba->get_StructuralVariationFeatureAdaptor;
            
            my $svfs = $svfa->fetch_all_by_Slice($sub_slice, 1);
            
            for my $svf (@$svfs) {
               
                # ignore CNV probes

                next if $svf->var_class eq 'CNV_PROBE';

                my $coords = join '-', $svf->seq_region_name, $svf->seq_region_start, $svf->seq_region_end; 
                
                if (my $prev_coords = $prev_svs->{$svf->variation_name}) {
                    warn "repeated SV: ".$svf->variation_name." coords 1: $prev_coords, coords 2: $coords\n" if $prev_coords ne $coords;
                    next;
                }

                # we can now have SVs that map to multiple locations so we can't use the
                # feature's own identifier and we have to use a file-wide count as for
                # normal variations

                my $gvf_line = $svf->to_gvf(extra_attrs => {ID => ++$id_count});
                    
                print GVF "$gvf_line\n" if $gvf_line;
                
                $prev_svs->{$svf->variation_name} = $coords if $gvf_line;
            }
        }
    }
}

close GVF;

if ($compress) {
    system("gzip -f $output_file") == 0 or die "Failed to gzip $output_file";
}

my $config = {};

# finds out memory usage
sub memory {
    my @mem;
    
    open IN, "ps -o rss,vsz $$ |";
    while(<IN>) {
        next if $_ =~ /rss/i;
        chomp;
        @mem = split;
    }
    close IN;
    
    return \@mem;
}


sub check_mem {
    eval q{use Devel::Size qw(total_size)};
    my $mem = memory();
    my $tot;
    $tot += $_ for @$mem;

    if($tot > 1000000) {
        $tot = sprintf("%.2fGB", $tot / (1024 * 1024));
    }

    elsif($tot > 1000) {
        $tot = sprintf("%.2fMB", $tot / 1024);
    }

    my $mem_diff = mem_diff($config);
    print "MEMORY $tot ", (join " ", @$mem),"\nDIFF ", (join " ", @$mem_diff), "\n";
}

sub mem_diff {
    my $config = shift;
    my $mem = memory();
    my @diffs;
    
    if(defined($config->{memory})) {
        for my $i(0..(scalar @{$config->{memory}} - 1)) {
            push @diffs, $mem->[$i] - $config->{memory}->[$i];
        }
    }
    else {
        @diffs = @$mem;
    }
    
    $config->{memory} = $mem;
    
    return \@diffs;
}

__END__

=head1 NAME 

dump_gvf.pl

=head1 DESCRIPTION

Create GVF files from the variation databases

=head1 SYNOPSIS

dump_gvf.pl --species NAME [options]

=head1 EXAMPLE COMMAND LINES

  dump_gvf.pl --species mouse --registry ensembl.registry --include_failed

    dump all variation features, including failed variations, from the mouse 
    database to the default mouse.gvf file

  dump_gvf.pl --species human --registry ensembl.registry --individual Watson \
    --seq_region 1 --seq_region 2 --output watson.gvf

    dump all variations, including genotypes, found in Watson's chromosomes 
    1 & 2 to the watson.gvf file

  dump_gvf.pl --species human --population '1000GENOMES:low_coverage:CEU' \
    --output 1kg_pilot1_CEU.gvf --compress

    dump all variations, including allele frequencies, found in the specified
    population and then compress the file to 1kg_pilot1_CEU.gvf.gz

  dump_gvf.pl --species human --somatic --output human_somatic.gvf

    dump somatic mutations from human

=head1 OPTIONS

=over 4

=item B<--species NAME>

Fetch data for this species, accepts latin name, or any alias set up in your 
registry file. This option is required.

=item B<--registry FILE>

Use database connection details from this registry file, will default to the
registry specified in the $ENSEMBL_REGISTRY environment variable if not
specified.

=item B<--output>

Dump GVF to the specified file. Defaults to <species_name>.gvf

=item B<--somatic>

Fetch somatic mutations rather than germline variations

=item B<--compress>

gzip the file when the dump has finished

=item B<--include_failed>

Include variations flagged as failed in the GVF file

=item B<--just_failed>

Only include variations flagged as failed

=item B<--include_consequences>

Include the consequences of variations on transcripts and regulatory regions etc.

=item B<--include_coding_details>

Include extra information for coding variations, including the amino acid alleles and any SIFT or PolyPhen predictions

=item B<--include_global_maf>

Include the global minor allele frequency information from dbSNP 

=item B<--individual NAME>

Only include variations found in this individual, including genotypes if available.
B<NB> You can't specify an individual and a population at the same time.

=item B<--population NAME>

Only include variations found in this population, including allele frequencies if available
B<NB> You can't specify a population and an individual at the same time.

=item B<--set NAME>

Only include variations found in this variation set

=item B<--chunk_size NUM_KILOBASES>

Fetch features from slices this number of kilobases long to avoid excessive memory usage

=item B<--seq_region NAME>

Only include variations from this toplevel seq_region (normally chromosome name, e.g. 4 or X).
You can specify several seq_regions by repeating this option for each name.

=item B<--include_structural_variations>

Include structural variations in the output

=item B<--just_structural_variations>

Only include structural variations in the output

=item B<--help>

Display this documentation

=head1

For help with this script address questions to http://lists.ensembl.org/mailman/listinfo/dev

