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


=head1 NAME
dump_vcf.pl - dumps variations from variation DB into VCF file

=cut

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Iterator;
use Bio::EnsEMBL::Utils::Sequence qw(expand);
use Bio::EnsEMBL::Utils::Slice qw(split_Slices);
use Bio::EnsEMBL::Utils::Sequence qw(expand);
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;


use Getopt::Long;
use FileHandle;
use List::Util qw(first);
$| = 1;

my $config = {};
GetOptions(
    $config,
    'help|h',
    'output_file|o=s',
    'species=s',
    'registry|r=s',
    'host=s',
    'user=s',
    'port=i',
    'db_version=i',
    'population=s',
    'individuals=s',
    'allele_freq_in_population=s',
    'seq_regions=s',
    'consequences',
    'protein_coding_details',
    'ancestral_allele',
    'global_maf',
    'evidence_values',
    'ref_fasta_file=s',
    'data_source=s',
) or die "ERROR: Failed to parse command line arguments";

if (defined($config->{help})) {
    &usage;
    exit(0);
}

die "species argument required, try --help for usage instructions\n" unless $config->{species};

die "Can't fetch for a population and an individual at once" 
    if $config->{population} && $config->{individuals};

# default to a sensible file name
my $species = $config->{species};
$config->{output_file} ||= "$species.vcf";

my $reg = 'Bio::EnsEMBL::Registry';

if ( defined($config->{host}) && defined($config->{user})) {
    $config->{port} ||= 3306;
    $reg->load_registry_from_db(-host => $config->{host}, -user => $config->{user}, -port => $config->{port}, -db_version => $config->{db_version});
} else {
    if (-e $config->{registry}) {
        $reg->load_all($config->{registry});
	} else {
        die "ERROR: could not read from registry file ".$config->{registry}."\n";
	}
}
	
# connect to DB
my $vdba = $reg->get_DBAdaptor($config->{species},'variation') || usage( "Cannot find variation db for ".$config->{species}." in ".$config->{registry_file} );
my $cdba = $reg->get_DBAdaptor($config->{species},'core') || usage( "Cannot find core db for ".$config->{species}." in ".$config->{registry_file} );
my $sliceAdaptor      = $cdba->get_SliceAdaptor;
my $vfAdaptor         = $vdba->get_VariationFeatureAdaptor;
my $populationAdaptor = $vdba->get_PopulationAdaptor; 
my $alleleAdaptor     = $vdba->get_AlleleAdaptor;
my $igtAdaptor        = $vdba->get_IndividualGenotypeAdaptor;

my $individuals;
my $population;
if ($config->{allele_freq_in_population}) {
    print $config->{allele_freq_in_population}, "\n";
    $population = $populationAdaptor->fetch_by_name($config->{allele_freq_in_population});
}
if ($config->{population}) {
    $population = $populationAdaptor->fetch_by_name($config->{population});
    $individuals = $population->get_all_Individuals();
    $config->{sample_data} = 1;
    $config->{format} = 'GT';
}
if ($config->{individuals}) {
    my $individualAdaptor = $vdba->get_IndividualAdaptor;
    if (uc $config->{individuals} eq 'ALL') {
      $individuals = $individualAdaptor->fetch_all();
    } else {
      my @names = split(',', $config->{individuals});
      foreach my $name (@names) {
          my $individual_objects = $individualAdaptor->fetch_all_by_name($name);
          die "More than one individual for name $name." if (scalar @$individual_objects > 1);
          push @$individuals, $individual_objects->[0];
      }
    }
    $config->{sample_data} = 1;
    $config->{format} = 'GT';
}

if ($config->{sample_data}) {
    foreach my $individual (@$individuals) {
        $config->{sample}->{$individual->name} = 1;
    } 
}
my $map_evidence = {
    'Multiple_observations' => 'E_MO',
    'Frequency'             => 'E_Freq',
    'HapMap'                => 'E_HM',
    '1000Genomes'           => 'E_1000G',
    'Cited'                 => 'E_C',
};

my $slices;
if ($config->{seq_regions}) {
    my @regions = split(',', $config->{seq_regions});
    foreach my $region (@regions) {
        my $slice = $sliceAdaptor->fetch_by_region('toplevel', $region);
        push @$slices, $slice;
    }
} else {
    # include only non-duplicate regions
    $slices = $sliceAdaptor->fetch_all('toplevel', undef, 0, 1);
}
my $max_length = 1e4;
my $overlap = 0;
my ($slice_pieces, $vf_it, $vfs);

my $file_handle = FileHandle->new("> ".$config->{output_file});

print_header($config);

foreach my $slice (@$slices) {
    $slice_pieces = split_Slices([$slice], $max_length, $overlap);
    foreach my $slice_piece (@$slice_pieces) {
        $vf_it = $vfAdaptor->fetch_Iterator_by_Slice($slice_piece);
        while (my $vf = $vf_it->next) {
            next if ($vf->seq_region_start <= $slice_piece->start); # avoid duplicated lines caused by vf overlapping two slice pieces
            my $vcf_line = {};
            $vcf_line->{CHROM} = $vf->seq_region_name;
            $vcf_line->{ID} = $vf->variation_name;
            my @alleles = split /\//, $vf->allele_string;
            map {expand(\$_)} @alleles;
            $vcf_line->{allele_string} = join('/', @alleles);
            next unless alleles($config, $vf, $vcf_line);
            my $source = $vf->source;
            $source .= '_' . $vf->source_version if defined $vf->source_version;
            $vcf_line->{INFO}->{DB} = $source;
            if ($config->{sample_data}) {
                next unless (genotypes($vcf_line, $vf));
            }
            if ($config->{consequences} || $config->{protein_coding_details}) {
                consequences($config, $vcf_line, $vf);
            }
            if ($config->{global_maf}) {
                global_maf($vcf_line, $vf);
            }
            if ($config->{allele_freq_in_population}) {
                allele_freq_in_population($vcf_line, $vf);
            }
            if ($config->{ancestral_allele}) {
                ancestral_allele($vcf_line, $vf);
            }
            if ($config->{evidence_values}) {
                evidence_values($vcf_line, $vf);
            }
            print_vcf_line($config, $vcf_line);
        }
    }
}
$file_handle->close();

sub alleles {
    my $config = shift;
    my $vf = shift;
    my $vcf_line = shift;
    my %allele_lengths;
    my @alleles = split /\//, $vcf_line->{allele_string};
    # quality control alleles
    my $qc = 1;
    foreach my $allele (@alleles) {
        unless ($allele=~/(^[ACTGNacgtn]+$)|(^-$)|(^<.*>$)/) {
            warn "$allele in $vcf_line->{ID} did not pass allele quality control checks";
            $qc = 0;
        }
    }
    if ($qc) {
        # look for imbalance in the allele string

        foreach my $allele(@alleles) {
            $allele =~ s/\-//g;
            $allele_lengths{length($allele)} = 1;
        }

        my $start = $vf->seq_region_start;
        my $end = $vf->seq_region_end;
        if ($start > $end) {
            $start = $end;
        }
        # in/del/unbalanced
        if(scalar keys %allele_lengths > 1) {
            
            # we need the ref base before the variation
            # default to N in case we can't get it
            my $prev_base = 'N';
            my $loc = $start - 1;
            my $seq_region_name = $vcf_line->{CHROM};
            my $slice = $sliceAdaptor->fetch_by_toplevel_location("$seq_region_name:$loc-$loc");
            $prev_base = $slice->seq if defined($slice);
            
            for my $i(0..$#alleles) {
                $alleles[$i] =~ s/\-//g;
                $alleles[$i] = $prev_base.$alleles[$i];
            }
            $vcf_line->{POS} = $start - 1;
            $vcf_line->{REF} = shift @alleles;
            $vcf_line->{ALT} = join(',', @alleles);
        }
        # balanced sub
        else {
            $vcf_line->{POS} = $start;
            $vcf_line->{REF} = shift @alleles;
            $vcf_line->{ALT} = join(',', @alleles);
        }
    }
    return $qc;
}

sub genotypes {
    my $vcf_line = shift;
    my $vf = shift;
    my $has_gt = 0;
    foreach my $individual (@$individuals) {
        my $igts = $igtAdaptor->fetch_all_by_Variation($vf->variation, $individual);
        if (scalar @$igts > 0) {
            my $igt = $igts->[0];
            my @alleles = split /\//, $vcf_line->{allele_string};
            my $allele1_idx = (first { $alleles[$_] eq $igt->genotype->[0] } 0..$#alleles) || '.';
            my $allele2_idx = (first { $alleles[$_] eq $igt->genotype->[1] } 0..$#alleles) || '.';
            my $genotype = "$allele1_idx/$allele2_idx";
            $vcf_line->{sample}->{$individual->name} = $genotype;
            $has_gt = 1;
        } else {
            $vcf_line->{sample}->{$individual->name} = '.';
        }
    }
    return $has_gt;
}

sub consequences {
    my $config = shift;
    my $vcf_line = shift;
    my $vf = shift;
    my $tvs = $vf->get_all_TranscriptVariations;
    my @alleles = split /\//, $vcf_line->{allele_string};
    #shift @alleles; # shift off reference allele?
    foreach my $tv (@$tvs) {
        #push @return, map {tva_to_line($config, $_)} @{$tv->get_all_alternate_TranscriptVariationAlleles};
        foreach my $tva (@{$tv->get_all_alternate_TranscriptVariationAlleles}) {
            my $allele_idx = first { $alleles[$_] eq $tva->variation_feature_seq } 0..$#alleles;
            if ($config->{consequences}) {
                for my $oc (@{$tva->get_all_OverlapConsequences}) {
                    push @{ $vcf_line->{INFO}->{Variant_effect} ||= [] },
                                    join('|',
                                        $oc->SO_term,
                                        $allele_idx,
                                        $oc->feature_SO_term,
                                        $tv->transcript_stable_id,);
                }    
            }
            if ($config->{protein_coding_details}) {
                if ($tva->pep_allele_string) {
                    push @{ $vcf_line->{INFO}->{variant_peptide} ||= [] },
                        join('|',
                            $allele_idx,
                            $tva->peptide,
                            $tv->transcript_stable_id,);

                    for my $tool (qw(sift polyphen)) {
                        my $pred_meth = $tool.'_prediction';
                        my $score_meth = $tool.'_score';
                        if (my $pred = $tva->$pred_meth) {
                            $pred =~ s/\s/_/g;
                            push @{ $vcf_line->{INFO}->{$tool."_prediction"} ||= [] }, 
                                join('|', 
                                    $allele_idx, 
                                    $pred, 
                                    $tva->$score_meth,
                                    $tv->transcript_stable_i,);
                        }
                    }
                }
            }
        }
    }
} 

sub global_maf {
    my $vcf_line = shift;
    my $vf = shift;
    my $variation = $vf->variation;
    if (defined($variation->minor_allele_frequency)) {
        my @alleles = split /\//, $vcf_line->{allele_string};
        my $allele_idx = first { $alleles[$_] eq $variation->minor_allele } 0..$#alleles;
        if (defined($allele_idx)) {
            $vcf_line->{INFO}->{GMAF} = 
                join('|',
                    $allele_idx,
                    $variation->minor_allele_frequency,
                    $variation->minor_allele_count,);
        }
    }
}

sub ancestral_allele {
    my $vcf_line = shift;
    my $vf = shift;
    my $variation = $vf->variation;
    if (defined($variation->ancestral_allele)) {
        $vcf_line->{INFO}->{AA} = $variation->ancestral_allele;
    }
}

sub evidence_values {
    my $vcf_line = shift;
    my $vf = shift;
    my $variation = $vf->variation;
    my $values = $variation->get_all_evidence_values();
    if (scalar @$values) {
        push @{$vcf_line->{FLAG}}, map {$map_evidence->{$_}} @$values;
    }
}

sub allele_freq_in_population {
    my $vcf_line = shift;
    my $vf = shift;
    my @alleles_in_population = @{$alleleAdaptor->fetch_all_by_Variation($vf->variation, $population)};
    return unless (@alleles_in_population);
    my $allele_freqs;
    foreach my $allele (@alleles_in_population) {
        $allele_freqs->{$allele->allele()} = $allele->frequency();
    }
    my @alleles = split /\//, $vcf_line->{allele_string};
    shift @alleles; # shift off reference allele?
    my @af = ();
    foreach my $allele (@alleles) {
        my $freq = $allele_freqs->{$allele};
        $freq ||= 0.0;
        push @af, $freq;
    }
    $vcf_line->{INFO}->{AF} = join(',', @af);
}

sub print_header {
    my $config = shift;
    my $version = '4.1';
    my ($sec, $min, $hr, $mday, $mon, $year) = localtime;
    $year += 1900; # correct the year
    $mon++;        # correct the month
    my $file_date = sprintf "%4d%02d%02d", $year, $mon, $mday;

    print $file_handle "##fileformat=VCFv$version\n";
    print $file_handle "##fileData=$file_date\n";

    if ($config->{ref_fasta_file}) {
        print $file_handle "##reference=$config->{ref_fasta_file}\n";
    }
    if ($config->{data_source}) {
        print $file_handle "##source=$config->{data_source}\n";
    } else {
        my $mca = $cdba->get_MetaContainerAdaptor;
        my $schema_version = $mca->get_schema_version;
        my $species_name = $mca->get_scientific_name;
        $species_name =~ s/ /_/g;
        my $url = 'http://e'.$schema_version.'.ensembl.org/'.$species_name;
        print $file_handle "##source=ensembl,version=$schema_version,url=$url\n";
    }
    print $file_handle "##INFO=<ID=DB,Number=1,Type=String,Description=\"Source for variant and if available source version\">\n";
    if ($config->{allele_freq_in_population}) {
        print $file_handle "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n";
    }
    if ($config->{ancestral_allele}) {
        print $file_handle "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n";
    }
    if ($config->{global_maf}) {
        print $file_handle "##INFO=<ID=MAIDX,Number=1,Type=String,Description=\"Minor Allele (starting with 0 for REF allele)\">\n";
        print $file_handle "##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor Allele Frequency\">\n";
        print $file_handle "##INFO=<ID=MAC,Number=1,Type=Integer,Description=\"Minor Allele Count\">\n";
        print $file_handle "##INFO=<ID=GMAF,Number=1,Type=ListOfString,Description=\"Global Minor Allele Frequency.\",Format=MAIDX|MAF|MAC>\n";
    }
    if ($config->{consequences}) {
        print $file_handle "##INFO=<ID=SV,Number=1,Type=String,Description=\"Sequence variant\">\n";
        print $file_handle "##INFO=<ID=IDX,Number=1,Type=Integer,Description=\"0-based (starting with first reported ALT allele) index value that identifies which variant sequence the effect is being described for.\">\n";
        print $file_handle "##INFO=<ID=FT,Number=1,Type=String,Description=\"Feature type that is being affected. This term must be the SO term sequence_feature or one of its children.\">\n";
        print $file_handle "##INFO=<ID=FID,Number=.,Type=ListOfString,Description=\"Feature IDs correspond to ID attributes in a GFF3 file that describe the sequence features (for example genes or mRNAs).\">\n";
        print $file_handle "##INFO=<ID=VE,Number=.,Type=ListOfString,Description=\"Effect that a sequence alteration has on a sequence feature that overlaps it.\",Format=SV|IDX|FT|FID>\n";
    }
    if ($config->{protein_coding_details}) {
        #variant peptide
        print $file_handle "##INFO=<ID=IDX,Number=1,Type=Integer,Description=\"0-based index value that identifies which variant sequence the effect is being described for.\">\n";
        print $file_handle "##INFO=<ID=AmAc,Number=1,Type=String,Description=\"Amino acid translated as result of missense variant.\">\n";
        print $file_handle "##INFO=<ID=FID,Number=.,Type=ListOfString,Description=\"Feature IDs correspond to ID attributes in a GFF3 file that describe the sequence features (for example genes or mRNAs).\">\n";
        print $file_handle "##INFO=<ID=VarPep,Number=.,Type=ListOfString,Description=\"Effect that a sequence alteration has on a sequence feature that overlaps it.\",Format=IDX|AmAc|FID>\n";
        #reference peptide
        print $file_handle "##INFO=<ID=RefPep,Number=1,Type=String,Description=\"Amino acid translated with reference allele.\">\n";
        #sift prediction
        print $file_handle "##INFO=<ID=IDX,Number=1,Type=Integer,Description=\"0-based index value that identifies which variant sequence the effect is being described for.\">\n";
        print $file_handle "##INFO=<ID=QP,Number=1,Type=String,Description=\"Qualitative prediction (either tolerated or deleterious).\">\n";
        print $file_handle "##INFO=<ID=NV,Number=1,Type=String,Description=\"Numerical value which is the normalized probability that the amino acid change is tolerated so scores nearer 0 are more likely to be deleterious.\">\n";
        print $file_handle "##INFO=<ID=FID,Number=.,Type=ListOfString,Description=\"Feature IDs correspond to ID attributes in a GFF3 file that describe the sequence features (for example genes or mRNAs).\">\n";
        print $file_handle "##INFO=<ID=Sift,Number=.,Type=ListOfString,Description=\"Sift prediction.\",Format=IDX|QP|NV|FID>\n";
        #polyphen prediction
        print $file_handle "##INFO=<ID=IDX,Number=1,Type=Integer,Description=\"0-based index value that identifies which variant sequence the effect is being described for.\">\n";
        print $file_handle "##INFO=<ID=QP,Number=1,Type=String,Description=\"Qualitative prediction (one of probably damaging, possibly damaging, benign or unknown).\">\n";
        print $file_handle "##INFO=<ID=NV,Number=1,Type=String,Description=\"Numerical value which is the probability that a substitution is damaging, so values nearer 1 are more confidently predicted to be deleterious.\">\n";
        print $file_handle "##INFO=<ID=FID,Number=.,Type=ListOfString,Description=\"Feature IDs correspond to ID attributes in a GFF3 file that describe the sequence features (for example genes or mRNAs).\">\n";
        print $file_handle "##INFO=<ID=Polyphen,Number=.,Type=ListOfString,Description=\"Polyphen prediction.\",Format=IDX|QP|NV|FID>\n";
    }
    if ($config->{evidence_values}) {
        print $file_handle "##INFO=<ID=E_MO,Number=0,Type=Flag,Description=\"Multiple_observations. The variant has multiple independent dbSNP submissions, i.e. submissions with a different submitter handles or different discovery samples\">\n";
        print $file_handle "##INFO=<ID=E_Freq,Number=0,Type=Flag,Description=\"Frequency. The variant is reported to be polymorphic in at least one sample.\">\n";
        print $file_handle "##INFO=<ID=E_HM,Number=0,Type=Flag,Description=\"HapMap. The variant is polymorphic in at least one HapMap panel (human only).\">\n";
        print $file_handle "##INFO=<ID=E_1000G,Number=0,Type=Flag,Description=\"1000Genomes. The variant was discovered in the 1000 genomes project (human only).\">\n";
        print $file_handle "##INFO=<ID=E_C,Number=0,Type=Flag,Description=\"Cited. dbSNP holds a citation from PubMed for the variant.\">\n";
    }
    if ($config->{sample_data}) {
        print $file_handle "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    }

    my @header_line = ('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO');
    if ($config->{sample_data}) {
        push @header_line, 'FORMAT';
        push @header_line, sort keys %{$config->{sample}};
    };

    print $file_handle join("\t", @header_line), "\n";
}

sub print_vcf_line {
    my $config = shift;
    my $vcf_line = shift;
    my ($info, @info_entries, $sample);
    # info
    if ($vcf_line->{INFO}) {
        for my $key (keys %{$vcf_line->{INFO}}) {
            my $val = $vcf_line->{INFO}->{$key};
            if (ref $val eq 'ARRAY') {
                    push @info_entries, $key . '=' . join(',', @$val);
            } else {
                push @info_entries, $key . '=' . $val;
            }
        }
        $info = "\t" . join(';', @info_entries);
    }
    if ($vcf_line->{FLAG}) {
        $info = $info . ';' . join(';', @{$vcf_line->{FLAG}});
    }

    # sample
    if ($vcf_line->{sample}) {
        $sample = "\t" . $config->{format} . "\t" . join("\t", map {$vcf_line->{sample}->{$_}} sort keys %{$vcf_line->{sample}});
    }
    $vcf_line->{QUAL} = '.';
    $vcf_line->{FILTER} = '.';
    $sample ||= ''; 
    my $line = join("\t", map {$vcf_line->{$_}} qw/CHROM POS ID REF ALT QUAL FILTER/);
    print $file_handle $line . $info . $sample . "\n";
}

sub usage {
    my $usage =<<END;
Usage:
perl dump_vcf.pl [arguments]

Options
-h | --help                 Display this message and quit.
-o | --output_file          Output file.
--species                   Species to use.
-r | --registry             Registry file used to build database connections.
--host                      Database host. 
--user                      Database user.
--port                      Database port.
--db_version                Database version.
--population                Name of population. Dump genotypes for all indivivduals in this population.
--individual                Comma-separated list of individual names. Dump genotypes for those individuals.
--seq_region                Comma-separated list of seq_regions. Dump only for these regions. If not specified all toplevel regions are used.
--allele_freq_in_population Add allele frequencies for this population.
--consequences              Add variation consequence data. 
--protein_coding_details    Add protein function prediction.
--ancestral_allele          Add ancestral allele.
--global_maf                Add global minor allele frequency data.
--evidence_values           Add evidence values supporting a variant as a guide to its potential reliability (available from ensembl 71)
--ref_fasta_file            Location of reference sequence backing the data contained in the VCF file.
--data_source               Data source for VCF dump. Default is ensembl with schema_version and url to species used for data dump. 
END
    print $usage;
}
