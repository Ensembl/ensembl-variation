=head1 LICENSE

 Copyright (c) 1999-2012 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

# EnsEMBL module for Bio::EnsEMBL::Variation::Utils::Sequence
#
#

=head1 NAME

Bio::EnsEMBL::Variation::Utils::VEP - Methods used by the Variant Effect Predictor

=head1 SYNOPSIS

  use Bio::EnsEMBL::Variation::Utils::VEP qw(configure);

  my $config = configure();

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::Utils::VEP;

# module list
use Getopt::Long;
use FileHandle;
use File::Path qw(make_path);
use Storable qw(nstore_fd fd_retrieve freeze thaw);
use Scalar::Util qw(weaken);
use Digest::MD5 qw(md5_hex);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(MAX_DISTANCE_FROM_TRANSCRIPT overlap);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(unambiguity_code);
use Bio::EnsEMBL::Variation::Utils::EnsEMBL2GFF3;
use Bio::EnsEMBL::Variation::StructuralVariationFeature;
use Bio::EnsEMBL::Variation::DBSQL::StructuralVariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::TranscriptStructuralVariation;

# we need to manually include all these modules for caching to work
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::ProteinFeature;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DBSQL::GeneAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::TranslationAdaptor;
use Bio::EnsEMBL::DBSQL::TranscriptAdaptor;
use Bio::EnsEMBL::DBSQL::MetaContainer;
use Bio::EnsEMBL::DBSQL::CoordSystemAdaptor;

use Exporter;
use vars qw(@ISA @EXPORT_OK);
@ISA = qw(Exporter);

# open socket pairs for cross-process comms
use Socket;
socketpair(CHILD, PARENT, AF_UNIX, SOCK_STREAM, PF_UNSPEC) or die "ERROR: Failed to open socketpair: $!";
CHILD->autoflush(1);
PARENT->autoflush(1);

@EXPORT_OK = qw(
    &parse_line
    &vf_to_consequences
    &validate_vf
    &read_cache_info
    &dump_adaptor_cache
    &load_dumped_adaptor_cache
    &load_dumped_variation_cache
    &get_all_consequences
    &get_slice
    &build_slice_cache
    &build_full_cache
    &regions_from_hash
    &get_time
    &debug
    &convert_to_vcf
    &progress
    &end_progress
    @REG_FEAT_TYPES
    @OUTPUT_COLS
    @VEP_WEB_CONFIG
    %FILTER_SHORTCUTS
);

our @OUTPUT_COLS = qw(
    Uploaded_variation
    Location
    Allele
    Gene
    Feature
    Feature_type
    Consequence
    cDNA_position
    CDS_position
    Protein_position
    Amino_acids
    Codons
    Existing_variation
    Extra
);

our @REG_FEAT_TYPES = qw(
    RegulatoryFeature
    MotifFeature
);

our @VEP_WEB_CONFIG = qw(
    format
    check_existing
    coding_only
    core_type
    hgnc
    protein
    hgvs
    terms
    check_frequency
    freq_filter
    freq_gt_lt
    freq_freq
    freq_pop
    sift
    polyphen
    regulatory
);

our %FILTER_SHORTCUTS = (
    upstream => {
        '5KB_upstream_variant' => 1,
        '2KB_upstream_variant' => 1,
    },
    downstream => {
        '5KB_downstream_variant'  => 1,
        '2KB_downstream_variant'  => 1,
        '500B_downstream_variant' => 1,
    },
    utr => {
        '5_prime_UTR_variant' => 1,
        '3_prime_UTR_variant' => 1,
    },
    splice => {
        splice_donor_variant    => 1,
        splice_acceptor_variant => 1,
        splice_region_variant   => 1,
    },
    coding_change => {
        stop_lost            => 1,
        stop_gained          => 1,
        missense_variant     => 1,
        frameshift_variant   => 1,
        inframe_insertion    => 1,
        inframe_deletion     => 1,
    },
    regulatory => {
        regulatory_region_variant => 1,
        TF_binding_site_variant   => 1,
    },
);

# parses a line of input, returns VF object(s)
sub parse_line {
    my $config = shift;
    my $line   = shift;
    
    # find out file format - will only do this on first line
    if(!defined($config->{format}) || (defined($config->{format}) && $config->{format} eq 'guess')) {
        $config->{format} = &detect_format($line);
        debug("Detected format of input file as ", $config->{format}) unless defined($config->{quiet});
        
        # HGVS and ID formats need DB
        die("ERROR: Can't use ".uc($config->{format})." format in offline mode") if $config->{format} =~ /id|hgvs/ && defined($config->{offline});
        
        # force certain options if format is VEP output
        if($config->{format} eq 'vep') {
            $config->{no_consequence} = 1;
            delete $config->{regulatory};
            debug("Forcing no consequence calculation") unless defined($config->{quiet});
        }
    }
    
    # check that format is vcf when using --individual
    die("ERROR: --individual only compatible with VCF input files\n") if defined($config->{individual}) && $config->{format} ne 'vcf';
    
    my $parse_method = 'parse_'.$config->{format};
    $parse_method =~ s/vep_//;
    my $method_ref   = \&$parse_method;
    
    my $vfs = &$method_ref($config, $line);
    
    $vfs = add_lrg_mappings($config, $vfs) if defined($config->{lrg});
    
    return $vfs;
}

# sub-routine to detect format of input
sub detect_format {
    my $line = shift;
    my @data = split /\s+/, $line;
    
    # HGVS: ENST00000285667.3:c.1047_1048insC
    if (
        scalar @data == 1 &&
        $data[0] =~ /^([^\:]+)\:.*?([cgmrp]?)\.?([\*\-0-9]+.*)$/i
    ) {
        return 'hgvs';
    }
    
    # variant identifier: rs123456
    elsif (
        scalar @data == 1
    ) {
        return 'id';
    }
    
    # VCF: 20  14370  rs6054257  G  A  29  0  NS=58;DP=258;AF=0.786;DB;H2  GT:GQ:DP:HQ
    elsif (
        $data[0] =~ /(chr)?\w+/ &&
        $data[1] =~ /^\d+$/ &&
        $data[3] =~ /^[ACGTN-]+$/i &&
        $data[4] =~ /^([\.ACGTN-]+\,?)+$/i
    ) {
        return 'vcf';
    }
    
    # pileup: chr1  60  T  A
    elsif (
        $data[0] =~ /(chr)?\w+/ &&
        $data[1] =~ /^\d+$/ &&
        $data[2] =~ /^[\*ACGTN-]+$/i &&
        $data[3] =~ /^[\*ACGTNRYSWKM\+\/-]+$/i
    ) {
        return 'pileup';
    }
    
    # ensembl: 20  14370  14370  A/G  +
    elsif (
        $data[0] =~ /\w+/ &&
        $data[1] =~ /^\d+$/ &&
        $data[2] =~ /^\d+$/ &&
        $data[3] =~ /[ACGTN-]+\/[ACGTN-]+/i
    ) {
        return 'ensembl';
    }
    
    # vep output: ID    1:142849179     -       -       -       -       INTERGENIC
    elsif (
        $data[0] =~ /\w+/ &&
        $data[1] =~ /^\w+?\:\d+(\-\d+)*$/ &&
        scalar @data == 14
    ) {
        return 'vep';
    }
    
    else {
        die("ERROR: Could not detect input file format\n");
    }
}

# parse a line of Ensembl format input into a variation feature object
sub parse_ensembl {
    my $config = shift;
    my $line = shift;
    
    my ($chr, $start, $end, $allele_string, $strand, $var_name) = split /\s+/, $line;
    
    my $vf = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
        start          => $start,
        end            => $end,
        allele_string  => $allele_string,
        strand         => $strand,
        map_weight     => 1,
        adaptor        => $config->{vfa},
        variation_name => $var_name,
        chr            => $chr,
    });
    
    return [$vf];
}

# parse a line of VCF input into a variation feature object
sub parse_vcf {
    my $config = shift;
    my $line = shift;
    
    my @data = split /\s+/, $line;
    
    # non-variant
    my $non_variant = 0;
    
    if($data[4] eq '.') {
        if(defined($config->{allow_non_variant})) {
            $non_variant = 1;
        }
        else {
            return [];
        }
    }
    
    # get relevant data
    my ($chr, $start, $end, $ref, $alt) = ($data[0], $data[1], $data[1], $data[3], $data[4]);
    
    # some VCF files have a GRCh37 pos defined in GP flag in INFO column
    # if user has requested, we can use that as the position instead
    if(defined $config->{gp}) {
        $chr = undef;
        $start = undef;
        
        foreach my $pair(split /\;/, $data[7]) {
            my ($key, $value) = split /\=/, $pair;
            if($key eq 'GP') {
                ($chr, $start) = split /\:/, $value;
                $end = $start;
            }
        }
        
        unless(defined($chr) and defined($start)) {
            warn "No GP flag found in INFO column" unless defined $config->{quiet};
            return [];
        }
    }
    
    # adjust end coord
    $end += (length($ref) - 1);
    
    # structural variation
    if((defined($data[7]) && $data[7] =~ /SVTYPE/) || $alt =~ /\<|\[|\]|\>/) {
        
        # parse INFO field
        my %info = ();
        
        foreach my $bit(split /\;/, $data[7]) {
            my ($key, $value) = split /\=/, $bit;
            $info{$key} = $value;
        }
        
        # like indels, SVs have the base before included for reference
        $start++;
        
        # work out the end coord
        if(defined($info{END})) {
            $end = $info{END};
        }
        elsif(defined($info{SVLEN})) {
            $end = $start + abs($info{SVLEN}) - 1;
        }
        
        # check for imprecise breakpoints
        my ($min_start, $max_start, $min_end, $max_end);
        
        if(defined($info{CIPOS})) {
            my ($low, $high) = split /\,/, $info{CIPOS};
            $min_start = $start + $low;
            $max_start = $start + $high;
        }
        
        if(defined($info{CIEND})) {
            my ($low, $high) = split /\,/, $info{CIEND};
            $min_end = $end + $low;
            $max_end = $end + $high;
        }
        
        # get type
        my $type = $info{SVTYPE};
        my $so_term;
        
        if(defined($type)) {
            # convert to SO term
            my %terms = (
                INS  => 'insertion',
                DEL  => 'deletion',
                TDUP => 'tandem_duplication',
                DUP  => 'duplication'
            );
            
            $so_term = defined $terms{$type} ? $terms{$type} : $type;
        }
        
        my $svf = Bio::EnsEMBL::Variation::StructuralVariationFeature->new_fast({
            start          => $start,
            inner_start    => $max_start,
            outer_start    => $min_start,
            end            => $end,
            inner_end      => $min_end,
            outer_end      => $max_end,
            strand         => 1,
            adaptor        => $config->{svfa},
            variation_name => $data[2] eq '.' ? undef : $data[2],
            chr            => $chr,
            class_SO_term  => $so_term,
        });
        
        return [$svf];
    }
    
    # normal variation
    else {
        # find out if any of the alt alleles make this an insertion or a deletion
        my ($is_indel, $is_sub, $ins_count, $total_count);
        foreach my $alt_allele(split /\,/, $alt) {
            $is_indel = 1 if $alt_allele =~ /D|I/;
            $is_indel = 1 if length($alt_allele) != length($ref);
            $is_sub = 1 if length($alt_allele) == length($ref);
            $ins_count++ if length($alt_allele) > length($ref);
            $total_count++;
        }
        
        # multiple alt alleles?
        if($alt =~ /\,/) {
            if($is_indel) {
                my @alts;
                
                if($alt =~ /D|I/) {
                    foreach my $alt_allele(split /\,/, $alt) {
                        # deletion (VCF <4)
                        if($alt_allele =~ /D/) {
                            push @alts, '-';
                        }
                        
                        elsif($alt_allele =~ /I/) {
                            $alt_allele =~ s/^I//g;
                            push @alts, $alt_allele;
                        }
                    }
                }
                
                else {
                    $ref = substr($ref, 1) || '-';
                    $start++;
                    
                    foreach my $alt_allele(split /\,/, $alt) {
                        $alt_allele = substr($alt_allele, 1);
                        $alt_allele = '-' if $alt_allele eq '';
                        push @alts, $alt_allele;
                    }
                }
                
                $alt = join "/", @alts;
            }
            
            else {
                # for substitutions we just need to replace ',' with '/' in $alt
                $alt =~ s/\,/\//g;
            }
        }
        
        elsif($is_indel) {
            # deletion (VCF <4)
            if($alt =~ /D/) {
                my $num_deleted = $alt;
                $num_deleted =~ s/\D+//g;
                $end += $num_deleted - 1;
                $alt = "-";
                $ref .= ("N" x ($num_deleted - 1)) unless length($ref) > 1;
            }
            
            # insertion (VCF <4)
            elsif($alt =~ /I/) {
                $ref = '-';
                $alt =~ s/^I//g;
                $start++;
            }
            
            # insertion or deletion (VCF 4+)
            elsif(substr($ref, 0, 1) eq substr($alt, 0, 1)) {
                
                # chop off first base
                $ref = substr($ref, 1) || '-';
                $alt = substr($alt, 1) || '-';
                
                $start++;
            }
        }
        
        # create VF object
        my $vf = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
            start          => $start,
            end            => $end,
            allele_string  => $non_variant ? $ref : $ref.'/'.$alt,
            strand         => 1,
            map_weight     => 1,
            adaptor        => $config->{vfa},
            variation_name => $data[2] eq '.' ? undef : $data[2],
            chr            => $chr,
        });
        
        # flag as non-variant
        $vf->{non_variant} = 1 if $non_variant;
        
        # individuals?
        if(defined($config->{individual})) {
            my @alleles = split /\//, $ref.'/'.$alt;
            
            my @return;
            
            foreach my $ind(keys %{$config->{ind_cols}}) {
                
                # get alleles present in this individual
                my @bits;
                my $gt = (split /\:/, $data[$config->{ind_cols}->{$ind}])[0];
                
                my $phased = ($gt =~ /\|/ ? 1 : 0);
                
                foreach my $bit(split /\||\/|\\/, $gt) {
                    push @bits, $alleles[$bit] unless $bit eq '.';
                }
                
                # shallow copy VF
                my $vf_copy = { %$vf };
                bless $vf_copy, ref($vf);
                
                # get non-refs
                my %non_ref = map {$_ => 1} grep {$_ ne $ref} @bits;
                
                # construct allele_string
                if(scalar keys %non_ref) {
                    $vf_copy->{allele_string} = $ref."/".(join "/", keys %non_ref);
                }
                else {
                    $vf_copy->{allele_string} = $ref;
                    $vf_copy->{non_variant} = 1;
                }
                
                # store phasing info
                $vf_copy->{phased} = defined($config->{phased} ? 1 : $phased);
                
                # store GT
                $vf_copy->{genotype} = \@bits;
                
                # store individual name
                $vf_copy->{individual} = $ind;
                
                push @return, $vf_copy;
            }
            
            return \@return;
        }
        else {
            return [$vf];
        }
    }
}

# parse a line of pileup input into variation feature objects
sub parse_pileup {
    my $config = shift;
    my $line = shift;
    
    my @data = split /\s+/, $line;
    
    # pileup can produce more than one VF per line
    my @return;
    
    # normal variant
    if($data[2] ne "*"){
        my $var;
        
        if($data[3] =~ /^[A|C|G|T]$/) {
            $var = $data[3];
        }
        else {
            ($var = unambiguity_code($data[3])) =~ s/$data[2]//ig;
        }
        
        for my $alt(split //, $var){
            push @return, Bio::EnsEMBL::Variation::VariationFeature->new_fast({
                start          => $data[1],
                end            => $data[1],
                allele_string  => $data[2].'/'.$alt,
                strand         => 1,
                map_weight     => 1,
                adaptor        => $config->{vfa},
                chr            => $data[0],
            });
        }
    }
    
    # in/del
    else {
        my %tmp_hash = map {$_ => 1} split /\//, $data[3];
        my @genotype = keys %tmp_hash;
        
        foreach my $allele(@genotype){
            if(substr($allele,0,1) eq "+") { #ins
                push @return, Bio::EnsEMBL::Variation::VariationFeature->new_fast({
                    start          => $data[1] + 1,
                    end            => $data[1],
                    allele_string  => '-/'.substr($allele, 1),
                    strand         => 1,
                    map_weight     => 1,
                    adaptor        => $config->{vfa},
                    chr            => $data[0],
                });
            }
            elsif(substr($allele,0,1) eq "-"){ #del
                push @return, Bio::EnsEMBL::Variation::VariationFeature->new_fast({
                    start          => $data[1] + 1,
                    end            => $data[1] + length(substr($allele, 1)),
                    allele_string  => substr($allele, 1).'/-',
                    strand         => 1,
                    map_weight     => 1,
                    adaptor        => $config->{vfa},
                    chr            => $data[0],
                });
            }
            elsif($allele ne "*"){
                warn("WARNING: invalid pileup indel genotype: $line\n") unless defined $config->{quiet};
            }
        }
    }
    
    return \@return;
}

# parse a line of HGVS input into a variation feature object
sub parse_hgvs {
    my $config = shift;
    my $line = shift;
    
    my $vf;
    
    # not all hgvs notations are supported yet, so we have to wrap it in an eval
    eval { $vf = $config->{vfa}->fetch_by_hgvs_notation($line, $config->{sa}, $config->{ta}) };
    
    if((!defined($vf) || (defined $@ && length($@) > 1)) && defined($config->{coordinator})) {
        eval { $vf = $config->{vfa}->fetch_by_hgvs_notation($line, $config->{ofsa}, $config->{ofta}) };
    }
    
    if(!defined($vf) || (defined $@ && length($@) > 1)) {
        warn("WARNING: Unable to parse HGVS notation \'$line\'\n$@") unless defined $config->{quiet};
        return [];
    }
    
    # get whole chromosome slice
    my $slice = $vf->slice->adaptor->fetch_by_region($vf->slice->coord_system->name, $vf->slice->seq_region_name);
    
    $vf = $vf->transfer($slice);
    
    # name it after the HGVS
    $vf->{variation_name} = $line;
    
    # add chr attrib
    $vf->{chr} = $vf->slice->seq_region_name;
    
    return [$vf];
}

# parse a variation identifier e.g. a dbSNP rsID
sub parse_id {
    my $config = shift;
    my $line = shift;
    
    my $v_obj = $config->{va}->fetch_by_name($line);
    
    return [] unless defined $v_obj;
    
    my @vfs = @{$v_obj->get_all_VariationFeatures};
    delete $_->{dbID} for @vfs;
    delete $_->{overlap_consequences} for @vfs;
    $_->{chr} = $_->seq_region_name for @vfs;
    
    return \@vfs;
}

# parse a line of VEP output
sub parse_vep {
    my $config = shift;
    my $line = shift;
    
    my @data = split /\t/, $line;
    
    my ($chr, $start, $end) = split /\:|\-/, $data[1];
    $end ||= $start;
    
    # might get allele string from ID
    my $allele_string;
    
    if($data[0] =~ /^\w\_\w\_\w$/) {
        my @split = split /\_/, $data[0];
        $allele_string = $split[-1] if $split[-1] =~ /[ACGTN-]+\/[ACGTN-]+/;
    }
    
    $allele_string ||= 'N/'.($data[6] =~ /intergenic/ ? 'N' : $data[2]);
    
    my $vf = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
        start          => $start,
        end            => $end,
        allele_string  => $allele_string,
        strand         => 1,
        map_weight     => 1,
        adaptor        => $config->{vfa},
        chr            => $chr,
        variation_name => $data[0],
    });
    
    return [$vf];
}



# converts to VCF format
sub convert_to_vcf {
    my $config = shift;
    my $vf = shift;
    
    # look for imbalance in the allele string
    my %allele_lengths;
    my @alleles = split /\//, $vf->allele_string;
    
    foreach my $allele(@alleles) {
        $allele =~ s/\-//g;
        $allele_lengths{length($allele)} = 1;
    }
    
    # in/del/unbalanced
    if(scalar keys %allele_lengths > 1) {
        
        # we need the ref base before the variation
        # default to N in case we can't get it
        my $prev_base = 'N';
        
        unless(defined($config->{cache})) {
            my $slice = $vf->slice->sub_Slice($vf->start - 1, $vf->start - 1);
            $prev_base = $slice->seq if defined($slice);
        }
        
        for my $i(0..$#alleles) {
            $alleles[$i] =~ s/\-//g;
            $alleles[$i] = $prev_base.$alleles[$i];
        }
        
        return [
            $vf->{chr} || $vf->seq_region_name,
            $vf->start - 1,
            $vf->variation_name,
            shift @alleles,
            (join ",", @alleles),
            '.', '.', '.'
        ];
        
    }
    
    # balanced sub
    else {
        return [
            $vf->{chr} || $vf->seq_region_name,
            $vf->start,
            $vf->variation_name,
            shift @alleles,
            (join ",", @alleles),
            '.', '.', '.'
        ];
    }
}


# tries to map a VF to the LRG coordinate system
sub add_lrg_mappings {
    my $config = shift;
    my $vfs = shift;
    
    my @new_vfs;
    
    foreach my $vf(@$vfs) {
        
        # add the unmapped VF to the array
        push @new_vfs, $vf;
        
        # make sure the VF has an attached slice
        $vf->{slice} ||= get_slice($config, $vf->{chr});
        next unless defined($vf->{slice});
        
        # transform LRG <-> chromosome
        my $new_vf = $vf->transform($vf->{slice}->coord_system->name eq 'lrg' ? 'chromosome' : 'lrg');
        
        # add it to the array if transformation worked
        if(defined($new_vf)) {
            
            # update new VF's chr entry
            $new_vf->{chr} = $new_vf->seq_region_name;
            push @new_vfs, $new_vf;
        }
    }
    
    return \@new_vfs;
}


# wrapper for whole_genome_fetch and vf_to_consequences
# takes config and a listref of VFs, returns listref of line hashes for printing
sub get_all_consequences {
    my $config     = shift;
    my $listref    = shift;
   
    if ($config->{extra}) {
        eval "use Plugin qw($config);"
    }
    
    # initialize caches
    $config->{$_.'_cache'} ||= {} for qw(tr rf slice);

    # build hash
    my %vf_hash;
    push @{$vf_hash{$_->{chr}}{int($_->{start} / $config->{chunk_size})}{$_->{start}}}, $_ for grep {!defined($_->{non_variant})} @$listref;
    
    my @non_variant = grep {defined($_->{non_variant})} @$listref;
    debug("Skipping ".(scalar @non_variant)." non-variant loci\n") unless defined($config->{quiet});
    
    # get regions
    my $regions = &regions_from_hash($config, \%vf_hash);
    my $trim_regions = $regions;
    
    # get trimmed regions - allows us to discard out-of-range transcripts
    # when using cache
    #if(defined($config->{cache})) {
    #    my $tmp = $config->{cache};
    #    delete $config->{cache};
    #    $trim_regions = regions_from_hash($config, \%vf_hash);
    #    $config->{cache} = $tmp;
    #}
    
    # prune caches
    prune_cache($config, $config->{tr_cache}, $regions, $config->{loaded_tr});
    prune_cache($config, $config->{rf_cache}, $regions, $config->{loaded_rf});
    
    # get chr list
    my %chrs = map {$_->{chr} => 1} @$listref;
    
    my $fetched_tr_count = 0;
    $fetched_tr_count = fetch_transcripts($config, $regions, $trim_regions)
        unless defined($config->{no_consequences});
    
    my $fetched_rf_count = 0;
    $fetched_rf_count = fetch_regfeats($config, $regions, $trim_regions)
        if defined($config->{regulatory})
        && !defined($config->{no_consequences});
    
    # check we can use MIME::Base64
    if(defined($config->{fork})) {
        eval q{ use MIME::Base64; };
        
        if($@) {
            debug("WARNING: Unable to load MIME::Base64, forking disabled") unless defined($config->{quiet});
            delete $config->{fork};
        }
    }
    
    my (@temp_array, @return, @pids);
    if(defined($config->{fork})) {
        my $size = scalar @$listref;
        
        while(my $tmp_vf = shift @$listref) {
            push @temp_array, $tmp_vf;
            
            # fork
            if(scalar @temp_array >= ($size / $config->{fork}) || scalar @$listref == 0) {
                
                my $pid = fork;
                
                if(!defined($pid)) {
                    debug("WARNING: Failed to fork - will attempt to continue without forking") unless defined($config->{quiet});
                    push @temp_array, @$listref;
                    push @return, @{vf_list_to_cons($config, \@temp_array, $regions)};
                    last;
                }
                elsif($pid) {
                    push @pids, $pid;
                    @temp_array = ();
                }
                elsif($pid == 0) {
                    
                    $config->{forked} = $$;
                    $config->{quiet} = 1;
                    
                    # redirect STDERR to PARENT so we can catch errors
                    *STDERR = *PARENT;
                    
                    my $cons = vf_list_to_cons($config, \@temp_array, $regions);
                    
                    # what we're doing here is sending a serialised hash of the
                    # results through to the parent process through the socket.
                    # This is then thawed by the parent process.
                    # $$, or the PID, is added so that the input can be sorted
                    # back into the correct order for output
                    print PARENT $$." ".encode_base64(freeze($_), "\t")."\n" for @$cons;
                    
                    # some plugins may cache stuff, check for this and try and
                    # reconstitute it into parent's plugin cache
                    foreach my $plugin(@{$config->{plugins}}) {
                        
                        next unless defined($plugin->{has_cache});
                        
                        # delete unnecessary stuff and stuff that can't be serialised
                        delete $plugin->{$_} for qw(config feature_types variant_feature_types version feature_types_wanted variant_feature_types_wanted params);
                        print PARENT $$." PLUGIN ".ref($plugin)." ".encode_base64(freeze($plugin), "\t")."\n";
                    }
                    
                    # we need to tell the parent this child is finished
                    # otherwise it keeps listening
                    print PARENT "DONE $$\n";
                    close PARENT;
                    
                    exit(0);
                }
            }
        }
        
        debug("Calculating consequences") unless defined($config->{quiet});
        
        my $fh = $config->{out_file_handle};
        my $done_processes = 0;
        my $done_vars = 0;
        my $total_size = $size;
        my $pruned_count;
        
        # create a hash to store the returned data by PID
        # this means we can sort it correctly on output
        my %by_pid;
        
        # read child input
        while(<CHILD>) {
            
            # child finished
            if(/^DONE/) {
                $done_processes++;
                last if $done_processes == scalar @pids;
            }
            
            # variant finished / progress indicator
            elsif(/^BUMP/) {
                progress($config, ++$done_vars, $total_size);# if $pruned_count == scalar @pids;
            }
            
            # output
            elsif(/^\-?\d+ /) {
                
                # plugin
                if(/^\-?\d+ PLUGIN/) {
                    
                    m/^(\-?\d+) PLUGIN (\w+) /;
                    my ($pid, $plugin) = ($1, $2);
                    
                    # remove the PID
                    s/^\-?\d+ PLUGIN \w+ //;
                    chomp;
                    
                    my $tmp = thaw(decode_base64($_));
                    
                    next unless defined($plugin);
                    
                    # copy data to parent plugin
                    my ($parent_plugin) = grep {ref($_) eq $plugin} @{$config->{plugins}};
                    
                    next unless defined($parent_plugin);
                    
                    merge_hashes($parent_plugin, $tmp);
                }
                
                else {
                    # grab the PID
                    m/^(\-?\d+)\s/;
                    my $pid = $1;
                    die "ERROR: Could not parse forked PID from line $_" unless defined($pid);
                    
                    # remove the PID
                    s/^\-?\d+\s//;
                    chomp;
                    
                    # decode and thaw "output" from forked process
                    push @{$by_pid{$pid}}, thaw(decode_base64($_));
                }
            }
            
            elsif(/^PRUNED/) {
                s/PRUNED //g;
                chomp;
                $pruned_count++;
                $total_size += $_;
            }
            
            elsif(/^DEBUG/) {
                print STDERR;
            }
            
            # something's wrong
            else {
                # kill the other pids
                kill(15, $_) for @pids;
                die("\nERROR: Forked process failed\n$_\n");
            }
        }
        
        end_progress($config);
        
        debug("Writing output") unless defined($config->{quiet});
        
        waitpid($_, 0) for @pids;
        
        # add the sorted data to the return array
        push @return, @{$by_pid{$_} || []} for @pids;
    }
    
    # no forking
    else {
        push @return, @{vf_list_to_cons($config, $listref, $regions)};
    }
    
    if(defined($config->{debug})) {
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
        debug(
            "LINES ", $config->{line_number},
            "\tMEMORY $tot ", (join " ", @$mem),
            "\tDIFF ", (join " ", @$mem_diff),
            "\tCACHE ", total_size($config->{tr_cache}).
            "\tRF ", total_size($config->{rf_cache}),
            "\tVF ", total_size(\%vf_hash),
        );
        #exit(0) if grep {$_ < 0} @$mem_diff;
    }
    
    return \@return;
}

sub vf_list_to_cons {
    my $config = shift;
    my $listref = shift;
    my $regions = shift;    
    
    # get non-variants
    my @non_variants = grep {$_->{non_variant}} @$listref;
    
    my %vf_hash = ();
    push @{$vf_hash{$_->{chr}}{int($_->{start} / $config->{chunk_size})}{$_->{start}}}, $_ for grep {!defined($_->{non_variant})} @$listref;
    
    # check existing VFs
    &check_existing_hash($config, \%vf_hash) if defined($config->{check_existing});# && scalar @$listref > 10;
    
    # get overlapping SVs
    &check_svs_hash($config, \%vf_hash) if defined($config->{check_svs});
    
    # if we are forked, we can trim off some stuff
    if(defined($config->{forked})) {
        my $tmp = $config->{cache};
        delete $config->{cache};
        $regions = regions_from_hash($config, \%vf_hash);
        $config->{cache} = $tmp;
        
        # prune caches
        my $new_count = 0;
        
        $new_count += prune_cache($config, $config->{tr_cache}, $regions, $config->{loaded_tr});
        $new_count += prune_cache($config, $config->{rf_cache}, $regions, $config->{loaded_rf});
        
        print PARENT "PRUNED $new_count\n";
    }
    
    my @return;
    
    foreach my $chr(sort {($a !~ /^\d+$/ || $b !~ /^\d+/) ? $a cmp $b : $a <=> $b} keys %vf_hash) {
        my $finished_vfs = whole_genome_fetch($config, $chr, \%vf_hash);
        
        # non-variants?
        if(scalar @non_variants) {
            push @$finished_vfs, grep {$_->{chr} eq $chr} @non_variants;
            
            # need to re-sort
            @$finished_vfs = sort {$a->{start} <=> $b->{start} || $a->{end} <=> $b->{end}} @$finished_vfs;
        }
        
        debug("Calculating consequences") unless defined($config->{quiet});
        
        my $vf_count = scalar @$finished_vfs;
        my $vf_counter = 0;
        
        while(my $vf = shift @$finished_vfs) {
            progress($config, $vf_counter++, $vf_count) unless $vf_count == 1;
            
            my $filter_ok = 1;
            
            # filtered output
            if(defined($config->{filter})) {
                $filter_ok = filter_by_consequence($config, $vf);
                $config->{filter_count} += $filter_ok;
            }
            
            # skip filtered lines
            next unless $filter_ok;
            
            # original output
            if(defined($config->{original})) {
                push @return, $vf->{_line};
            }
            
            # GVF output
            elsif(defined($config->{gvf})) {
                $vf->source("User");
                
                $config->{gvf_id} ||= 1;
                
                # get custom annotation
                my $custom_annotation = defined($config->{custom}) ? get_custom_annotation($config, $vf) : {};
                $custom_annotation->{ID} = $config->{gvf_id}++;
                
                
                my $tmp = $vf->to_gvf(
                    include_consequences => defined($config->{no_consequences}) ? 0 : 1,
                    extra_attrs          => $custom_annotation,
                );
                push @return, \$tmp;
            }
            
            # VCF output
            elsif(defined($config->{vcf})) {
                
                # convert to VCF, otherwise get line
                my $line = $config->{format} eq 'vcf' ? [split /\s+/, $vf->{_line}] : convert_to_vcf($config, $vf);
                
                if(!defined($line->[7]) || $line->[7] eq '.') {
                    $line->[7] = '';
                }
                
                # get all the lines the normal way                
                # and process them into VCF-compatible string
                my $string = 'CSQ=';
                
                foreach my $line(@{vf_to_consequences($config, $vf)}) {
                    
                    # use the field list (can be user-defined by setting --fields)
                    for my $col(@{$config->{fields}}) {
                        
                        # skip fields already represented in the VCF
                        next if $col eq 'Uploaded_variation' or $col eq 'Location' or $col eq 'Extra';
                        
                        # search for data in main line hash as well as extra field
                        my $data = defined $line->{$col} ? $line->{$col} : $line->{Extra}->{$col};
                        
                        # "-" means null for everything except the Allele field (confusing...)
                        $data = undef if defined($data) and $data eq '-' and $col ne 'Allele';
                        $data =~ s/\,/\&/g if defined $data;
                        $string .= defined($data) ? $data : '';
                        $string .= '|';
                    }
                    
                    $string =~ s/\|$//;
                    $string .= ',';
                }
                
                $string =~ s/\,$//;
                
                if(!defined($config->{no_consequences}) && $string ne 'CSQ=') {
                    $line->[7] .= ($line->[7] ? ';' : '').$string;
                }
                
                # get custom annotation
                if(defined($config->{custom}) && scalar @{$config->{custom}}) {
                    my $custom_annotation = get_custom_annotation($config, $vf);
                    foreach my $key(keys %{$custom_annotation}) {
                        $line->[7] .= ($line->[7] ? ';' : '').$key.'='.$custom_annotation->{$key};
                    }
                }
                
                my $tmp = join "\t", @$line;
                push @return, \$tmp;
            }
            
            # no consequence output from vep input
            elsif(defined($config->{no_consequences}) && $config->{format} eq 'vep') {
                
                my $line = [split /\s+/, $vf->{_line}];
                
                if($line->[13] eq '-') {
                    $line->[13] = '';
                }
                
                # get custom annotation
                if(defined($config->{custom})) {
                    my $custom_annotation = get_custom_annotation($config, $vf);
                    foreach my $key(keys %{$custom_annotation}) {
                        $line->[13] .= ($line->[13] ? ';' : '').$key.'='.$custom_annotation->{$key};
                    }
                }
                
                my $tmp = join "\t", @$line;
                push @return, \$tmp;
            }
            
            # normal output
            else {
                push @return, @{vf_to_consequences($config, $vf)};
            }
            
            print PARENT "BUMP\n" if defined($config->{forked});
        }
        
        end_progress($config) unless scalar @$listref == 1;
    }
    
    return \@return;
}

# takes a variation feature and returns ready to print consequence information
sub vf_to_consequences {
    my $config = shift;
    my $vf = shift;
    
    # use a different method for SVs
    return svf_to_consequences($config, $vf) if $vf->isa('Bio::EnsEMBL::Variation::StructuralVariationFeature'); 
    
    my @return = ();
    
    # method name for consequence terms
    my $term_method = $config->{terms}.'_term';
        
    # find any co-located existing VFs
    $vf->{existing} ||= find_existing($config, $vf) if defined $config->{check_existing};
    
    # skip based on frequency checks?
    if(defined($config->{check_frequency}) && defined($vf->{existing}) && $vf->{existing} ne '-' && defined($config->{va})) {
        return [] unless grep {$_} map {check_frequencies($config, $_)} reverse split(/\,/, $vf->{existing});
        $vf->{freqs} = $config->{filtered_freqs};
    }
    
    # force empty hash into object's transcript_variations if undefined from whole_genome_fetch
    # this will stop the API trying to go off and fill it again
    $vf->{transcript_variations} ||= {} if defined $config->{whole_genome};
    
    # regulatory stuff
    if(!defined $config->{coding_only} && defined $config->{regulatory}) {
        
        for my $rfv (@{ $vf->get_all_RegulatoryFeatureVariations }) {
            
            my $rf = $rfv->regulatory_feature;
            
            my $base_line = {
                Feature_type => 'RegulatoryFeature',
                Feature      => $rf->stable_id,
            };
            
            if(defined($config->{cell_type}) && scalar(@{$config->{cell_type}})) {
                $base_line->{Extra}->{CELL_TYPE} = join ",",
                    map {$_.':'.$rf->{cell_types}->{$_}}
                    grep {$rf->{cell_types}->{$_}}
                    @{$config->{cell_type}};
                    
                $base_line->{Extra}->{CELL_TYPE} =~ s/\s+/\_/g;
            }
            
            # this currently always returns 'RegulatoryFeature', so we ignore it for now
            #$base_line->{Extra}->{REG_FEAT_TYPE} = $rf->feature_type->name;
            
            for my $rfva (@{ $rfv->get_all_alternate_RegulatoryFeatureVariationAlleles }) {
                
                my $line = init_line($config, $vf, $base_line);
                
                $line->{Allele}         = $rfva->variation_feature_seq;
                $line->{Consequence}    = join ',', 
                    map { $_->$term_method || $_->SO_term } 
                        @{ $rfva->get_all_OverlapConsequences };

                $line = run_plugins($rfva, $line, $config);
                        
                push @return, $line;
            }
        }
        
        for my $mfv (@{ $vf->get_all_MotifFeatureVariations }) {
            
            my $mf = $mfv->motif_feature;
           
            # check that the motif has a binding matrix, if not there's not 
            # much we can do so don't return anything

            next unless defined $mf->binding_matrix;

            my $matrix = $mf->binding_matrix->description.' '.$mf->display_label;
            $matrix =~ s/\s+/\_/g;
            
            my $base_line = {
                Feature_type => 'MotifFeature',
                Feature      => $mf->binding_matrix->name,
                Extra        => {
                    MOTIF_NAME  => $matrix,
                }
            };
            
            if(defined($config->{cell_type}) && scalar(@{$config->{cell_type}})) {
                $base_line->{Extra}->{CELL_TYPE} = join ",",
                    map {$_.':'.$mf->{cell_types}->{$_}}
                    grep {$mf->{cell_types}->{$_}}
                    @{$config->{cell_type}};
                    
                $base_line->{Extra}->{CELL_TYPE} =~ s/\s+/\_/g;
            }
            
            for my $mfva (@{ $mfv->get_all_alternate_MotifFeatureVariationAlleles }) {
                
                my $line = init_line($config, $vf, $base_line);
                
                $line->{Extra}->{MOTIF_POS}          = $mfva->motif_start if defined $mfva->motif_start;
                $line->{Extra}->{HIGH_INF_POS}       = ($mfva->in_informative_position ? 'Y' : 'N');
                
                my $delta = $mfva->motif_score_delta if $mfva->variation_feature_seq =~ /^[ACGT]+$/;

                $line->{Extra}->{MOTIF_SCORE_CHANGE} = sprintf("%.3f", $delta) if defined $delta;

                $line->{Allele}         = $mfva->variation_feature_seq;
                $line->{Consequence}    = join ',', 
                    map { $_->$term_method || $_->SO_term } 
                        @{ $mfva->get_all_OverlapConsequences };

                $line = run_plugins($mfva, $line, $config);
                
                push @return, $line;
            }
        }
    }
    
    # get TVs
    my $tvs = $vf->get_all_TranscriptVariations;
    
    # only most severe
    if(defined($config->{most_severe}) || defined($config->{summary})) {
       
        my $line = init_line($config, $vf);
        
        if(defined($config->{summary})) {
            $line->{Consequence} = join ",", @{$vf->consequence_type($config->{terms}) || $vf->consequence_type};
        }
        else {
            $line->{Consequence} = $vf->display_consequence($config->{terms}) || $vf->display_consequence;
        }
        
        push @return, $line;
    }
    
    # pass a true argument to get_IntergenicVariation to stop it doing a reference allele check
    # (to stay consistent with the rest of the VEP)
    elsif ((my $iv = $vf->get_IntergenicVariation(1)) && !defined($config->{no_intergenic})) {
      
        for my $iva (@{ $iv->get_all_alternate_IntergenicVariationAlleles }) {
            
            my $line = init_line($config, $vf);
           
            $line->{Allele} = $iva->variation_feature_seq;
    
            my $cons = $iva->get_all_OverlapConsequences->[0];

            $line->{Consequence} = $cons->$term_method || $cons->SO_term;

            $line = run_plugins($iva, $line, $config);
                    
            push @return, $line;
        }
    }
    
    # user wants only one conseqeunce per gene
    elsif(defined($config->{per_gene})) {
        
        # sort the TVA objects into a hash by gene
        my %by_gene;
        
        foreach my $tv(@$tvs) {
            next if(defined $config->{coding_only} && !($tv->affects_cds));
            
            my $gene = $tv->transcript->{_gene_stable_id} || $config->{ga}->fetch_by_transcript_stable_id($tv->transcript->stable_id)->stable_id;
            
            push @{$by_gene{$gene}}, @{$tv->get_all_alternate_TranscriptVariationAlleles};
        }
        
        foreach my $gene(keys %by_gene) {
            my ($lowest, $lowest_tva);
            
            # at the moment this means the one that comes out last will be picked
            # if there is more than one TVA with the same rank of consequence
            foreach my $tva(@{$by_gene{$gene}}) {
                foreach my $oc(@{$tva->get_all_OverlapConsequences}) {
                    if(!defined($lowest) || $oc->rank < $lowest) {
                        $lowest = $oc->rank;
                        $lowest_tva = $tva;
                    }
                }
            }
            
            push @return, tva_to_line($config, $lowest_tva);
        }
    }
    
    else {
        foreach my $tv(@$tvs) {
            next if(defined $config->{coding_only} && !($tv->affects_cds));
            
            push @return, map {tva_to_line($config, $_)} @{$tv->get_all_alternate_TranscriptVariationAlleles};
            
            undef $tv->{$_} for keys %$tv;
        }
    }
    
    return \@return;
}

# get consequences for a structural variation feature
sub svf_to_consequences {
    my $config = shift;
    my $svf    = shift;
    
    my @return = ();
    
    my $term_method = $config->{terms}.'_term';
    
    if(defined $config->{whole_genome}) {
        $svf->{transcript_structural_variations} ||= [];
        $svf->{regulation_structural_variations}->{$_} ||= [] for @REG_FEAT_TYPES;
    }
    
    if ((my $iv = $svf->get_IntergenicStructuralVariation(1)) && !defined($config->{no_intergenic})) {
      
        for my $iva (@{ $iv->get_all_alternate_IntergenicStructuralVariationAlleles }) {
            
            my $line = init_line($config, $svf);
           
            $line->{Allele} = '-';
            
            my $cons = $iva->get_all_OverlapConsequences->[0];
            
            $line->{Consequence} = $cons->$term_method || $cons->SO_term;
            
            $line = run_plugins($iva, $line, $config);
            
            push @return, $line;
        }
    }
    
    foreach my $svo(@{$svf->get_all_StructuralVariationOverlaps}) {
        
        next if $svo->isa('Bio::EnsEMBL::Variation::IntergenicStructuralVariation');
        
        my $feature = $svo->feature;
        
        # get feature type
        my $feature_type = (split '::', ref($feature))[-1];
        
        my $base_line = {
            Feature_type     => $feature_type,
            Feature          => $feature->stable_id,
            Allele           => $svf->class_SO_term,
        };
        
        if($svo->isa('Bio::EnsEMBL::Variation::BaseTranscriptVariation')) {
            $base_line->{cDNA_position}    = format_coords($svo->cdna_start, $svo->cdna_end);
            $base_line->{CDS_position}     = format_coords($svo->cds_start, $svo->cds_end);
            $base_line->{Protein_position} = format_coords($svo->translation_start, $svo->translation_end);
        }
        
        foreach my $svoa(@{$svo->get_all_StructuralVariationOverlapAlleles}) {
            my $line = init_line($config, $svf, $base_line);
            
            $line->{Consequence} = join ",",
                #map {s/feature/$feature_type/e; $_}
                map {$_->$term_method}
                sort {$a->rank <=> $b->rank}
                @{$svoa->get_all_OverlapConsequences};
            
            # work out overlap amounts            
            my $overlap_start  = (sort {$a <=> $b} ($svf->start, $feature->start))[-1];
            my $overlap_end    = (sort {$a <=> $b} ($svf->end, $feature->end))[0];
            my $overlap_length = ($overlap_end - $overlap_start) + 1;
            my $overlap_pc     = 100 * ($overlap_length / (($feature->end - $feature->start) + 1));
            
            $line->{Extra}->{OverlapBP} = $overlap_length if $overlap_length > 0;
            $line->{Extra}->{OverlapPC} = sprintf("%.2f", $overlap_pc) if $overlap_pc > 0;
            
            add_extra_fields($config, $line, $svoa);
            
            push @return, $line;
        }
    }
    
    return \@return;
}

# run all of the configured plugins on a VariationFeatureOverlapAllele instance
# and store any results in the provided line hash
sub run_plugins {

    my ($bvfoa, $line_hash, $config) = @_;

    my $skip_line = 0;

    for my $plugin (@{ $config->{plugins} || [] }) {

        # check that this plugin is interested in this type of variation feature

        if ($plugin->check_variant_feature_type(ref $bvfoa->base_variation_feature)) {

            # check that this plugin is interested in this type of feature
            
            if ($plugin->check_feature_type(ref $bvfoa->feature || 'Intergenic')) {
                
                eval {
                    my $plugin_results = $plugin->run($bvfoa, $line_hash);
                    
                    if (defined $plugin_results) {
                        if (ref $plugin_results eq 'HASH') {
                            for my $key (keys %$plugin_results) {
                                $line_hash->{Extra}->{$key} = $plugin_results->{$key};
                            }
                        }
                        else {
                            warn "Plugin '".(ref $plugin)."' did not return a hashref, output ignored!\n";
                        }
                    }
                    else {
                        # if a plugin returns undef, that means it want to filter out this line
                        $skip_line = 1;
                    }
                };
                if ($@) {
                    warn "Plugin '".(ref $plugin)."' went wrong: $@";
                }
                
                # there's no point running any other plugins if we're filtering this line, 
                # because the first plugin to skip the line wins, so we might as well last 
                # out of the loop now and avoid any unnecessary computation

                last if $skip_line;
            }
        }
    }

    return $skip_line ? undef : $line_hash;
}

# turn a TranscriptVariationAllele into a line hash
sub tva_to_line {
    my $config = shift;
    my $tva = shift;
    
    my $tv = $tva->transcript_variation;
    my $t  = $tv->transcript;
    
    # method name for consequence terms
    my $term_method = $config->{terms}.'_term';
    
    my $base_line = {
        Feature_type     => 'Transcript',
        Feature          => (defined $t ? $t->stable_id : undef),
        cDNA_position    => format_coords($tv->cdna_start, $tv->cdna_end),
        CDS_position     => format_coords($tv->cds_start, $tv->cds_end),
        Protein_position => format_coords($tv->translation_start, $tv->translation_end),
        Allele           => $tva->variation_feature_seq,
        Amino_acids      => $tva->pep_allele_string,
        Codons           => $tva->display_codon_allele_string,
        Consequence      => join ",", map {$_->$term_method || $_->SO_term} sort {$a->rank <=> $b->rank} @{$tva->get_all_OverlapConsequences},
    };
    
    my $line = init_line($config, $tva->variation_feature, $base_line);
    
    # HGVS
    if(defined $config->{hgvs}) {my $hgvs_t = $tva->hgvs_transcript;
        my $hgvs_p = $tva->hgvs_protein;
        
        $line->{Extra}->{HGVSc} = $hgvs_t if $hgvs_t;
        $line->{Extra}->{HGVSp} = $hgvs_p if $hgvs_p;
    }
    
    foreach my $tool (qw(SIFT PolyPhen)) {
        my $lc_tool = lc($tool);
        
        if (my $opt = $config->{$lc_tool}) {
            my $want_pred   = $opt =~ /^p/i;
            my $want_score  = $opt =~ /^s/i;
            my $want_both   = $opt =~ /^b/i;
            
            if ($want_both) {
                $want_pred  = 1;
                $want_score = 1;
            }
            
            next unless $want_pred || $want_score;
            
            my $pred_meth   = $lc_tool.'_prediction';
            my $score_meth  = $lc_tool.'_score';
            
            my $pred = $tva->$pred_meth;
            
            if($pred) {
                
                if ($want_pred) {
                    $pred =~ s/\s+/\_/;
                    $line->{Extra}->{$tool} = $pred;
                }
                    
                if ($want_score) {
                    my $score = $tva->$score_meth;
                    
                    if(defined $score) {
                        if($want_pred) {
                            $line->{Extra}->{$tool} .= "($score)";
                        }
                        else {
                            $line->{Extra}->{$tool} = $score;
                        }
                    }
                }
            }
        }
    }
    
    $line = add_extra_fields($config, $line, $tva);
    
    return $line;
}

sub add_extra_fields {
    my $config = shift;
    my $line   = shift;
    my $bvfoa  = shift;
    
    # overlapping SVs
    if(defined $config->{check_svs} && defined $bvfoa->base_variation_feature->{overlapping_svs}) {
        $line->{Extra}->{SV} = $bvfoa->base_variation_feature->{overlapping_svs};
    }
    
    # add transcript-specific fields
    $line = add_extra_fields_transcript($config, $line, $bvfoa) if $bvfoa->isa('Bio::EnsEMBL::Variation::BaseTranscriptVariationAllele');
    
    # run plugins
    $line = run_plugins($bvfoa, $line, $config);
    
    return $line;
}

sub add_extra_fields_transcript {
    my $config = shift;
    my $line = shift;
    my $tva = shift;
    
    my $tv = $tva->base_variation_feature_overlap;
    my $tr = $tva->transcript;
    
    # get gene
    my $gene;
    
    $line->{Gene} = $tr->{_gene_stable_id};
    
    if(!defined($line->{Gene})) {
        $gene = $config->{ga}->fetch_by_transcript_stable_id($tr->stable_id);
        $line->{Gene} = $gene ? $gene->stable_id : '-';
    }
    
    # exon/intron numbers
    if ($config->{numbers}) {
        $line->{Extra}->{EXON} = $tv->exon_number if defined $tv->exon_number;
        $line->{Extra}->{INTRON} = $tv->intron_number if defined $tv->intron_number;
    }

    if ($config->{domains}) {
        my $feats = $tv->get_overlapping_ProteinFeatures;

        my @strings;

        for my $feat (@$feats) {
            my $label = $feat->analysis->display_label.':'.$feat->hseqname;

            # replace any special characters
            $label =~ s/[\s;=]/_/g;

            push @strings, $label;
        }

        $line->{Extra}->{DOMAINS} = join ',', @strings if @strings;
    }
    
    # distance to transcript
    if($line->{Consequence} =~ /(up|down)stream/i) {
        $line->{Extra}->{DISTANCE} = $tv->distance_to_transcript;
    }

    # HGNC
    if(defined $config->{hgnc}) {
        my $hgnc;
        $hgnc = $tr->{_gene_hgnc};
        
        if(!defined($hgnc)) {
            if(!defined($gene)) {
                $gene = $config->{ga}->fetch_by_transcript_stable_id($tr->stable_id);
            }
            
            my @entries = grep {$_->database eq 'HGNC'} @{$gene->get_all_DBEntries()};
            if(scalar @entries) {
                $hgnc = $entries[0]->display_id;
            }
        }
        
        $hgnc = undef if defined($hgnc) && $hgnc eq '-';
        
        $line->{Extra}->{HGNC} = $hgnc if defined($hgnc);
    }
    
    # CCDS
    if(defined($config->{ccds})) {
        my $ccds = $tr->{_ccds};
        
        if(!defined($ccds)) {
            my @entries = grep {$_->database eq 'CCDS'} @{$tr->get_all_DBEntries};
            $ccds = $entries[0]->display_id if scalar @entries;
        }
        
        $ccds = undef if defined($ccds) && $ccds eq '-';
        
        $line->{Extra}->{CCDS} = $ccds if defined($ccds);
    }   
    
    # refseq xref
    if(defined($config->{xref_refseq})) {
        my $refseq = $tr->{_refseq};
        
        if(!defined($refseq)) {
            my @entries = grep {$_->database eq 'RefSeq_mRNA'} @{$tr->get_all_DBEntries};
            if(scalar @entries) {
                $refseq = join ",", map {$_->display_id."-".$_->database} @entries;
            }
        }
        
        $refseq = undef if defined($refseq) && $refseq eq '-';
        
        $line->{Extra}->{RefSeq} = $refseq if defined($refseq);
    }
    
    # protein ID
    if(defined $config->{protein}) {
        my $protein = $tr->{_protein};
        
        if(!defined($protein)) {
            $protein = $tr->translation->stable_id if defined($tr->translation);
        }
        
        $protein = undef if defined($protein) && $protein eq '-';
        
        $line->{Extra}->{ENSP} = $protein if defined($protein);
    }
    
    # canonical transcript
    if(defined $config->{canonical}) {
        $line->{Extra}->{CANONICAL} = 'YES' if $tr->is_canonical;
    }
    
    return $line;
}

# initialize a line hash
sub init_line {
    my $config = shift;
    my $vf = shift;
    my $base_line = shift;
    
    my $line = {
        Uploaded_variation  => $vf->variation_name,
        Location            => ($vf->{chr} || $vf->seq_region_name).':'.format_coords($vf->start, $vf->end),
        Existing_variation  => $vf->{existing},
        Extra               => {},
    };
    
    # add custom info
    if(defined($config->{custom}) && scalar @{$config->{custom}}) {
        # merge the custom hash with the extra hash
        my $custom = get_custom_annotation($config, $vf);
        
        for my $key (keys %$custom) {
            $line->{Extra}->{$key} = $custom->{$key};
        }
    }
    
    # individual?
    $line->{Extra}->{IND} = $vf->{individual} if defined($vf->{individual});
    
    # frequencies?
    $line->{Extra}->{FREQS} = join ",", @{$vf->{freqs}} if defined($vf->{freqs});
    
    # gmaf?
    $line->{Extra}->{GMAF} = $vf->{gmaf} if defined($config->{gmaf}) && defined($vf->{gmaf});
    
    # copy entries from base_line
    if(defined($base_line)) {
        $line->{$_} = $base_line->{$_} for keys %$base_line;
    }
    
    return $line;
}


# get custom annotation for a single VF
sub get_custom_annotation {
    my $config = shift;
    my $vf = shift;
    my $cache = shift;
    
    return $vf->{custom} if defined($vf->{custom});
    
    my $annotation = {};
    
    my $chr = $vf->{chr};
    
    if(!defined($cache)) {
        # spoof regions
        my $regions;
        $regions->{$chr} = [$vf->{start}.'-'.$vf->{end}];
        $cache = cache_custom_annotation($config, $regions, $chr);
    }
    
    foreach my $custom(@{$config->{custom}}) {
        next unless defined($cache->{$chr}->{$custom->{name}});
        
        # exact type must match coords of variant exactly
        if($custom->{type} eq 'exact') {
            foreach my $feature(values %{$cache->{$chr}->{$custom->{name}}->{$vf->{start}}}) {
                
                next unless
                    $feature->{chr}   eq $chr &&
                    $feature->{start} eq $vf->{start} &&
                    $feature->{end}   eq $vf->{end};
                    
                $annotation->{$custom->{name}} .= $feature->{name}.',';
            }
        }
        
        # overlap type only needs to overlap, but we need to search the whole range
        elsif($custom->{type} eq 'overlap') {
            foreach my $pos(keys %{$cache->{$chr}->{$custom->{name}}}) {
                foreach my $feature(values %{$cache->{$chr}->{$custom->{name}}->{$pos}}) {
                    
                    next unless
                        $feature->{chr}   eq $chr &&
                        $feature->{end}   >= $vf->{start} &&
                        $feature->{start} <= $vf->{end};
                        
                    $annotation->{$custom->{name}} .= $feature->{name}.',';
                }
            }
        }
        
        # trim off trailing commas
        $annotation->{$custom->{name}} =~ s/\,$//g if defined($annotation->{$custom->{name}});
    }
    
    return $annotation;
}

# decides whether to print a VF based on user defined consequences
sub filter_by_consequence {
    my $config = shift;
    my $vf = shift;
    my $filters = $config->{filter};
    
    # find it if we only have "no"s
    my $only_nos = 0;
    $only_nos = 1 if (sort {$a <=> $b} values %$filters)[-1] == 0;
    
    my ($yes, $no) = (0, 0);
    
    # get all consequences across all term types
    my @types = ('SO', 'display');
    
    my @cons;
    push @cons, @{$vf->consequence_type($_)} for @types;
    
    # add regulatory consequences
    if(defined($config->{regulatory})) {
        foreach my $term_type(@types) {
            my $term_method = $term_type.'_term';
            
            for my $rfv (@{ $vf->get_all_RegulatoryFeatureVariations }) {
                for my $rfva(@{$rfv->get_all_alternate_RegulatoryFeatureVariationAlleles}) {
                    push @cons, map {$_->$term_method} @{ $rfva->get_all_OverlapConsequences };
                }
            }   
            for my $mfv (@{ $vf->get_all_MotifFeatureVariations }) {
                for my $mfva(@{$mfv->get_all_alternate_MotifFeatureVariationAlleles}) {
                    push @cons, map {$_->$term_method} @{ $mfva->get_all_OverlapConsequences };
                }
            }
        }
    }
    
    foreach my $con(grep {defined($_) && defined($filters->{$_})} @cons) {
        if($filters->{$con} == 1) {
            $yes = 1;
        }
        else {
            $no = 1;
        }
    }
    
    # check special case, coding
    if(defined($filters->{coding})) {
        if(grep {$_->affects_cds} @{$vf->get_all_TranscriptVariations}) {
            if($filters->{coding} == 1) {
                $yes = 1;
            }
            else {
                $no = 1;
            }
        }
    }
    
    my $ok = 0;
    if($only_nos) {
        $ok = 1 if !$no;
    }
    else {
        $ok = 1 if $yes && !$no;
    }
    
    return $ok;
}


# takes VFs created from input, fixes and checks various things
sub validate_vf {
    my $config = shift;
    my $vf = shift;
    
    # user specified chr skip list
    return 0 if defined($config->{chr}) && !$config->{chr}->{$vf->{chr}};
    
    # fix inputs
    $vf->{chr} =~ s/chr//ig unless $vf->{chr} =~ /^chromosome$/i;
    $vf->{chr} = 'MT' if $vf->{chr} eq 'M';
    $vf->{strand} ||= 1;
    $vf->{strand} = ($vf->{strand} =~ /\-/ ? "-1" : "1");
    
    # sanity checks
    unless($vf->{start} =~ /^\d+$/ && $vf->{end} =~ /^\d+$/) {
      warn("WARNING: Start ".$vf->{start}." or end ".$vf->{end}." coordinate invalid on line ".$config->{line_number}."\n") unless defined $config->{quiet};
      return 0;
    }
    
    # structural variation?
    return validate_svf($config, $vf) if $vf->isa('Bio::EnsEMBL::Variation::StructuralVariationFeature');
    
    # uppercase allele string
    $vf->{allele_string} =~ tr/[a-z]/[A-Z]/;
    
    unless($vf->{allele_string} =~ /([ACGT-]+\/*)+/) {
      warn("WARNING: Invalid allele string ".$vf->{allele_string}." on line ".$config->{line_number}." or possible parsing error\n") unless defined $config->{quiet};
      return 0;
    }
    
    # insertion should have start = end + 1
    if($vf->{allele_string} =~ /^\-\// && $vf->{start} != $vf->{end} + 1) {
        warn(
            "WARNING: Alleles look like an insertion (".
            $vf->{allele_string}.
            ") but coordinates are not start = end + 1 (START=".
            $vf->{start}.", END=".$vf->{end}.
            ") on line ".$config->{line_number}."\n"
        ) unless defined($config->{quiet});
        return 0;
    }
    
    # check length of reference matches seq length spanned
    my @alleles = split /\//, $vf->{allele_string};
    my $ref_allele = shift @alleles;
    my $tmp_ref_allele = $ref_allele;
    $tmp_ref_allele =~ s/\-//g;
    
    #if(($vf->{end} - $vf->{start}) + 1 != length($tmp_ref_allele)) {
    #    warn(
    #        "WARNING: Length of reference allele (".$ref_allele.
    #        " length ".length($tmp_ref_allele).") does not match co-ordinates ".$vf->{start}."-".$vf->{end}.
    #        " on line ".$config->{line_number}
    #    ) unless defined($config->{quiet});
    #    return 0;
    #}
    
    # flag as unbalanced
    foreach my $allele(@alleles) {
        $allele =~ s/\-//g;
        $vf->{indel} = 1 unless length($allele) == length($tmp_ref_allele);
    }
    
    # check reference allele if requested
    if(defined $config->{check_ref}) {
        my $ok = 0;
        my $slice_ref_allele;
        
        # insertion, therefore no ref allele to check
        if($ref_allele eq '-') {
            $ok = 1;
        }
        else {
            my $slice_ref = $vf->{slice}->sub_Slice($vf->{start}, $vf->{end}, $vf->{strand});
            
            if(!defined($slice_ref)) {
                warn "WARNING: Could not fetch sub-slice from ".$vf->{start}."\-".$vf->{end}."\(".$vf->{strand}."\) on line ".$config->{line_number} unless defined $config->{quiet};
            }
            
            else {
                $slice_ref_allele = $slice_ref->seq;
                $ok = ($slice_ref_allele eq $ref_allele ? 1 : 0);
            }
        }
        
        if(!$ok) {
            warn
                "WARNING: Specified reference allele $ref_allele ",
                "does not match Ensembl reference allele",
                ($slice_ref_allele ? " $slice_ref_allele" : ""),
                " on line ".$config->{line_number} unless defined $config->{quiet};
            return 0;
        }
    }
    
    return 1;
}


# validate a structural variation
sub validate_svf {
    my $config = shift;
    my $svf = shift;
    
    return 1;
}


# takes a hash of VFs and fetches consequences by pre-fetching overlapping transcripts
# from database and/or cache
sub whole_genome_fetch {
    my $config = shift;
    my $chr = shift;
    my $vf_hash = shift;
    
    my (%vf_done, @finished_vfs, %seen_rfs);
    
    if(defined($config->{offline}) && !-e $config->{dir}.'/'.$chr) {
        debug("No cache found for chromsome $chr") unless defined($config->{quiet});
        
        foreach my $chunk(keys %{$vf_hash->{$chr}}) {
            foreach my $pos(keys %{$vf_hash->{$chr}{$chunk}}) {
                push @finished_vfs, @{$vf_hash->{$chr}{$chunk}{$pos}};
            }
        }
        
        return \@finished_vfs;
    }
    
    my $slice_cache = $config->{slice_cache};
    build_slice_cache($config, $config->{tr_cache}) unless defined($slice_cache->{$chr});
    
    debug("Analyzing chromosome $chr") unless defined($config->{quiet});
    
    # custom annotations
    whole_genome_fetch_custom($config, $vf_hash, $chr) if defined($config->{custom});
    
    # split up normal variations from SVs
    my ($tmp_vf_hash, @svfs);
    
    foreach my $chunk(keys %{$vf_hash->{$chr}}) {
        foreach my $pos(keys %{$vf_hash->{$chr}{$chunk}}) {
            foreach my $vf(@{$vf_hash->{$chr}{$chunk}{$pos}}) {
                if($vf->isa('Bio::EnsEMBL::Variation::StructuralVariationFeature')) {
                    push @svfs, $vf;
                }
                else {
                    push @{$tmp_vf_hash->{$chr}{$chunk}{$pos}}, $vf;
                }
            }
        }
    }
    
    $vf_hash = $tmp_vf_hash;
    
    # transcript annotations
    whole_genome_fetch_transcript($config, $vf_hash, $chr)
        unless defined($config->{no_consequences});
    
    # regulatory annotations
    whole_genome_fetch_reg($config, $vf_hash, $chr)
        if defined($config->{regulatory})
        && !defined($config->{no_consequences});
        
    # structural variations
    @finished_vfs = @{whole_genome_fetch_sv($config, \@svfs, $chr)}
        if scalar @svfs;
    
    # sort results into @finished_vfs array
    foreach my $chunk(keys %{$vf_hash->{$chr}}) {
        foreach my $pos(keys %{$vf_hash->{$chr}{$chunk}}) {
            
            # pinch slice from slice cache if we don't already have it
            $_->{slice} ||= $slice_cache->{$chr} for @{$vf_hash->{$chr}{$chunk}{$pos}};
            
            if(defined($config->{regulatory})) {
                foreach my $type(@REG_FEAT_TYPES) {
                    $_->{regulation_variations}->{$type} ||= [] for @{$vf_hash->{$chr}{$chunk}{$pos}};
                }
            }
            
            if(defined($config->{custom})) {
                $_->{custom} ||= {} for @{$vf_hash->{$chr}{$chunk}{$pos}};
            }
            
            $_->{transcript_variations} ||= {} for @{$vf_hash->{$chr}{$chunk}{$pos}};
            
            # add to final array
            push @finished_vfs, @{$vf_hash->{$chr}{$chunk}{$pos}};
        }
    }
    
    # sort
    @finished_vfs = sort {$a->{start} <=> $b->{start} || $a->{end} <=> $b->{end}} @finished_vfs;
    
    # clean hash
    delete $vf_hash->{$chr};
    
    return \@finished_vfs;
}

sub whole_genome_fetch_custom {
    my $config = shift;
    my $vf_hash = shift;
    my $chr = shift;
    
    return unless scalar @{$config->{custom}};
    
    # create regions based on VFs instead of chunks
    my $tmp_regions;
    
    foreach my $chunk(keys %{$vf_hash->{$chr}}) {
        foreach my $pos(keys %{$vf_hash->{$chr}{$chunk}}) {
            foreach my $vf(@{$vf_hash->{$chr}{$chunk}{$pos}}) {
                push @{$tmp_regions->{$chr}}, ($vf->{start}-1).'-'.($vf->{end}+1);
            }
        }
    }
    
    return unless defined($tmp_regions->{$chr});
    
    # cache annotations
    my $annotation_cache = cache_custom_annotation($config, $tmp_regions, $chr);
    
    # count and report
    my $total_annotations = 0;
    $total_annotations += scalar keys %{$annotation_cache->{$chr}->{$_}} for keys %{$annotation_cache->{$chr}};
    debug("Retrieved $total_annotations custom annotations (", (join ", ", map {(scalar keys %{$annotation_cache->{$chr}->{$_}}).' '.$_} keys %{$annotation_cache->{$chr}}), ")");
    
    # compare annotations to variations in hash
    debug("Analyzing custom annotations") unless defined($config->{quiet});
    my $total = scalar keys %{$vf_hash->{$chr}};
    my $i = 0;
    
    foreach my $chunk(keys %{$vf_hash->{$chr}}) {
        progress($config, $i++, $total);
        
        foreach my $pos(keys %{$vf_hash->{$chr}{$chunk}}) {
            foreach my $vf(@{$vf_hash->{$chr}{$chunk}{$pos}}) {
                
                $vf->{custom} = get_custom_annotation($config, $vf, $annotation_cache);
            }
        }
    }
    
    end_progress($config);
}

sub whole_genome_fetch_transcript {
    my $config = shift;
    my $vf_hash = shift;
    my $chr = shift;
    
    my $tr_cache = $config->{tr_cache};
    my $slice_cache = $config->{slice_cache};
    
    my $up_down_size = MAX_DISTANCE_FROM_TRANSCRIPT;
    
    # check we have defined regions
    return unless defined($vf_hash->{$chr}) && defined($tr_cache->{$chr});
    
    # copy slice from transcript to slice cache
    $slice_cache = build_slice_cache($config, $tr_cache) unless defined($slice_cache->{$chr});
    
    debug("Analyzing variants") unless defined($config->{quiet});
    
    my $tr_counter = 0;
    my $tr_count   = scalar @{$tr_cache->{$chr}};
    
    while($tr_counter < $tr_count) {
        
        progress($config, $tr_counter, $tr_count);
        print PARENT "BUMP\n" if defined($config->{forked});
        
        my $tr = $tr_cache->{$chr}->[$tr_counter++];
        
        # do each overlapping VF
        my $s = $tr->start - $up_down_size;
        my $e = $tr->end + $up_down_size;
        
        # get the chunks this transcript overlaps
        my %chunks;
        $chunks{$_} = 1 for (int($s/$config->{chunk_size})..int($e/$config->{chunk_size}));
        map {delete $chunks{$_} unless defined($vf_hash->{$chr}{$_})} keys %chunks;
        
        # pointer to previous VF
        # used to tell plugins this is the last variant analysed in this transcript
        my $previous_vf;
        
        foreach my $chunk(keys %chunks) {
            foreach my $vf(
                grep {$_->{start} <= $e && $_->{end} >= $s}
                map {@{$vf_hash->{$chr}{$chunk}{$_}}}
                keys %{$vf_hash->{$chr}{$chunk}}
            ) {
                # pinch slice from slice cache if we don't already have it
                $vf->{slice} ||= $slice_cache->{$chr};
                $vf->{slice} = $slice_cache->{$chr} if defined($vf->{slice}->{is_fake});
                
                my $tv = Bio::EnsEMBL::Variation::TranscriptVariation->new(
                    -transcript        => $tr,
                    -variation_feature => $vf,
                    -adaptor           => $config->{tva},
                    -no_ref_check      => 1,
                    -no_transfer       => 1
                );
                
                # prefetching stuff here prevents doing loads at the
                # end and makes progress reporting more useful
                $tv->_prefetch_for_vep;
                
                $vf->add_TranscriptVariation($tv);
                
                # cache VF on the transcript if it is an unbalanced sub
                push @{$tr->{indels}}, $vf if defined($vf->{indel});
                
                if(defined($config->{individual})) {
                    
                    # store VF on transcript, weaken reference to avoid circularity
                    push @{$tr->{vfs}}, $vf;
                    weaken($tr->{vfs}->[-1]);
                    
                    delete $previous_vf->{last_in_transcript}->{$tr->stable_id};
                    $vf->{last_in_transcript}->{$tr->stable_id} = 1;
                }
                
                $previous_vf = $vf;
            }
        }
    }
    
    end_progress($config);
}

sub whole_genome_fetch_reg {
    my $config = shift;
    my $vf_hash = shift;
    my $chr = shift;
    
    my $rf_cache = $config->{rf_cache};
    my $slice_cache = $config->{slice_cache};
    
    foreach my $type(keys %{$rf_cache->{$chr}}) {
        debug("Analyzing ".$type."s") unless defined($config->{quiet});
        
        my $constructor = 'Bio::EnsEMBL::Variation::'.$type.'Variation';
        
        my $rf_counter = 0;
        my $rf_count = scalar @{$rf_cache->{$chr}->{$type}};
        
        while($rf_counter < $rf_count) {
            
            progress($config, $rf_counter, $rf_count);
            print PARENT "BUMP\n" if defined($config->{forked});
            
            my $rf = $rf_cache->{$chr}->{$type}->[$rf_counter++];
            
            # do each overlapping VF
            my $s = $rf->{start};
            my $e = $rf->{end};
            
            # get the chunks this transcript overlaps
            my %chunks;
            $chunks{$_} = 1 for (int($s/$config->{chunk_size})..int($e/$config->{chunk_size}));
            map {delete $chunks{$_} unless defined($vf_hash->{$chr}{$_})} keys %chunks;
            
            foreach my $chunk(keys %chunks) {
                foreach my $vf(
                    grep {$_->{start} <= $e && $_->{end} >= $s}
                    map {@{$vf_hash->{$chr}{$chunk}{$_}}}
                    keys %{$vf_hash->{$chr}{$chunk}}
                ) {
                    push @{$vf->{regulation_variations}->{$type}}, $constructor->new(
                        -variation_feature  => $vf,
                        -feature            => $rf,
                        -no_ref_check       => 1,
                        -no_transfer        => 1
                    );
                }
            }
        }
        
        end_progress($config);
    }
}

sub whole_genome_fetch_sv {
    my $config = shift;
    my $svfs = shift;
    my $chr = shift;
    
    my $tr_cache = $config->{tr_cache};
    my $rf_cache = $config->{rf_cache};
    my $slice_cache = $config->{slice_cache};
    
    debug("Analyzing structural variations") unless defined($config->{quiet});
    
    my($i, $total) = (0, scalar @$svfs);
    
    my @finished_vfs;
    
    foreach my $svf(@$svfs) {
        progress($config, $i++, $total);
        
        my %done_genes = ();
        
        if(defined($tr_cache->{$chr})) {
            foreach my $tr(grep {overlap($_->{start} - MAX_DISTANCE_FROM_TRANSCRIPT, $_->{end} + MAX_DISTANCE_FROM_TRANSCRIPT, $svf->{start}, $svf->{end})} @{$tr_cache->{$chr}}) {
                my $svo = Bio::EnsEMBL::Variation::TranscriptStructuralVariation->new(
                    -transcript                   => $tr,
                    -structural_variation_feature => $svf,
                    -no_transfer                  => 1
                );
                
                $svf->add_TranscriptStructuralVariation($svo);
            }
        }
        
        $svf->{transcript_structural_variations} ||= {};
        
        # do regulatory features
        if(defined($config->{regulatory}) && defined($rf_cache->{$chr})) {
            foreach my $rf_type(qw/RegulatoryFeature/) {#keys %{$rf_cache->{$chr}}) {
                foreach my $rf(grep {$_->{start} <= $svf->{end} && $_->end >= $svf->{end}} @{$rf_cache->{$chr}->{$rf_type}}) {
                    my $svo = Bio::EnsEMBL::Variation::StructuralVariationOverlap->new(
                        -feature                      => $rf,
                        -structural_variation_feature => $svf,
                        -no_transfer                  => 1
                    );
                    
                    push @{$svf->{regulation_structural_variations}->{$rf_type}}, $svo;
                }
                
                $svf->{regulation_structural_variations}->{$rf_type} ||= [];
            }
        }
        
        # sort them
        #$svf->_sort_svos;
        push @finished_vfs, $svf;
    }
    
    end_progress($config);
    
    return \@finished_vfs;
}

# retrieves transcripts given region list
sub fetch_transcripts {
    my $config = shift;
    my $regions = shift;
    my $trim_regions = shift;
    
    my $tr_cache = $config->{tr_cache};
    my $slice_cache = $config->{slice_cache};
    
    my ($count_from_mem, $count_from_db, $count_from_cache, $count_duplicates, $count_trimmed) = (0, 0, 0, 0, 0);
    
    my %seen_trs;
    
    $count_from_mem = 0;
    my $region_count = 0;
    foreach my $chr(keys %{$regions}) {
        $count_from_mem += scalar @{$tr_cache->{$chr}} if defined($tr_cache->{$chr}) && ref($tr_cache->{$chr}) eq 'ARRAY';
        $region_count += scalar @{$regions->{$chr}};
    }
    
    my $counter;
    
    debug("Reading transcript data from cache and/or database") unless defined($config->{quiet});
    
    foreach my $chr(keys %{$regions}) {
        foreach my $region(sort {(split /\-/, $a)[0] <=> (split /\-/, $b)[1]} @{$regions->{$chr}}) {
            progress($config, $counter++, $region_count);
            
            # skip regions beyond the end of the chr
            next if defined($slice_cache->{$chr}) && (split /\-/, $region)[0] > $slice_cache->{$chr}->length;
            
            next if defined($config->{loaded_tr}->{$chr}->{$region});
            
            # force quiet so other methods don't mess up the progress bar
            my $quiet = $config->{quiet};
            $config->{quiet} = 1;
            
            # try and load cache from disk if using cache
            my $tmp_cache;
            if(defined($config->{cache})) {
                #$tmp_cache = (
                #    defined($config->{tabix}) ?
                #    load_dumped_transcript_cache_tabix($config, $chr, $region, $trim_regions) :
                #    load_dumped_transcript_cache($config, $chr, $region)
                #);
                $tmp_cache = load_dumped_transcript_cache($config, $chr, $region);
                $count_from_cache += scalar @{$tmp_cache->{$chr}} if defined($tmp_cache->{$chr});
                $config->{loaded_tr}->{$chr}->{$region} = 1;
            }
            
            # no cache found on disk or not using cache
            if(!defined($tmp_cache->{$chr})) {
                
                if(defined($config->{offline})) {
                    # restore quiet status
                    $config->{quiet} = $quiet;
                    
                    debug("WARNING: Could not find cache for $chr\:$region") unless defined($config->{quiet});
                    next;
                }
                
                # spoof temporary region hash
                my $tmp_hash;
                push @{$tmp_hash->{$chr}}, $region;
                
                $tmp_cache = cache_transcripts($config, $tmp_hash);
                
                # make it an empty arrayref that gets cached
                # so we don't get confused and reload next time round
                $tmp_cache->{$chr} ||= [];
                
                $count_from_db += scalar @{$tmp_cache->{$chr}};
                
                # dump to disk if writing to cache
                (defined($config->{tabix}) ? dump_transcript_cache_tabix($config, $tmp_cache, $chr, $region) : dump_transcript_cache($config, $tmp_cache, $chr, $region)) if defined($config->{write_cache});
                
                $config->{loaded_tr}->{$chr}->{$region} = 1;
            }
            
            # add loaded transcripts to main cache
            if(defined($tmp_cache->{$chr})) {
                while(my $tr = shift @{$tmp_cache->{$chr}}) {
                    
                    # track already added transcripts by dbID
                    my $dbID = $tr->dbID;
                    if($seen_trs{$dbID}) {
                        $count_duplicates++;
                        next;
                    }
                    
                    # trim out?
                    #if(defined($trim_regions) && defined($trim_regions->{$chr})) {
                    #    my $tmp_count = scalar grep {
                    #        overlap(
                    #            (split /\-/, $_)[0], (split /\-/, $_)[1],
                    #            $tr->{start}, $tr->{end}
                    #        )
                    #    } @{$trim_regions->{$chr}};
                    #    
                    #    if(!$tmp_count) {
                    #        $count_trimmed++;
                    #        next;
                    #    }
                    #}
                    
                    $seen_trs{$dbID} = 1;
                    
                    push @{$tr_cache->{$chr}}, $tr;
                }
            }
            
            $tr_cache->{$chr} ||= [];
            
            undef $tmp_cache;
            
            # restore quiet status
            $config->{quiet} = $quiet;
            
            # build slice cache
            $slice_cache = build_slice_cache($config, $tr_cache) unless defined($slice_cache->{$chr});
        }
    }
    
    end_progress($config);
    
    my $tr_count = 0;
    $tr_count += scalar @{$tr_cache->{$_}} for keys %$tr_cache;
    
    debug("Retrieved $tr_count transcripts ($count_from_mem mem, $count_from_cache cached, $count_from_db DB, $count_duplicates duplicates)") unless defined($config->{quiet});
    
    return $tr_count;
}

sub fetch_regfeats {
    my $config = shift;
    my $regions = shift;
    my $trim_regions = shift;
    
    my $rf_cache = $config->{rf_cache};
    my $slice_cache = $config->{slice_cache};
    
    my ($count_from_mem, $count_from_db, $count_from_cache, $count_duplicates, $count_trimmed) = (0, 0, 0, 0, 0);
    
    my %seen_rfs;
    
    $count_from_mem = 0;
    my $region_count = 0;
    
    foreach my $chr(keys %$regions) {
        if(defined($rf_cache->{$chr}) && ref($rf_cache->{$chr}) eq 'HASH') {
            $count_from_mem += scalar @{$rf_cache->{$chr}->{$_}} for keys %{$rf_cache->{$chr}};
        }
        $region_count += scalar @{$regions->{$chr}};
    }
    
    my $counter = 0;
    
    debug("Reading regulatory data from cache and/or database") unless defined($config->{quiet});
    
    foreach my $chr(keys %$regions) {
        foreach my $region(sort {(split /\-/, $a)[0] cmp (split /\-/, $b)[1]} @{$regions->{$chr}}) {
            progress($config, $counter++, $region_count);
            
            next if defined($config->{loaded_rf}->{$chr}->{$region});
            
            # skip regions beyond the end of the chr
            next if defined($slice_cache->{$chr}) && (split /\-/, $region)[0] > $slice_cache->{$chr}->length;
            
            # force quiet so other methods don't mess up the progress bar
            my $quiet = $config->{quiet};
            $config->{quiet} = 1;
            
            # try and load cache from disk if using cache
            my $tmp_cache;
            if(defined($config->{cache})) {
                $tmp_cache = load_dumped_reg_feat_cache($config, $chr, $region);
                
                #$tmp_cache = 
                #    defined($config->{tabix}) ?
                #    load_dumped_reg_feat_cache_tabix($config, $chr, $region, $trim_regions) :
                #    load_dumped_reg_feat_cache($config, $chr, $region);
                
                
                if(defined($tmp_cache->{$chr})) {
                    $count_from_cache += scalar @{$tmp_cache->{$chr}->{$_}} for keys %{$tmp_cache->{$chr}};
                }
                
                # flag as loaded
                $config->{loaded_rf}->{$chr}->{$region} = 1;
            }
            
            # no cache found on disk or not using cache
            if(!defined($tmp_cache->{$chr})) {
                
                if(defined($config->{offline})) {
                    
                    # restore quiet status
                    $config->{quiet} = $quiet;
                    
                    debug("WARNING: Could not find cache for $chr\:$region") unless defined($config->{quiet});
                    next;
                }
                
                # spoof temporary region hash
                my $tmp_hash;
                push @{$tmp_hash->{$chr}}, $region;
                
                $tmp_cache = cache_reg_feats($config, $tmp_hash);
                
                # make it an empty arrayref that gets cached
                # so we don't get confused and reload next time round
                $tmp_cache->{$chr} ||= {};
                
                $count_from_db += scalar @{$tmp_cache->{$chr}->{$_}} for keys %{$tmp_cache->{$chr}};
                
                # dump to disk if writing to cache
                #dump_reg_feat_cache($config, $tmp_cache, $chr, $region) if defined($config->{write_cache});
                (defined($config->{tabix}) ? dump_reg_feat_cache_tabix($config, $tmp_cache, $chr, $region) : dump_reg_feat_cache($config, $tmp_cache, $chr, $region)) if defined($config->{write_cache});
                
                # restore deleted coord_system adaptor
                foreach my $type(keys %{$tmp_cache->{$chr}}) {
                    $_->{slice}->{coord_system}->{adaptor} = $config->{csa} for @{$tmp_cache->{$chr}->{$type}};
                }
                
                # flag as loaded
                $config->{loaded_rf}->{$chr}->{$region} = 1;                        
            }
            
            # add loaded reg_feats to main cache
            if(defined($tmp_cache->{$chr})) {
                foreach my $type(keys %{$tmp_cache->{$chr}}) {
                    while(my $rf = shift @{$tmp_cache->{$chr}->{$type}}) {
                        
                        # filter on cell type
                        if(defined($config->{cell_type}) && scalar(@{$config->{cell_type}})) {
                            next unless grep {$rf->{cell_types}->{$_}} @{$config->{cell_type}};
                        }
                     
                        # trim out?
                        #if(defined($trim_regions) && defined($trim_regions->{$chr})) {
                        #    my $tmp_count = scalar grep {
                        #        overlap(
                        #            (split /\-/, $_)[0], (split /\-/, $_)[1],
                        #            $rf->{start}, $rf->{end}
                        #        )
                        #    } @{$trim_regions->{$chr}};
                        #    
                        #    if(!$tmp_count) {
                        #        $count_trimmed++;
                        #        next;
                        #    }
                        #}
                        
                        # track already added reg_feats by dbID
                        my $dbID = $rf->{dbID};
                        if($seen_rfs{$dbID}) {
                            $count_duplicates++;
                            next;
                        }
                        $seen_rfs{$dbID} = 1;
                        
                        push @{$rf_cache->{$chr}->{$type}}, $rf;
                    }
                }
            }
            
            undef $tmp_cache;
            
            # restore quiet status
            $config->{quiet} = $quiet;
        }
    }
    
    end_progress($config);
    
    my $rf_count = 0;
    
    foreach my $chr(keys %$rf_cache) {
        foreach my $type(keys %{$rf_cache->{$chr}}) {
            $rf_count += scalar @{$rf_cache->{$chr}->{$type}};
        }
    }
    
    debug("Retrieved $rf_count regulatory features ($count_from_mem mem, $count_from_cache cached, $count_from_db DB, $count_duplicates duplicates)") unless defined($config->{quiet});
    
    return $rf_count;
}

# gets existing VFs for a vf_hash
sub check_existing_hash {
    my $config = shift;
    my $vf_hash = shift;
    my $variation_cache;
    
    # we only care about non-SVs here
    my %new_hash;
    
    foreach my $chr(keys %{$vf_hash}) {
        foreach my $chunk(keys %{$vf_hash->{$chr}}) {
            foreach my $pos(keys %{$vf_hash->{$chr}->{$chunk}}) {
                foreach my $var(grep {$_->isa('Bio::EnsEMBL::Variation::VariationFeature')} @{$vf_hash->{$chr}->{$chunk}->{$pos}}) {
                    push @{$new_hash{$chr}->{$chunk}->{$pos}}, $var;
                }
            }
        }
    }
    
    $vf_hash = \%new_hash;
    
    debug("Checking for existing variations") unless defined($config->{quiet});
    
    my ($chunk_count, $counter);
    $chunk_count += scalar keys %{$vf_hash->{$_}} for keys %{$vf_hash};
    
    foreach my $chr(keys %{$vf_hash}) {
        
        my %loaded_regions;
        
        foreach my $chunk(keys %{$vf_hash->{$chr}}) {
            progress($config, $counter++, $chunk_count);
            
            # get the VFs for this chunk
            my ($start, $end);
            
            # work out start and end using chunk_size
            $start = $config->{chunk_size} * $chunk;
            $end = $config->{chunk_size} * ($chunk + 1);
            
            # using cache?
            if(defined($config->{cache})) {
                my $tmp_regions;
                push @{$tmp_regions->{$chr}}, $start.'-'.$end;
                
                my $converted_regions = convert_regions($config, $tmp_regions);
                
                foreach my $region(@{$converted_regions->{$chr}}) {
                
                    unless($loaded_regions{$region}) {
                        my $tmp_cache = load_dumped_variation_cache($config, $chr, $region);
                        
                        # load from DB if not found in cache
                        if(!defined($tmp_cache->{$chr})) {
                            if(defined($config->{offline})) {
                                debug("WARNING: Could not find variation cache for $chr\:$region") unless defined($config->{quiet});
                                next;
                            }
                            
                            $tmp_cache->{$chr} = get_variations_in_region($config, $chr, $region);
                            dump_variation_cache($config, $tmp_cache, $chr, $region) if defined($config->{write_cache});
                        }
                        
                        # merge tmp_cache with the main cache
                        foreach my $key(keys %{$tmp_cache->{$chr}}) {
                            $variation_cache->{$chr}->{$key} = $tmp_cache->{$chr}->{$key};
                            delete $tmp_cache->{$chr}->{$key};
                        }
                        
                        # clear memory
                        undef $tmp_cache;
                        
                        # record this region as fetched
                        $loaded_regions{$region} = 1;
                    }
                }
            }
            
            # no cache, get all variations in region from DB
            else {
                
                my ($min, $max);
                
                # we can fetch smaller region when using DB
                foreach my $pos(keys %{$vf_hash->{$chr}->{$chunk}}) {
                    foreach my $var(@{$vf_hash->{$chr}->{$chunk}->{$pos}}) {
                        foreach my $coord(qw(start end)) {
                            $min = $var->{$coord} if !defined($min) || $var->{$coord} < $min;
                            $max = $var->{$coord} if !defined($max) || $var->{$coord} > $max;
                        }
                    }
                }
                
                $variation_cache->{$chr} = get_variations_in_region($config, $chr, $min.'-'.$max);
            }
            
            # now compare retrieved vars with vf_hash
            foreach my $pos(keys %{$vf_hash->{$chr}->{$chunk}}) {
                foreach my $var(@{$vf_hash->{$chr}->{$chunk}->{$pos}}) {
                    my @found;
                    my @gmafs;
                    
                    if(defined($variation_cache->{$chr})) {
                        if(my $existing_vars = $variation_cache->{$chr}->{$pos}) {
                            foreach my $existing_var(grep {$_->[1] <= $config->{failed}} @$existing_vars) {
                                unless(is_var_novel($config, $existing_var, $var)) {
                                    push @found, $existing_var->[0];
                                    push @gmafs, $existing_var->[6].":".$existing_var->[7] if defined($existing_var->[6]) && defined($existing_var->[7]);
                                }
                            }
                        }
                    }
                    
                    $var->{existing} = join ",", @found;
                    $var->{existing} ||= '-';
                    
                    $var->{gmaf} = join ",", @gmafs;
                    $var->{gmaf} ||= undef;
                }
            }
        }
        
        delete $variation_cache->{$chr};
    }
    
    end_progress($config);
}

# gets overlapping SVs for a vf_hash
sub check_svs_hash {
    my $config = shift;
    my $vf_hash = shift;
    
    debug("Checking for overlapping structural variations") unless defined($config->{quiet});
    
    my ($chunk_count, $counter);
    $chunk_count += scalar keys %{$vf_hash->{$_}} for keys %{$vf_hash};
    
    foreach my $chr(keys %{$vf_hash}) {
        foreach my $chunk(keys %{$vf_hash->{$chr}}) {
            
            progress($config, $counter++, $chunk_count);
            
            # work out start and end using chunk_size
            my ($start, $end);
            $start = $config->{chunk_size} * $chunk;
            $end = $config->{chunk_size} * ($chunk + 1);
            
            # check for structural variations
            if(defined($config->{sa})) {
                my $slice = $config->{sa}->fetch_by_region(undef, $chr, $start, $end);
                
                if(defined($slice)) {
                    my $svs = $config->{svfa}->fetch_all_by_Slice($slice);
                    
                    foreach my $pos(keys %{$vf_hash->{$chr}->{$chunk}}) {
                        foreach my $var(@{$vf_hash->{$chr}->{$chunk}->{$pos}}) {
                            my $string = join ",",
                                map {$_->variation_name}
                                grep {$_->seq_region_start <= $var->end && $_->seq_region_end >= $var->start}
                                @$svs;
                            
                            $var->{overlapping_svs} = $string if $string;
                        }
                    }
                }
            }
        }
    }
    
    end_progress($config);
}

# gets a slice from the slice adaptor
sub get_slice {
    my $config = shift;
    my $chr = shift;
    my $otherfeatures = shift;
    $otherfeatures ||= '';
    
    my $slice;
    
    # with a FASTA DB we can just spoof slices
    if(defined($config->{fasta_db})) {
        $slice = Bio::EnsEMBL::Slice->new(
			-COORD_SYSTEM      => $config->{coord_system},
			-START             => 1,
			-END               => $config->{fasta_db}->length($chr),
			-SEQ_REGION_NAME   => $chr,
			-SEQ_REGION_LENGTH => $config->{fasta_db}->length($chr)
		);
        
        $slice->{is_fake} = 1;
        
        return $slice;
    }
    
    return undef unless defined($config->{sa}) && defined($chr);
    
    # first try to get a chromosome
    eval { $slice = $config->{$otherfeatures.'sa'}->fetch_by_region(undef, $chr); };
    
    return $slice;
}




# METHODS THAT DEAL WITH "REGIONS"
##################################

# gets regions from VF hash
sub regions_from_hash {
    my $config = shift;
    my $vf_hash = shift;
    
    my %include_regions;
    
    # if using cache we just want the regions of cache_region_size
    # since that's what we'll get from the cache (or DB if no cache found)
    if(defined($config->{cache})) {
        
        my $region_size = $config->{cache_region_size};
        
        foreach my $chr(keys %$vf_hash) {
            $include_regions{$chr} = [];
            my %temp_regions;
            
            foreach my $chunk(keys %{$vf_hash->{$chr}}) {
                foreach my $pos(keys %{$vf_hash->{$chr}{$chunk}}) {
                    my ($s, $e) = ($pos - MAX_DISTANCE_FROM_TRANSCRIPT, $pos + MAX_DISTANCE_FROM_TRANSCRIPT);
                    
                    my $low = int ($s / $region_size);
                    my $high = int ($e / $region_size) + 1;
                    
                    for my $i($low..($high - 1)) {
                        $temp_regions{(($i * $region_size) + 1).'-'.(($i + 1) * $region_size)} = 1;
                    }
                }
            }
            
            @{$include_regions{$chr}} = keys %temp_regions;
        }
    }
    
    # if no cache we don't want to fetch more than is necessary, so find the
    # minimum covered region of the variations in the hash
    else {
        foreach my $chr(keys %$vf_hash) {
            $include_regions{$chr} = [];
            
            foreach my $chunk(keys %{$vf_hash->{$chr}}) {
                foreach my $pos(keys %{$vf_hash->{$chr}{$chunk}}) {
                    add_region($_->start, $_->end, $include_regions{$chr}) for @{$vf_hash->{$chr}{$chunk}{$pos}};
                }
            }
        }
        
        # merge regions
        merge_regions(\%include_regions, $config);
    }
    
    return \%include_regions;
}

# adds a region to region list, expanding existing one if overlaps
sub add_region {
    my $start = shift;
    my $end = shift;
    my $region_list = shift;
    
    # fix end for insertions
    $end = $start if $end < $start;
    
    my $added = 0;
    my $i = 0;
    
    while ($i < scalar @$region_list) {
        my ($region_start, $region_end) = split /\-/, $region_list->[$i];
        
        if($start <= $region_end && $end >= $region_start) {
            my $new_region_start = ($start < $end ? $start : $end) - MAX_DISTANCE_FROM_TRANSCRIPT;
            my $new_region_end = ($start > $end ? $start : $end) + MAX_DISTANCE_FROM_TRANSCRIPT;
            
            $new_region_start = 1 if $new_region_start < 1;
            
            $region_start = $new_region_start if $new_region_start < $region_start;
            $region_end = $new_region_end if $new_region_end > $region_end;
            
            $region_list->[$i] = $region_start.'-'.$region_end;
            $added = 1;
        }
        
        $i++;
    }
    
    unless($added) {
        my $s = $start - MAX_DISTANCE_FROM_TRANSCRIPT;
        $s = 1 if $s < 1;
        
        push @{$region_list}, $s.'-'.($end + MAX_DISTANCE_FROM_TRANSCRIPT);
    }
}

# merges overlapping regions from scans
sub merge_regions {
    my $include_regions = shift;
    my $config = shift;
    my $consecutive = shift;
    $consecutive ||= 0;
    
    # now merge overlapping regions
    foreach my $chr(keys %$include_regions) {
        my $max_index = $#{$include_regions->{$chr}};
        my (@new_regions, %skip);
        
        for my $i(0..$max_index) {
            next if $skip{$i};
            my ($s, $e) = split /\-/, $include_regions->{$chr}[$i];
            
            for my $j(($i+1)..$max_index) {
                next if $skip{$j};
                my ($ns, $ne) = split /\-/, $include_regions->{$chr}[$j];
                
                if($s <= ($ne + $consecutive) && $e >= ($ns - $consecutive)) {
                    $s = $ns if $ns < $s;
                    $e = $ne if $ne > $e;
                    
                    $skip{$j} = 1;
                }
            }
            
            push @new_regions, $s.'-'.$e;
        }
        
        # replace original
        $include_regions->{$chr} = \@new_regions;
        
        $config->{region_count} += scalar @new_regions;
    }
    
    return $include_regions;
}

# converts regions as determined by scan_file to regions loadable from cache
sub convert_regions {
    my $config = shift;
    my $regions = shift;
    
    return undef unless defined $regions;
    
    my $region_size = $config->{cache_region_size};
    
    my %new_regions;
    
    foreach my $chr(keys %$regions) {
        my %temp_regions;
        
        foreach my $region(@{$regions->{$chr}}) {
            my ($s, $e) = split /\-/, $region;
            
            my $low = int ($s / $region_size);
            my $high = int ($e / $region_size) + 1;
            
            for my $i($low..($high - 1)) {
                $temp_regions{(($i * $region_size) + 1).'-'.(($i + 1) * $region_size)} = 1;
            }
        }
        
        @{$new_regions{$chr}} = keys %temp_regions;
    }
    
    return \%new_regions;
}





# CACHE METHODS
###############

# prunes a cache to get rid of features not in regions in use
sub prune_cache {
    my $config  = shift;
    my $cache   = shift;
    my $regions = shift;
    my $loaded  = shift;
    
    # delete no longer in use chroms
    foreach my $chr(keys %$cache) {
        delete $cache->{$chr} unless defined $regions->{$chr} && scalar @{$regions->{$chr}};
    }
    
    my $new_count = 0;
    
    foreach my $chr(keys %$cache) {
        
        # get total area spanned by regions    
        my ($min, $max);
        foreach my $region(@{$regions->{$chr}}) {
            my ($s, $e) = split /\-/, $region;
            $min = $s if !defined($min) or $s < $min;
            $max = $e if !defined($max) or $e > $max;
        }
        
        # transcript cache
        if(ref($cache->{$chr}) eq 'ARRAY') {
            $cache->{$chr} = prune_min_max($cache->{$chr}, $min, $max);
            $new_count += scalar @{$cache->{$chr}};
        }
        # regfeat cache
        elsif(ref($cache->{$chr}) eq 'HASH') {
            for(keys %{$cache->{$chr}}) {
                $cache->{$chr}->{$_} = prune_min_max($cache->{$chr}->{$_}, $min, $max);
                $new_count += scalar @{$cache->{$chr}->{$_}};
            }
        }
        
        # update loaded regions
        my %have_regions = map {$_ => 1} @{$regions->{$chr}};
        
        foreach my $region(keys %{$loaded->{$chr}}) {
            delete $loaded->{$chr}->{$region} unless defined $have_regions{$region};
        }
    }
    
    return $new_count;
}

# does the actual pruning
sub prune_min_max {
    my $array = shift;
    my $min   = shift;
    my $max   = shift;
    
    # splice out features not in area spanned by min/max
    my $i = 0;
    my $f_count = scalar @{$array};
    my @new_cache;
    
    while($i < $f_count) {
        my $f = $array->[$i];
        
        $i++;
        
        if($max - $f->start() > 0 && $f->end - $min > 0) {
            push @new_cache, $f;
        }
        
        # do some cleaning for transcripts
        elsif(defined $f->{translation}) {
            delete $f->{translation}->{transcript};
            delete $f->{translation};
        }
    }
    
    undef $array;
    return \@new_cache;
}

# get transcripts for slices
sub cache_transcripts {
    my $config = shift;
    my $include_regions = shift;
    
    my $tr_cache;
    my $i;
    
    debug("Caching transcripts") unless defined($config->{quiet});
    
    foreach my $chr(keys %$include_regions) {
        
        my $slice = get_slice($config, $chr);
        
        next unless defined $slice;
        
        # prefetch some things
        $slice->is_circular;
        
        # trim bumf off the slice
        delete $slice->{coord_system}->{adaptor} if defined($config->{write_cache});
        
        # no regions?
        if(!scalar @{$include_regions->{$chr}}) {
            my $start = 1;
            my $end = $config->{cache_region_size};
            
            while($start < $slice->end) {
                push @{$include_regions->{$chr}}, $start.'-'.$end;
                $start += $config->{cache_region_size};
                $end += $config->{cache_region_size};
            }
        }
        
        my $region_count;
        
        if(scalar keys %$include_regions == 1) {
            my ($chr) = keys %$include_regions;
            $region_count = scalar @{$include_regions->{$chr}};
            debug("Caching transcripts for chromosome $chr") unless defined($config->{quiet});
        }
        
        foreach my $region(@{$include_regions->{$chr}}) {
            progress($config, $i++, $region_count || $config->{region_count});
            
            my ($s, $e) = split /\-/, $region;
            
            # sanity check start and end
            $s = 1 if $s < 1;
            $e = $slice->end if $e > $slice->end;
            
            # get sub-slice
            my $sub_slice = $slice->sub_Slice($s, $e);
            
            # for some reason unless seq is called here the sequence becomes Ns later
            $sub_slice->seq;
            
            # add transcripts to the cache, via a transfer to the chrom's slice
            if(defined($sub_slice)) {
                foreach my $gene(map {$_->transfer($slice)} @{$sub_slice->get_all_Genes(undef, undef, 1)}) {
                    my $gene_stable_id = $gene->stable_id;
                    my $canonical_tr_id = $gene->{canonical_transcript_id};
                    
                    my @trs;
                    
                    foreach my $tr(@{$gene->get_all_Transcripts}) {
                        $tr->{_gene_stable_id} = $gene_stable_id;
                        $tr->{_gene} = $gene;
                        
                        # indicate if canonical
                        $tr->{is_canonical} = 1 if defined $canonical_tr_id and $tr->dbID eq $canonical_tr_id;
                        
                        if(defined($config->{prefetch})) {
                            prefetch_transcript_data($config, $tr);
                        }
                        
                        # CCDS
                        elsif(defined($config->{ccds})) {
                            my @entries = grep {$_->database eq 'CCDS'} @{$tr->get_all_DBEntries};
                            $tr->{_ccds} = $entries[0]->display_id if scalar @entries;
                        }
                        
                        # strip some unnecessary data from the transcript object
                        clean_transcript($tr) if defined($config->{write_cache});
                        
                        push @trs, $tr;
                    }
                    
                    # sort the transcripts by translation so we can share sift/polyphen stuff
                    # between transcripts and save cache space
                    if(defined($config->{write_cache}) && (defined($config->{sift}) || defined($config->{polyphen}))) {
                        
                        my $prev_tr;
                        
                        # sort them by peptide seqeuence as transcripts with identical peptides
                        # will have identical SIFT/PolyPhen prediction strings
                        foreach my $tr(sort {$a->{_variation_effect_feature_cache}->{peptide} cmp $b->{_variation_effect_feature_cache}->{peptide}} grep {$_->{_variation_effect_feature_cache}->{peptide}} @trs)  {
                            
                            if(
                                defined($prev_tr) &&
                                $prev_tr->{_variation_effect_feature_cache}->{peptide}
                                    eq $tr->{_variation_effect_feature_cache}->{peptide}
                            ) {
                                
                                foreach my $analysis(qw(sift polyphen)) {
                                    next unless defined($config->{$analysis});
                                    $tr->{_variation_effect_feature_cache}->{protein_function_predictions}->{$analysis} = $prev_tr->{_variation_effect_feature_cache}->{protein_function_predictions}->{$analysis};
                                }
                            }
                            
                            $prev_tr = $tr;
                        }
                    }
                    
                    # clean the gene
                    clean_gene($gene);
                    
                    push @{$tr_cache->{$chr}}, @trs;
                }
            }
        }
    }
    
    end_progress($config);
    
    return $tr_cache;
}

# gets rid of extra bits of info attached to the transcript that we don't need
sub clean_transcript {
    my $tr = shift;
    
    foreach my $key(qw(display_xref external_db external_display_name external_name external_status created_date status description edits_enabled modified_date dbentries is_current analysis transcript_mapper)) {
        delete $tr->{$key} if defined($tr->{$key});
    }
    
    # clean all attributes but miRNA
    if(defined($tr->{attributes})) {
        my @new_atts;
        foreach my $att(@{$tr->{attributes}}) {
            push @new_atts, $att if $att->{code} eq 'miRNA';
        }
        $tr->{attributes} = \@new_atts;
    }
    
    # clean the translation
    if(defined($tr->translation)) {
        
        # sometimes the translation points to a different transcript?
        $tr->{translation}->{transcript} = $tr;
        weaken($tr->{translation}->{transcript});
        
        for my $key(qw(attributes protein_features created_date modified_date)) {
            delete $tr->translation->{$key};
        }
    }
}

# gets rid of extra bits of info attached to genes. At the moment this is almost
# everything as genes are only used for their locations when looking at
# structural variations
sub clean_gene {
    my $gene = shift;
    
    # delete almost everything in the gene
    map {delete $gene->{$_}}
        grep {
            $_ ne 'start' &&
            $_ ne 'end' &&
            $_ ne 'strand' && 
            $_ ne 'stable_id'
        }
    keys %{$gene};
}

# build slice cache from transcript cache
sub build_slice_cache {
    my $config = shift;
    my $tr_cache = shift;
    
    my %slice_cache;
    
    foreach my $chr(keys %$tr_cache) {
        $slice_cache{$chr} = scalar @{$tr_cache->{$chr}} ? $tr_cache->{$chr}[0]->slice : &get_slice($config, $chr);
        
        if(!defined($slice_cache{$chr})) {
            delete $slice_cache{$chr}
        }
        
        else {
            # reattach adaptor to the coord system
            $slice_cache{$chr}->{coord_system}->{adaptor} ||= $config->{csa};
        }
    }
    
    return \%slice_cache;
}

# pre-fetches per-transcript data
sub prefetch_transcript_data {
    my $config = shift;
    my $tr = shift;
    
    # introns
    my $introns = $tr->get_all_Introns;
    
    if(defined($introns)) {
        foreach my $intron(@$introns) {
            foreach my $key(qw(adaptor analysis dbID next prev seqname)) {
                delete $intron->{$key};
            }
        }
    }
    
    $tr->{_variation_effect_feature_cache}->{introns} ||= $introns;
    
    # translateable_seq, mapper
    $tr->{_variation_effect_feature_cache}->{translateable_seq} ||= $tr->translateable_seq;
    $tr->{_variation_effect_feature_cache}->{mapper} ||= $tr->get_TranscriptMapper;
    
    # peptide
    unless ($tr->{_variation_effect_feature_cache}->{peptide}) {
        my $translation = $tr->translate;
        $tr->{_variation_effect_feature_cache}->{peptide} = $translation ? $translation->seq : undef;
    }
    
    # protein features
    if(defined($config->{domains}) || defined($config->{write_cache})) {
        my $pfs = $tr->translation ? $tr->translation->get_all_ProteinFeatures : [];
        
        # clean them to save cache space
        foreach my $pf(@$pfs) {
            
            # remove everything but the coord, analysis and ID fields
            foreach my $key(keys %$pf) {
                delete $pf->{$key} unless
                    $key eq 'start' ||
                    $key eq 'end' ||
                    $key eq 'analysis' ||
                    $key eq 'hseqname';
            }
            
            # remove everything from the analysis but the display label
            foreach my $key(keys %{$pf->{analysis}}) {
                delete $pf->{analysis}->{$key} unless $key eq '_display_label';
            }
        }
        
        $tr->{_variation_effect_feature_cache}->{protein_features} = $pfs;
    }
    
    # codon table
    unless ($tr->{_variation_effect_feature_cache}->{codon_table}) {
        # for mithocondrial dna we need to to use a different codon table
        my $attrib = $tr->slice->get_all_Attributes('codon_table')->[0];
        
        $tr->{_variation_effect_feature_cache}->{codon_table} = $attrib ? $attrib->value : 1;
    }
    
    # sift/polyphen
    if(defined($config->{pfpma}) && defined($tr->{_variation_effect_feature_cache}->{peptide})) {
        foreach my $analysis(qw(sift polyphen)) {
            next unless defined($config->{$analysis});
            my $a = $analysis;
            $a .= '_humvar' if $a eq 'polyphen';
            $tr->{_variation_effect_feature_cache}->{protein_function_predictions}->{$a} ||= $config->{pfpma}->fetch_by_analysis_translation_md5($a, md5_hex($tr->{_variation_effect_feature_cache}->{peptide}))
        }
    }
    
    # gene
    $tr->{_gene} ||= $config->{ga}->fetch_by_transcript_stable_id($tr->stable_id);
    
    # gene HGNC
    if(defined $config->{hgnc}) {
        
        # get from gene cache if found already
        if(defined($tr->{_gene}->{_hgnc})) {
            $tr->{_gene_hgnc} = $tr->{_gene}->{_hgnc};
        }
        else {
            my @entries = grep {$_->database eq 'HGNC'} @{$tr->{_gene}->get_all_DBEntries()};
            if(scalar @entries) {
                $tr->{_gene_hgnc} = $entries[0]->display_id;
            }
            
            $tr->{_gene_hgnc} ||= '-';
            
            # cache it on the gene object too
            $tr->{_gene}->{_hgnc} = $tr->{_gene_hgnc};
        }
    }
    
    # CCDS
    my @entries = grep {$_->database eq 'CCDS'} @{$tr->get_all_DBEntries};
    $tr->{_ccds} = $entries[0]->display_id if scalar @entries;
    $tr->{_ccds} ||= '-';
    
    # refseq
    @entries = grep {$_->database eq 'RefSeq_mRNA'} @{$tr->get_all_DBEntries};
    if(scalar @entries) {
        $tr->{_refseq} = join ",", map {$_->display_id} @entries;
    }
    else {
        $tr->{_refseq} = '-';
    }
    
    # protein stable ID
    $tr->{_protein} = $tr->translation ? $tr->translation->stable_id : '-';
    
    return $tr;
}

sub get_dump_file_name {
    my $config = shift;
    my $chr    = shift;
    my $region = shift;
    my $type   = shift;
    
    $type ||= 'transcript';
    
    if($type eq 'transcript') {
        $type = '';
    }
    else {
        $type = '_'.$type;
    }
    
    #my ($s, $e) = split /\-/, $region;
    #my $subdir = int($s / 1e6);
    #
    #my $dir = $config->{dir}.'/'.$chr.'/'.$subdir;
    
    my $dir = $config->{dir}.'/'.$chr;
    my $dump_file = $dir.'/'.$region.$type.(defined($config->{tabix}) ? '_tabix' : '').'.gz';
    
    # make directory if it doesn't exist
    if(!(-e $dir) && defined($config->{write_cache})) {
        make_path($dir);
    }
    
    return $dump_file;
}

# dumps out transcript cache to file
sub dump_transcript_cache {
    my $config = shift;
    my $tr_cache = shift;
    my $chr = shift;
    my $region = shift;
    
    debug("Dumping cached transcript data") unless defined($config->{quiet});
    
    # clean the slice adaptor before storing
    clean_slice_adaptor($config);
    
    strip_transcript_cache($config, $tr_cache);
    
    $config->{reg}->disconnect_all;
    delete $config->{sa}->{dbc}->{_sql_helper};
    
    my $dump_file = get_dump_file_name($config, $chr, $region, 'transcript');
    
    debug("Writing to $dump_file") unless defined($config->{quiet});
    
    # storable
    open my $fh, "| gzip -9 -c > ".$dump_file or die "ERROR: Could not write to dump file $dump_file";
    nstore_fd($tr_cache, $fh);
    close $fh;
}

#sub dump_transcript_cache_tabix {
#    my $config = shift;
#    my $tr_cache = shift;
#    my $chr = shift;
#    my $region = shift;
#    
#    debug("Dumping cached transcript data") unless defined($config->{quiet});
#    
#    # clean the slice adaptor before storing
#    clean_slice_adaptor($config);
#    
#    strip_transcript_cache($config, $tr_cache);
#    
#    $config->{reg}->disconnect_all;
#    
#    my $dir = $config->{dir}.'/'.$chr;
#    my $dump_file = $dir.'/'.($region || "dump").'_tabix.gz';
#    
#    # make directory if it doesn't exist
#    if(!(-e $dir)) {
#        make_path($dir);
#    }
#    
#    debug("Writing to $dump_file") unless defined($config->{quiet});
#    
#    use Storable qw(nfreeze);
#    use MIME::Base64 qw(encode_base64);
#    #open NEW, "| bgzip -c > ".$dump_file or die "ERROR: Could not write to dump file $dump_file";
#    #
#    #foreach my $tr(sort {$a->start <=> $b->start} @{$tr_cache->{$chr}}) {
#    #    print NEW join "\t", (
#    #        $chr,
#    #        $tr->start,
#    #        $tr->end,
#    #        encode_base64(freeze($tr), "")
#    #    );
#    #    print NEW "\n";
#    #}
#    #close NEW;
#    #
#    ## tabix it
#    #my $output = `tabix -s 1 -b 2 -e 3 -f $dump_file 2>&1`;
#    #die("ERROR: Failed during tabix indexing\n$output\n") if $output;
#    open NEW, "| gzip -9 -c > ".$dump_file or die "ERROR: Could not write to dump file $dump_file";
#    
#    foreach my $tr(sort {$a->start <=> $b->start} @{$tr_cache->{$chr}}) {
#        print NEW join "\t", (
#            $chr,
#            $tr->start,
#            $tr->end,
#            encode_base64(freeze($tr), "")
#        );
#        print NEW "\n";
#    }
#    close NEW;
#}

# loads in dumped transcript cache to memory
sub load_dumped_transcript_cache {
    my $config = shift;
    my $chr = shift;
    my $region = shift;
    
    my $dump_file = get_dump_file_name($config, $chr, $region, 'transcript');
    
    return undef unless -e $dump_file;
    
    debug("Reading cached transcript data for chromosome $chr".(defined $region ? "\:$region" : "")." from dumped file") unless defined($config->{quiet});
    
    open my $fh, $config->{compress}." ".$dump_file." |" or return undef;
    my $tr_cache;
    $tr_cache = fd_retrieve($fh);
    close $fh;
    
    # reattach adaptors
    foreach my $t(@{$tr_cache->{$chr}}) {
        if(defined($t->{translation})) {
            $t->{translation}->{adaptor} = $config->{tra} if defined $t->{translation}->{adaptor};
            $t->{translation}->{transcript} = $t;
            weaken($t->{translation}->{transcript});
        }
        
        $t->{slice}->{adaptor} = $config->{sa};
    }
    
    return $tr_cache;
}

#sub load_dumped_transcript_cache_tabix {
#    my $config = shift;
#    my $chr = shift;
#    my $region = shift;
#    my $trim_regions = shift;
#    
#    my $dir = $config->{dir}.'/'.$chr;
#    my $dump_file = $dir.'/'.($region || "dump").'_tabix.gz';
#    
#    #print STDERR "Reading from $dump_file\n";
#    
#    return undef unless -e $dump_file;
#    
#    debug("Reading cached transcript data for chromosome $chr".(defined $region ? "\:$region" : "")." from dumped file") unless defined($config->{quiet});
#    
#    my $tr_cache;
#    
#    use MIME::Base64 qw(decode_base64);
#    use Storable qw(thaw);
#    
#    my ($s, $e) = split /\-/, $region;
#    my @regions = grep {overlap($s, $e, (split /\-/, $_))} @{$trim_regions->{$chr}};
#    my $regions = "";
#    $regions .= " $chr\:$_" for @regions;
#    
#    #print STDERR "tabix $dump_file $regions |\n";
#    #open IN, "tabix $dump_file $regions |";
#    open IN, "gzip -dc $dump_file |";
#    while(<IN>) {
#    
#    #$DB::single = 1;
#        my ($chr, $start, $end, $blob) = split /\t/, $_;
#        next unless grep {overlap($start, $end, (split /\-/, $_))} @regions;
#        my $tr = thaw(decode_base64($blob));
#        push @{$tr_cache->{$chr}}, $tr;
#    }
#    close IN;
#    
#    # reattach adaptors
#    foreach my $t(@{$tr_cache->{$chr}}) {
#        if(defined($t->{translation})) {
#            $t->{translation}->{adaptor} = $config->{tra} if defined $t->{translation}->{adaptor};
#            $t->{translation}->{transcript} = $t;
#            weaken($t->{translation}->{transcript});
#        }
#        
#        $t->{slice}->{adaptor} = $config->{sa};
#    }
#    
#    # add empty array ref so code doesn't try and fetch from DB too
#    $tr_cache->{$chr} ||= [];
#    
#    return $tr_cache;
#}

# strips cache before writing to disk
sub strip_transcript_cache {
    my $config = shift;
    my $cache = shift;
    
    foreach my $chr(keys %$cache) {
        foreach my $tr(@{$cache->{$chr}}) {
            foreach my $exon(@{$tr->{_trans_exon_array}}) {
                delete $exon->{slice}->{adaptor};
                
                for(qw(adaptor created_date modified_date is_current version is_constitutive _seq_cache dbID slice)) {
                    delete $exon->{$_};
                }
            }
            
            delete $tr->{adaptor};
            delete $tr->{slice}->{adaptor};
        }
    }
}

# cleans slice adaptor before storing in cache
sub clean_slice_adaptor{
    my $config = shift;
    
    # clean some stuff off the slice adaptor
    delete $config->{sa}->{asm_exc_cache};
    $config->{sa}->{sr_name_cache} = {};
    $config->{sa}->{sr_id_cache} = {};
    delete $config->{sa}->{db}->{seq_region_cache};
    delete $config->{sa}->{db}->{name_cache};
}


# dump adaptors to cache
sub dump_adaptor_cache {
    my $config = shift;
    
    $config->{reg}->disconnect_all;
    delete $config->{sa}->{dbc}->{_sql_helper};
    
    my $dir = $config->{dir};
    my $dump_file = $dir.'/adaptors.gz';
    
    # make directory if it doesn't exist
    if(!(-e $dir)) {
        make_path($dir);
	}
    
    open my $fh, "| gzip -9 -c > ".$dump_file or die "ERROR: Could not write to dump file $dump_file";
    nstore_fd($config, $fh);
    close $fh;
}

# load dumped adaptors
sub load_dumped_adaptor_cache {
    my $config = shift;
    
    my $dir = $config->{dir};
    my $dump_file = $dir.'/adaptors.gz';
    
    return undef unless -e $dump_file;
    
    debug("Reading cached adaptor data") unless defined($config->{quiet});
    
    open my $fh, $config->{compress}." ".$dump_file." |" or return undef;
    my $cached_config;
    $cached_config = fd_retrieve($fh);
    close $fh;
    
    $config->{$_} = $cached_config->{$_} for qw(sa ga ta vfa svfa tva pfpma mca csa RegulatoryFeature_adaptor MotifFeature_adaptor);
    
    return 1;
}

# dumps cached variations to disk
sub dump_variation_cache {
    my $config = shift;
    my $v_cache = shift;
    my $chr = shift;
    my $region = shift;
    
    my $dump_file = get_dump_file_name($config, $chr, $region, 'var');
    
    open DUMP, "| gzip -9 -c > ".$dump_file or die "ERROR: Could not write to adaptor dump file $dump_file";
    
    foreach my $pos(keys %{$v_cache->{$chr}}) {
        foreach my $v(@{$v_cache->{$chr}->{$pos}}) {
            my ($name, $failed, $start, $end, $as, $strand, $ma, $maf) = @$v;
            
            print DUMP join " ", (
                $name,
                $failed == 0 ? '' : $failed,
                $start,
                $end == $start ? '' : $end,
                $as,
                $strand == 1 ? '' : $strand,
                $ma || '',
                defined($maf) ? sprintf("%.2f", $maf) : '',
            );
            print DUMP "\n";
        }
    }
    
    close DUMP;    
}

# loads dumped variation cache
sub load_dumped_variation_cache {
    my $config = shift;
    my $chr = shift;
    my $region = shift;
    
    my $dump_file = get_dump_file_name($config, $chr, $region, 'var');
    
    return undef unless -e $dump_file;
    
    open DUMP, $config->{compress}." ".$dump_file." |" or return undef;
    
    my $v_cache;
    
    while(<DUMP>) {
        chomp;
        my ($name, $failed, $start, $end, $as, $strand, $ma, $maf) = split / /, $_;
        $failed ||= 0;
        $end ||= $start;
        $strand ||= 1;
        $ma ||= undef;
        $maf ||= undef;
        
        my @v = ($name, $failed, $start, $end, $as, $strand, $ma, $maf);
        push @{$v_cache->{$chr}->{$start}}, \@v;
    }
    
    close DUMP;
    
    return $v_cache;
}

# caches regulatory features
sub cache_reg_feats {
    my $config = shift;
    my $include_regions = shift;
    
    my $rf_cache;
    my $i;
    
    debug("Caching regulatory features") unless defined($config->{quiet});
    
    foreach my $chr(keys %$include_regions) {
        
        my $slice = get_slice($config, $chr);
        
        next unless defined $slice;
        
        # prefetch some things
        $slice->is_circular;
        
        # no regions?
        if(!scalar @{$include_regions->{$chr}}) {
            my $start = 1;
            my $end = $config->{cache_region_size};
            
            while($start < $slice->end) {
                push @{$include_regions->{$chr}}, $start.'-'.$end;
                $start += $config->{cache_region_size};
                $end += $config->{cache_region_size};
            }
        }
        
        my $region_count;
        
        if(scalar keys %$include_regions == 1) {
            my ($chr) = keys %$include_regions;
            $region_count = scalar @{$include_regions->{$chr}};
            debug("Caching transcripts for chromosome $chr") unless defined($config->{quiet});
        }
        
        foreach my $region(@{$include_regions->{$chr}}) {
            progress($config, $i++, $region_count || $config->{region_count});
            
            my ($s, $e) = split /\-/, $region;
            
            # sanity check start and end
            $s = 1 if $s < 1;
            $e = $slice->end if $e > $slice->end;
            
            # get sub-slice
            my $sub_slice = $slice->sub_Slice($s, $e);
            $sub_slice->{coord_system}->{adaptor} = $config->{csa};
            
            next unless defined($sub_slice);
            
            foreach my $type(@REG_FEAT_TYPES) {
                my $features = $config->{$type.'_adaptor'}->fetch_all_by_Slice($sub_slice);
                next unless defined($features);
                
                # cell types
                if(defined($config->{cell_type}) && scalar(@{$config->{cell_type}})) {
                    foreach my $rf(@$features) {
                        
                        my %cl;
                        
                        # get cell type by fetching all from stable ID
                        if($type eq 'RegulatoryFeature') {
                            %cl = map {
                                $_->feature_set->cell_type->name => $_->feature_type->name
                            } @{$rf->adaptor->fetch_all_by_stable_ID($rf->stable_id)};
                        }
                        
                        # get cell type by fetching regfeats that contain this MotifFeature
                        elsif($type eq 'MotifFeature') {
                            %cl = map {
                                $_->feature_set->cell_type->name => $_->feature_type->name
                            } @{$config->{'RegulatoryFeature_adaptor'}->fetch_all_by_attribute_feature($rf)};
                        }
                        
                        $rf->{cell_types} = \%cl;
                    }
                }
                
                push @{$rf_cache->{$chr}->{$type}},
                    map { clean_reg_feat($_) }
                    map { $_->transfer($slice) }
                    @{$features};
            }
        }
    }
    
    end_progress($config);
    
    return $rf_cache;
}


# cleans reg feats for caching
sub clean_reg_feat {
    my $rf = shift;
    
    foreach my $key(qw/adaptor binary_string bound_start bound_end attribute_cache feature_type feature_set analysis/) {
        delete $rf->{$key};
    }
    
    if(defined($rf->{binding_matrix})) {
        $rf->{_variation_effect_feature_cache}->{seq} = $rf->seq;
        
        foreach my $key(qw/adaptor feature_type analysis dbID/) {
            delete $rf->{binding_matrix}->{$key};
        }
    }
    
    return $rf;
}


# dumps out reg feat cache to file
sub dump_reg_feat_cache {
    my $config = shift;
    my $rf_cache = shift;
    my $chr = shift;
    my $region = shift;
    
    debug("Dumping cached reg feat data for $chr:$region") unless defined($config->{quiet});
    
    # clean the slice adaptor before storing
    clean_slice_adaptor($config);
    
    $config->{reg}->disconnect_all;
    delete $config->{sa}->{dbc}->{_sql_helper};
    
    foreach my $chr(keys %{$rf_cache}) {
        foreach my $type(keys %{$rf_cache->{$chr}}) {
            delete $_->{slice}->{coord_system}->{adaptor} for @{$rf_cache->{$chr}->{$type}};
        }
    }
    
    my $dump_file = get_dump_file_name($config, $chr, $region, 'reg');
    
    debug("Writing to $dump_file") unless defined($config->{quiet});
    
    # storable
    open my $fh, "| gzip -9 -c > ".$dump_file or die "ERROR: Could not write to dump file $dump_file";
    nstore_fd($rf_cache, $fh);
    close $fh;
}

#sub dump_reg_feat_cache_tabix {
#    my $config = shift;
#    my $rf_cache = shift;
#    my $chr = shift;
#    my $region = shift;
#    
#    debug("Dumping cached reg feat data") unless defined($config->{quiet});
#    
#    # clean the slice adaptor before storing
#    clean_slice_adaptor($config);
#    
#    $config->{reg}->disconnect_all;
#    delete $config->{sa}->{dbc}->{_sql_helper};
#    
#    $config->{reg}->disconnect_all;
#    
#    my $dump_file = get_dump_file_name($config, $chr, $region, 'reg');
#    
#    debug("Writing to $dump_file") unless defined($config->{quiet});
#    
#    use Storable qw(nfreeze);
#    use MIME::Base64 qw(encode_base64);
#    open NEW, "| gzip -9 -c > ".$dump_file or die "ERROR: Could not write to dump file $dump_file";
#    
#    foreach my $type(keys %{$rf_cache->{$chr}}) {
#        foreach my $rf(sort {$a->start <=> $b->start} @{$rf_cache->{$chr}->{$type}}) {
#            print NEW join "\t", (
#                $chr,
#                $rf->start,
#                $rf->end,
#                $type,
#                encode_base64(freeze($rf), "")
#            );
#            print NEW "\n";
#        }
#    }
#    close NEW;
#}

# loads in dumped transcript cache to memory
sub load_dumped_reg_feat_cache {
    my $config = shift;
    my $chr = shift;
    my $region = shift;
    
    my $dump_file = get_dump_file_name($config, $chr, $region, 'reg');
    
    return undef unless -e $dump_file;
    
    debug("Reading cached reg feat data for chromosome $chr".(defined $region ? "\:$region" : "")." from dumped file") unless defined($config->{quiet});
    
    open my $fh, $config->{compress}." ".$dump_file." |" or return undef;
    my $rf_cache;
    $rf_cache = fd_retrieve($fh);
    close $fh;
    
    return $rf_cache;
}



#sub load_dumped_reg_feat_cache_tabix {
#    my $config = shift;
#    my $chr = shift;
#    my $region = shift;
#    my $trim_regions = shift;
#    
#    my $dump_file = get_dump_file_name($config, $chr, $region, 'reg');
#    
#    #print STDERR "Reading from $dump_file\n";
#    
#    return undef unless -e $dump_file;
#    
#    debug("Reading cached reg feat data for chromosome $chr".(defined $region ? "\:$region" : "")." from dumped file") unless defined($config->{quiet});
#    
#    my $rf_cache;
#    
#    use MIME::Base64 qw(decode_base64);
#    use Storable qw(thaw);
#    
#    my ($s, $e) = split /\-/, $region;
#    my @regions = grep {overlap($s, $e, (split /\-/, $_))} @{$trim_regions->{$chr}};
#    my $regions = "";
#    $regions .= " $chr\:$_" for @regions;
#    
#    #print STDERR "tabix $dump_file $regions |\n";
#    #open IN, "tabix $dump_file $regions |";
#    open IN, "gzip -dc $dump_file |";
#    while(<IN>) {
#        my ($chr, $start, $end, $type, $blob) = split /\t/, $_;
#        next unless grep {overlap($start, $end, (split /\-/, $_))} @regions;
#        my $rf = thaw(decode_base64($blob));
#        push @{$rf_cache->{$chr}->{$type}}, $rf;
#    }
#    close IN;
#    
#    $rf_cache->{$chr}->{$_} ||= [] for @REG_FEAT_TYPES;
#    
#    return $rf_cache;
#}


# get custom annotation for a region
sub cache_custom_annotation {
    my $config = shift;
    my $include_regions = shift;
    my $chr = shift;
    
    #$include_regions = merge_regions($include_regions, $config, 1);
    
    my $annotation = {};
    
    my $total = scalar @{$config->{custom}} * scalar @{$include_regions->{$chr}}; 
    my $counter = 0;
    
    my $max_regions_per_tabix = 1000;
    
    debug("Caching custom annotations") unless defined($config->{quiet});
    
    foreach my $custom(@{$config->{custom}}) {
      
        my @regions = @{$include_regions->{$chr}};
        
        while(scalar @regions) {
            my $got_features = 0;
            
            my @tmp_regions = splice @regions, 0, $max_regions_per_tabix;
            
            progress($config, $counter, $total);
            $counter += scalar @tmp_regions;
            
            # some files may have e.g. chr10 instead of 10
            for my $tmp_chr($chr, 'chr'.$chr) {
                
                # bigwig needs to use bigWigToWig utility
                if($custom->{format} eq 'bigwig') {
                    foreach my $region(@tmp_regions) {
                        my ($s, $e) = split /\-/, $region;
                        my $tmp_file = $config->{tmpdir}.'/vep_tmp_'.$$.'_'.$tmp_chr.'_'.$s.'_'.$e;
                        my $bigwig_file = $custom->{file};
                        my $bigwig_output = `bigWigToWig -chrom=$tmp_chr -start=$s -end=$e $bigwig_file $tmp_file 2>&1`;
                        
                        die "\nERROR: Problem using bigwig file $bigwig_file\n$bigwig_output" if $bigwig_output;
                    }
                    
                    # concatenate all the files together
                    my $string = $config->{tmpdir}.'/vep_tmp_'.$$.'_*';
                    my $tmp_file = $config->{tmpdir}.'/vep_tmp_'.$$;
                    `cat $string > $tmp_file; rm $string`;
                    open CUSTOM, $tmp_file
                        or die "\nERROR: Could not read from temporary WIG file $tmp_file\n";
                }
                
                # otherwise use tabix
                else {
                    # tabix can fetch multiple regions, so construct a string
                    my $region_string = join " ", map {$tmp_chr.':'.$_} @tmp_regions;
                    
                    open CUSTOM, "tabix ".$custom->{file}." $region_string 2>&1 |"
                        or die "\nERROR: Could not open tabix pipe for ".$custom->{file}."\n";
                }
                
                # set an error flag so we don't have to check every line
                my $error_flag = 1;
                
                while(<CUSTOM>) {
                    chomp;
                    
                    # check for errors
                    if($error_flag) {
                        die "\nERROR: Problem using annotation file ".$custom->{file}."\n$_\n" if /invalid pointer|tabix|get_intv/;
                        $error_flag = 0;
                    }
                    
                    my @data = split /\t/, $_;
                    
                    my $feature;
                    
                    if($custom->{format} eq 'bed') {
                        $feature = {
                            chr    => $chr,
                            start  => $data[1],
                            end    => $data[2],
                            name   => $data[3],
                        };
                    }
                    
                    elsif($custom->{format} eq 'vcf') {                        
                        my $tmp_vf = parse_vcf($config, $_)->[0];
                        
                        $feature = {
                            chr    => $chr,
                            start  => $tmp_vf->{start},
                            end    => $tmp_vf->{end},
                            name   => $tmp_vf->{variation_name},
                        };
                    }
                    
                    elsif($custom->{format} eq 'gff' || $custom->{format} eq 'gtf') {
                        
                        my $name;
                        
                        # try and get a feature name from the attributes column
                        foreach my $attrib(split /\s*\;\s*/, $data[8]) {
                            my ($key, $value) = split /\=/, $attrib;
                            $name = $value if $key eq 'ID';
                        }
                        
                        $name ||= $data[2]."_".$data[0].":".$data[3]."-".$data[4];
                        
                        $feature = {
                            chr   => $chr,
                            start => $data[3],
                            end   => $data[4],
                            name  => $name,
                        };
                    }
                    
                    elsif($custom->{format} eq 'bigwig') {
                        $feature = {
                            chr   => $chr,
                            start => $data[0],
                            end   => $data[0],
                            name  => $data[1],
                        };
                    }
                    
                    if(defined($feature)) {
                        $got_features = 1;
                        
                        if(!defined($feature->{name}) || $custom->{coords}) {
                            $feature->{name} = $feature->{chr}.":".$feature->{start}."-".$feature->{end};
                        }
                        
                        # add the feature to the cache
                        $annotation->{$chr}->{$custom->{name}}->{$feature->{start}}->{$feature->{name}} = $feature;
                    }
                }
                close CUSTOM;
                
                # unlink temporary wig files
                unlink($config->{tmpdir}.'/vep_tmp_'.$$) if $custom->{format} eq 'bigwig';
                
                # no need to fetch e.g. "chr21" features if just "21" worked
                last if $got_features;
            }
        }
    }
    
    end_progress($config);
    
    return $annotation;
}

# builds a full cache for this species
sub build_full_cache {
    my $config = shift;
    
    my @slices;
    
    if($config->{build} =~ /all/i) {
        @slices = @{$config->{sa}->fetch_all('toplevel')};
        push @slices, @{$config->{sa}->fetch_all('lrg', undef, 1, undef, 1)} if defined($config->{lrg});
    }
    else {
        foreach my $val(split /\,/, $config->{build}) {
            my @nnn = split /\-/, $val;
            
            foreach my $chr($nnn[0]..$nnn[-1]) {
                my $slice = get_slice($config, $chr);
                push @slices, $slice if defined($slice);
            }
        }
    }
    
    foreach my $slice(@slices) {
        my $chr = $slice->seq_region_name;
        
        # check for features, we don't want a load of effectively empty dirs
        my $dbc = $config->{sa}->db->dbc;
        my $sth = $dbc->prepare("SELECT COUNT(*) FROM transcript WHERE seq_region_id = ?");
        $sth->execute($slice->get_seq_region_id);
        
        my $count;
        $sth->bind_columns(\$count);
        $sth->fetch;
        $sth->finish;
        
        next unless $count > 0;
        
        my $regions;
        
        # for progress
        my $region_count = int($slice->end / $config->{cache_region_size}) + 1;
        my $counter = 0;
        
        # initial region
        my ($start, $end) = (1, $config->{cache_region_size});
        
        debug((defined($config->{rebuild}) ? "Rebuild" : "Creat")."ing cache for chromosome $chr") unless defined($config->{quiet});
        
        while($start < $slice->end) {
            
            progress($config, $counter++, $region_count);
            
            # store quiet status
            my $quiet = $config->{quiet};
            $config->{quiet} = 1;
            
            # spoof regions
            $regions->{$chr} = [$start.'-'.$end];
            
            # store transcripts
            my $tmp_cache = (defined($config->{rebuild}) ? load_dumped_transcript_cache($config, $chr, $start.'-'.$end) : cache_transcripts($config, $regions));
            $tmp_cache->{$chr} ||= [];
            
            #(defined($config->{tabix}) ? dump_transcript_cache_tabix($config, $tmp_cache, $chr, $start.'-'.$end) : dump_transcript_cache($config, $tmp_cache, $chr, $start.'-'.$end));
            dump_transcript_cache($config, $tmp_cache, $chr, $start.'-'.$end);
            undef $tmp_cache;
            
            # store reg feats
            if(defined($config->{regulatory})) {
                my $rf_cache = cache_reg_feats($config, $regions);
                $rf_cache->{$chr} ||= {};
                
                dump_reg_feat_cache($config, $rf_cache, $chr, $start.'-'.$end);
                #(defined($config->{tabix}) ? dump_reg_feat_cache_tabix($config, $rf_cache, $chr, $start.'-'.$end) : dump_reg_feat_cache($config, $rf_cache, $chr, $start.'-'.$end));
                undef $rf_cache;
                
                # this gets cleaned off but needs to be there for the next loop
                $slice->{coord_system}->{adaptor} = $config->{csa};
            }
            
            # store variations
            my $variation_cache;
            $variation_cache->{$chr} = get_variations_in_region($config, $chr, $start.'-'.$end);
            $variation_cache->{$chr} ||= {};
            
            dump_variation_cache($config, $variation_cache, $chr, $start.'-'.$end);
            undef $variation_cache;
            
            # restore quiet status
            $config->{quiet} = $quiet;
            
            # increment by cache_region_size to get next region
            $start += $config->{cache_region_size};
            $end += $config->{cache_region_size};
        }
        
        end_progress($config);
        
        undef $regions;
    }
    
    write_cache_info($config);
}

# write an info file that defines what is in the cache
sub write_cache_info {
    my $config = shift;
    
    my $info_file = $config->{dir}.'/info.txt';
    
    open OUT, ">>$info_file" or die "ERROR: Could not write to cache info file $info_file\n";
    
    print OUT "# CACHE UPDATED ".get_time()."\n";
    
    foreach my $param(qw(
        host
        port
        user
        build
        regulatory
        sift
        polyphen
    )) {
        print OUT "$param\t".(defined $config->{$param} ? $config->{$param} : '-')."\n";
    }
    
    # cell types
    if(defined($config->{cell_type}) && scalar(@{$config->{cell_type}})) {
        my $cta = $config->{RegulatoryFeature_adaptor}->db->get_CellTypeAdaptor();
        print OUT "cell_types\t".(join ",", map {$_->name} @{$cta->fetch_all});
        print OUT "\n";
    }
    
    # sift/polyphen versions
    foreach my $tool(qw(sift polyphen)) {
        if(defined($config->{$tool})) {
            my $var_mca = $config->{reg}->get_adaptor($config->{species}, 'variation', 'metacontainer');
            
            if(defined($var_mca)) {
                my $version = $var_mca->list_value_by_key($tool.'_version');
                print OUT "$tool\_version\t".$version->[0]."\n" if defined($version) and scalar @$version;
            }
        }
    }
    
    close OUT;
}

# reads in cache info file
sub read_cache_info {
    my $config = shift;
    
    my $info_file = $config->{dir}.'/info.txt';
    
    open IN, $info_file or return 0;
    
    while(<IN>) {
        next if /^#/;
        chomp;
        my ($param, $value) = split /\t/;
        $config->{'cache_'.$param} = $value unless defined $value && $value eq '-';
    }
    
    close IN;
    
    return 1;
}

# format coords for printing
sub format_coords {
    my ($start, $end) = @_;
    
    if(defined($start)) {
        if(defined($end)) {
            if($start > $end) {
                return $end.'-'.$start;
            }
            elsif($start == $end) {
                return $start;
            }
            else {
                return $start.'-'.$end;
            }
        }
        else {
            return $start.'-?';
        }
    }
    elsif(defined($end)) {
        return '?-'.$end;
    }
    else  {
        return '-';
    }
}




# METHODS TO FIND CO-LOCATED / EXISTING VARIATIONS
##################################################

# finds an existing VF in the db
sub find_existing {
    my $config = shift;
    my $new_vf = shift;
    
    if(defined($config->{vfa}->db)) {
        
        my $maf_cols = have_maf_cols($config) ? 'vf.minor_allele, vf.minor_allele_freq' : 'NULL, NULL';
        
        my $sth = $config->{vfa}->db->dbc->prepare(qq{
            SELECT variation_name, IF(fv.variation_id IS NULL, 0, 1), seq_region_start, seq_region_end, allele_string, seq_region_strand, $maf_cols
            FROM variation_feature vf LEFT JOIN failed_variation fv
    	    ON vf.variation_id = fv.variation_id
            WHERE vf.seq_region_id = ?
            AND vf.seq_region_start = ?
            AND vf.seq_region_end = ?
            ORDER BY vf.source_id ASC
        });
        
        $sth->execute($new_vf->slice->get_seq_region_id, $new_vf->start, $new_vf->end);
        
        my @v;
        for my $i(0..7) {
            $v[$i] = undef;
        }
        
        $sth->bind_columns(\$v[0], \$v[1], \$v[2], \$v[3], \$v[4], \$v[5], \$v[6], \$v[7]);
        
        my @found;
        
        while($sth->fetch) {
            push @found, $v[0] unless is_var_novel($config, \@v, $new_vf) || $v[1] > $config->{failed};
        }
        
        $sth->finish();
        
        return (scalar @found ? join ",", @found : undef);
    }
    
    return undef;
}

# compare a new vf to one from the cache / DB
sub is_var_novel {
    my $config = shift;
    my $existing_var = shift;
    my $new_var = shift;
    
    my $is_novel = 1;
    
    $is_novel = 0 if $existing_var->[2] == $new_var->start && $existing_var->[3] == $new_var->end;
    
    if(defined($config->{check_alleles})) {
        my %existing_alleles;
        
        $existing_alleles{$_} = 1 for split /\//, $existing_var->[4];
        
        my $seen_new = 0;
        foreach my $a(split /\//, ($new_var->allele_string || "")) {
            reverse_comp(\$a) if $new_var->strand ne $existing_var->[5];
            $seen_new = 1 unless defined $existing_alleles{$a};
        }
        
        $is_novel = 1 if $seen_new;
    }
    
    return $is_novel;
}

# check frequencies of existing var against requested params
sub check_frequencies {
    my $config = shift;
    my $var_name = shift;
    
    my $v = $config->{va}->fetch_by_name($var_name);
    
    my $pass = 0;
    
    my $freq_pop      = $config->{freq_pop};
    my $freq_freq     = $config->{freq_freq};
    my $freq_gt_lt    = $config->{freq_gt_lt};
    
    my $freq_pop_name = (split /\_/, $freq_pop)[-1];
    $freq_pop_name = undef if $freq_pop_name =~ /1kg|hap/;
    
    delete $config->{filtered_freqs};
    
    foreach my $a(@{$v->get_all_Alleles}) {
        next unless defined $a->{population} || defined $a->{'_population_id'};
        next unless defined $a->frequency;
        next if $a->frequency > 0.5;
        
        my $pop_name = $a->population->name;
        
        if($freq_pop =~ /1kg/) { next unless $pop_name =~ /^1000.+(low|phase).+/i; }
        if($freq_pop =~ /hap/) { next unless $pop_name =~ /^CSHL-HAPMAP/i; }
        if($freq_pop =~ /any/) { next unless $pop_name =~ /^(CSHL-HAPMAP)|(1000.+(low|phase).+)/i; }
        if(defined $freq_pop_name) { next unless $pop_name =~ /$freq_pop_name/i; }
        
        $pass = 1 if $a->frequency >= $freq_freq and $freq_gt_lt eq 'gt';
        $pass = 1 if $a->frequency <= $freq_freq and $freq_gt_lt eq 'lt';
        
        $pop_name =~ s/\:/\_/g;
        push @{$config->{filtered_freqs}}, $pop_name.':'.$a->frequency;
        
        #warn "Comparing allele ", $a->allele, " ", $a->frequency, " for $var_name in population ", $a->population->name, " PASS $pass";
    }
    
    return 0 if $config->{freq_filter} eq 'exclude' and $pass == 1;
    return 0 if $config->{freq_filter} eq 'include' and $pass == 0;
    return 1;
}

# gets all variations in a region
sub get_variations_in_region {
    my $config = shift;
    my $chr = shift;
    my $region = shift;
    
    my ($start, $end) = split /\-/, $region;
    
    my %variations;
    
    if(defined($config->{vfa}->db)) {
        my $sr_cache = $config->{seq_region_cache};
        
        if(!defined($sr_cache)) {
            $sr_cache = cache_seq_region_ids($config);
            $config->{seq_region_cache} = $sr_cache;
        }
        
        # no seq_region_id?
        return {} unless defined($sr_cache) && defined($sr_cache->{$chr});
        
        my $maf_cols = have_maf_cols($config) ? 'vf.minor_allele, vf.minor_allele_freq' : 'NULL, NULL';
        
        my $sth = $config->{vfa}->db->dbc->prepare(qq{
            SELECT vf.variation_name, IF(fv.variation_id IS NULL, 0, 1), vf.seq_region_start, vf.seq_region_end, vf.allele_string, vf.seq_region_strand, $maf_cols
            FROM variation_feature vf
            LEFT JOIN failed_variation fv ON fv.variation_id = vf.variation_id
            WHERE vf.seq_region_id = ?
            AND vf.seq_region_start >= ?
            AND vf.seq_region_start <= ?
        });
        
        $sth->execute($sr_cache->{$chr}, $start, $end);
        
        my @v;
        for my $i(0..7) {
            $v[$i] = undef;
        }
        
        $sth->bind_columns(\$v[0], \$v[1], \$v[2], \$v[3], \$v[4], \$v[5], \$v[6], \$v[7]);
        
        while($sth->fetch) {
            my @v_copy = @v;
            push @{$variations{$v[2]}}, \@v_copy;
        }
        
        $sth->finish();
    }
    
    return \%variations;
}

sub cache_seq_region_ids {
    my $config = shift;
    
    my (%cache, $chr, $id);
    
    my $sth = $config->{vfa}->db->dbc->prepare(qq{
        SELECT seq_region_id, name FROM seq_region
    });
    
    $sth->execute();
    $sth->bind_columns(\$id, \$chr);
    $cache{$chr} = $id while $sth->fetch();
    $sth->finish;
    
    return \%cache;
}

sub have_maf_cols {
    my $config = shift;
    
    if(!defined($config->{have_maf_cols})) {
        my $sth = $config->{vfa}->db->dbc->prepare(qq{
            DESCRIBE variation_feature
        });
        
        $sth->execute();
        my @cols = map {$_->[0]} @{$sth->fetchall_arrayref};
        $sth->finish();
        
        $config->{have_maf_cols} = 0;
        $config->{have_maf_cols} = 1 if grep {$_ eq 'minor_allele'} @cols;
    }
    
    return $config->{have_maf_cols};
}

sub merge_hashes {
    my ($x, $y) = @_;

    foreach my $k (keys %$y) {
        if (!defined($x->{$k})) {
            $x->{$k} = $y->{$k};
        } else {
            if(ref($x->{$k}) eq 'ARRAY') {
                $x->{$k} = merge_arrays($x->{$k}, $y->{$k});
            }
            elsif(ref($x->{$k}) eq 'HASH') {
                $x->{$k} = merge_hashes($x->{$k}, $y->{$k});
            }
            else {
                $x->{$k} = $y->{$k};
            }
        }
    }
    return $x;
}

sub merge_arrays {
    my ($x, $y) = @_;
    
    my %tmp = map {$_ => 1} (@$x, @$y);
    
    return [keys %tmp];
}




# DEBUG AND STATUS METHODS
##########################

# gets time
sub get_time() {
    my @time = localtime(time());

    # increment the month (Jan = 0)
    $time[4]++;

    # add leading zeroes as required
    for my $i(0..4) {
        $time[$i] = "0".$time[$i] if $time[$i] < 10;
    }

    # put the components together in a string
    my $time =
        ($time[5] + 1900)."-".
        $time[4]."-".
        $time[3]." ".
        $time[2].":".
        $time[1].":".
        $time[0];

    return $time;
}

# prints debug output with time
sub debug {
    my $text = (@_ ? (join "", @_) : "No message");
    my $time = get_time;
    
    print $time." - ".$text.($text =~ /\n$/ ? "" : "\n");
}

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

# update or initiate progress bar
sub progress {
    my ($config, $i, $total) = @_;
    
    return if defined($config->{quiet}) || defined($config->{no_progress});
    
    my $width = $config->{terminal_width} || 60;
    my $percent = int(($i/$total) * 100);
    my $numblobs = int((($i/$total) * $width) - 2);
    
    # this ensures we're not writing to the terminal too much
    return if(defined($config->{prev_prog})) && $numblobs.'-'.$percent eq $config->{prev_prog};
    $config->{prev_prog} = $numblobs.'-'.$percent;
    
    printf("\r% -${width}s% 1s% 10s", '['.('=' x $numblobs).($numblobs == $width - 2 ? '=' : '>'), ']', "[ " . $percent . "% ]");
}

# end progress bar
sub end_progress {
    my $config = shift;
    return if defined($config->{quiet}) || defined($config->{no_progress});
    progress($config, 1,1);
    print "\n";
    delete $config->{prev_prog};
}

1;
