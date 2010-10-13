
package Bio::EnsEMBL::Variation::Utils::EnsEMBL2GFF3;

use strict;
use warnings;

# This module allows conversion of ensembl objects to GFF3 and GVF by inserting
# to_gff (and supporting _gff_hash) methods into the necessary feature classes

{
    package Bio::EnsEMBL::Slice;
    
    sub gff_header {
        my $self = shift;
      
        # build up a date string in the format specified by the GFF spec
    
        my ( $sec, $min, $hr, $mday, $mon, $year ) = localtime;
        $year += 1900;    # correct the year
        $mon++;           # correct the month
        
        my $date = "$year-$mon-$mday";
        
        my $region      = $self->seq_region_name;
        my $start       = $self->start;
        my $end         = $self->end;
        my $assembly    = $self->coord_system->version;
        
        my $hdr =
            "##gff-version 3\n"
          . "##date $date\n"
          . "##sequence-region $region $start $end\n"
          . "##genome-build $assembly\n";
    
        return $hdr;
    }
    
    sub gvf_header {
        my $self = shift;
        
        my $hdr = $self->gff_header(@_);
        
        $hdr .= "##gvf-version 1.02\n";
        $hdr .= "##feature-ontology http://song.cvs.sourceforge.net/viewvc/song/ontology/so.obo?revision=1.263\n";
    
        return $hdr;
    }
}

{

    package Bio::EnsEMBL::Feature;

    sub to_gff {

        my $self = shift;
        
        my %args = @_;
        
        # This parameter is assumed to be a hashref which includes extra attributes you'd
        # like to have appended onto the gff line for the feature
        my $extra_attrs = $args{extra_attrs};

        my $gff = $self->_gff_hash(@_);

        $gff->{score}  = '.' unless defined $gff->{score};
        $gff->{strand} = '.' unless defined $gff->{strand};
        $gff->{phase}  = '.' unless defined $gff->{phase};

        # order as per GFF3 spec: http://www.sequenceontology.org/gff3.shtml 

        my $gff_str = join( "\t",
            $gff->{seqid}, $gff->{source}, $gff->{type}, $gff->{start},
            $gff->{end}, $gff->{score}, $gff->{strand}, $gff->{phase},
        );

        if ($extra_attrs) {

            # combine the extra attributes with any existing ones (duplicate keys will get squashed!)
            $gff->{attributes} = {} unless defined $gff->{attributes};
            @{ $gff->{attributes} }{ keys %$extra_attrs } = values %$extra_attrs;
        }

        if ( $gff->{attributes} ) {

            my @attrs;

            for my $key (keys %{ $gff->{attributes} }) {

                my $val = $gff->{attributes}->{$key};

                if (ref $val eq 'ARRAY') {
                    push @attrs, map { $key . '='. $_ } @$val;
                }
                else {
                    push @attrs, $key . '=' . $val;
                }
            }
                
            $gff_str .= "\t" . join( ';', @attrs );
        }

        return $gff_str;
    }

    sub _gff_hash {
        my $self = shift;
        
        my %args = @_;
        
        my $rebase      = $args{rebase};
        my $gff_seqid   = $args{gff_seqid} || $self->slice->seq_region_name;
        my $gff_source  = $args{gff_source} || $self->_gff_source;

        my $seqid = $rebase ? $gff_seqid.'_'.$self->slice->start.'-'.$self->slice->end : $gff_seqid;
        my $start = $rebase ? $self->start : $self->seq_region_start;
        my $end = $rebase ? $self->end : $self->seq_region_end;

        my $gff = {
            seqid   => $gff_seqid,
            source  => $gff_source,
            type    => $self->_gff_type,
            start   => $start,
            end     => $end,
            strand  => (
                $self->strand == 1 ? '+' : ( $self->strand == -1 ? '-' : '.' )
            )
        };

        return $gff;
    }

    sub _gff_source {
        my $self = shift;
        
        if ($self->analysis) {
            return
                $self->analysis->gff_source
                || $self->analysis->logic_name;
        }
        else {
            return ref($self);
        }
    }

    sub _gff_type {
        my $self = shift;

        return
            ( $self->analysis && $self->analysis->gff_feature )
            || 'misc_feature';
    }
}

{
    package Bio::EnsEMBL::Variation::VariationFeature;

    my %EnsEMBL2SO_consequence = (
        'ESSENTIAL_SPLICE_SITE' => { 
            term            => 'splice_site_variant', 
            id              => 'SO:0001629', 
            feature_type    => 'mRNA',
        },
        'STOP_GAINED' => { 
            term            => 'stop_gained', 
            id              => 'SO:0001587',
            feature_type    => 'mRNA',
        },
        'STOP_LOST' => { 
            term            => 'stop_lost', 
            id              => 'SO:0001578',
            feature_type    => 'mRNA',
        },
        'COMPLEX_INDEL' => { 
            term            => 'Complex_change_in_transcript', 
            id              => 'SO:0001577',
            feature_type    => 'mRNA',
        },
        'FRAMESHIFT_CODING' => { 
            term            => 'Frameshift_variant', 
            id              => 'SO:0001589', 
            feature_type    => 'mRNA',
        },
        'NON_SYNONYMOUS_CODING' => { 
            term            => 'Non_synonymous_codon', 
            id              => 'SO:0001583', 
            feature_type    => 'mRNA',
        },
        'SPLICE_SITE' => { 
            term            => 'splice_region_variant', 
            id              => 'SO:0001630', 
            feature_type    => 'mRNA',
        },
        'PARTIAL_CODON' => { 
            term            => 'incomplete_terminal_codon_variant', 
            id              => 'SO:0001626',
            feature_type    => 'mRNA',
        },
        'SYNONYMOUS_CODING' => { 
            term            => 'Synonymous_codon', 
            id              => 'SO:0001588',
            feature_type    => 'mRNA',
        },
        'REGULATORY_REGION' => { 
            term            => 'regulatory_region_variant', 
            id              => 'SO:0001566',
            feature_type    => 'genomic',
        },
        'WITHIN_MATURE_miRNA' => { 
            term            => 'Mature_miRNA_variant', 
            id              => 'SO:0001620',
            feature_type    => 'miRNA',
        },
        '5PRIME_UTR' => { 
            term            => '5_prime_UTR_variant', 
            id              => 'SO:0001623',
            feature_type    => 'mRNA'
        },
        '3PRIME_UTR' => { 
            term            => '3_prime_UTR_variant', 
            id              => 'SO:0001624',
            feature_type    => 'mRNA',
        },
        'UTR' => {},
        'INTRONIC' => { 
            term            => 'Intron_variant', 
            id              => 'SO:0001627',
            feature_type    => 'mRNA',
        },
        'NMD_TRANSCRIPT' => { 
            term            => 'NMD_transcript_variant', 
            id              => 'SO:0001621',
            feature_type    => 'mRNA',
        },
        'WITHIN_NON_CODING_GENE' => { 
            term            => 'Non_coding_RNA_variant', 
            id              => 'SO:0001619',
            feature_type    => 'ncRNA'
        },
        'UPSTREAM' => { 
            term            => '5KB_upstream_variant', 
            id              => 'SO:0001635',
            feature_type    => 'mRNA',
        },
        'DOWNSTREAM' => { 
            term            => '5KB_downstream_variant',
            id              => 'SO:0001633',
            feature_type    => 'mRNA',
        },
        'HGMD_MUTATION' => {},
        'NO_CONSEQUENCE' => { 
            term            => 'Sequence_variant', 
            id              => 'SO:0001060',
            feature_type    => 'genomic',
        },
        'INTERGENIC' => { 
            term            => 'Intergenic_variant', 
            id              => 'SO:0001628',
            feature_type    => 'genomic',
        },
    );
    
    sub to_gvf {
        my $self = shift;
        return $self->to_gff(@_);
    }
    
    sub _gff_hash {
        my $self = shift;
        
        my $gff = $self->SUPER::_gff_hash(@_);

        $gff->{source} = $self->source;
        
        $gff->{type} = $self->var_class;

        # Use the variation name (rsID etc.) as the ID
        
        my $gvf_id = $self->variation_name;

        $gff->{attributes}->{ID} = $gvf_id;

        $gff->{attributes}->{Variant_effect} = [];

        # CNV probes and HGMD mutations are special cases, so deal with them first

        if ($self->allele_string eq 'CNV_PROBE') {
            
            $gff->{attributes}->{Reference_seq} = '~';
            $gff->{attributes}->{Variant_seq} = '~';
            
            for my $ens_cons (@{ $self->get_consequence_type }) {
                my $effect          = $EnsEMBL2SO_consequence{$ens_cons}->{term};
                my $feature_type    = $EnsEMBL2SO_consequence{$ens_cons}->{feature_type};
                my $features        = $gff->{seqid};
                        
                my $effect_string = join ' ', ($effect, 0, $feature_type, $features);
                push @{ $gff->{attributes}->{Variant_effect} }, $effect_string;
            }
            
            return $gff;
        }
        
        if ($self->allele_string eq 'HGMD_MUTATION') {
            
            $gff->{attributes}->{Reference_seq} = '~';
            $gff->{attributes}->{Variant_seq} = '~';
            
            # there's not really anything useful to put in the Variant_effect 
            # attribute, so we don't add one
            
            return $gff;
        }
        
        # the Variant_seq attribute requires a comma separated list of alleles

        my @alleles = split '/', $self->allele_string;
        my $ref_seq = shift @alleles; # shift off the reference allele
        
        $gff->{attributes}->{Variant_seq} = join ',', @alleles;

        # the reference sequence should be set to '~' if the sequence is longer than 50 nucleotides

        $ref_seq = '~' if length($ref_seq) > 50;
        $gff->{attributes}->{Reference_seq} = $ref_seq;

        my $allele_index = 0;
        
        for my $allele (@alleles) {
            
            my @tvs = @{ $self->get_all_TranscriptVariations };
            
            if (@tvs) {
                
                # this variation affects some transcripts
        
                for my $tv (@tvs) {
                    
                    for my $ens_cons (@{ $tv->consequence_type }) {
                    
                        my $effect          = $EnsEMBL2SO_consequence{$ens_cons}->{term};
                        my $feature_type    = $EnsEMBL2SO_consequence{$ens_cons}->{feature_type};
                        my $features        = $tv->transcript->stable_id;
                        
                        my $effect_string = join ' ', ($effect, $allele_index, $feature_type, $features);
                
                        push @{ $gff->{attributes}->{Variant_effect} }, $effect_string;
                    }
                }
            }
            else { 
                
                # no TranscriptVariations, so the effect of all alleles will just be the
                # effect of the VariationFeature
                
                for my $ens_cons (@{ $self->get_consequence_type }) {
                    
                    my $effect          = $EnsEMBL2SO_consequence{$ens_cons}->{term};
                    my $feature_type    = $EnsEMBL2SO_consequence{$ens_cons}->{feature_type};
                    # XXX: I guess the feature affected by this variation is the slice?
                    my $features        = $gff->{seqid};
                    
                    my $effect_string = join ' ', ($effect, $allele_index, $feature_type, $features);
                
                    push @{ $gff->{attributes}->{Variant_effect} }, $effect_string;
                }
            }
            
            $allele_index++;
        }
        
        return $gff;
    }
}

{
    package Bio::EnsEMBL::Variation::StructuralVariation;
    
    sub _gff_hash {
        
        my $self = shift;
        
        my $gff = $self->SUPER::_gff_hash(@_);
        
        $gff->{attributes}->{ID} = $self->variation_name;
        
        $gff->{source} = $self->source;
        
        $gff->{type} = $self->class;
        
        $gff->{attributes}->{Reference_seq} = $self->end > $self->start+50 ? '~' : $self->get_reference_sequence;
        
        return $gff;
    }
    
    sub to_gvf {
        my $self = shift;
        return $self->to_gff(@_);
    }
    
}


1;


