=head1 LICENSE

 Copyright (c) 1999-2011 The European Bioinformatics Institute and
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
    
    use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(transcript_effect);
    use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
    
    sub to_gvf {
        my $self = shift;
        return $self->to_gff(@_);
    }
    
    sub _gff_hash {
   
        my $self = shift;
        
        my $gff = $self->SUPER::_gff_hash(@_);
        
        my %args = @_;
        
        my $consequences = $args{consequences};

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
            
            # Pontus says I can't meaningfully calculate the consequence of these variations
            
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

        my $index = 0;
        my %allele_index = map { $_ => $index++ } @alleles;

        # the reference sequence should be set to '~' if the sequence is longer than 50 nucleotides

        $ref_seq = '~' if length($ref_seq) > 50;
        $gff->{attributes}->{Reference_seq} = $ref_seq;
        
        return $gff;

        if (my $db = $self->{adaptor}->db) {
            print "got a dba\n";
            my $tvna = $db->get_TranscriptVariationNewAdaptor;
            print "got a tvna\n" if $tvna;
            my $tvns = $tvna->fetch_all_by_VariationFeatures([$self]);
            print "Got ".scalar(@$tvns)." tvns\n";
            die;
        }

        my @transcripts = map { $_->transcript } @{ $self->get_all_TranscriptVariations };
        
        for my $transcript (@transcripts) {
        
            my $tv = transcript_effect($self, $transcript, $consequences);
           
            if ($tv) {

                for my $tva (@{ $tv->alt_alleles }) {

                    for my $cons (@{ $tva->consequences }) {

                        my $allele = $tva->seq;

                        reverse_comp(\$allele) unless $self->strand == $transcript->strand;

                        my $effect_string = join ' ', (
                            $cons->SO_term, 
                            $allele_index{$allele},
                            $cons->feature_SO_term,
                            $transcript->stable_id,
                        );
    
                        push @{ $gff->{attributes}->{Variant_effect} }, $effect_string;
                    }
                }
            }
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


