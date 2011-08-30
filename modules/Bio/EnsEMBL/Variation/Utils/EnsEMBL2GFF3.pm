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
        
        my %args = @_;
        
        # build up a date string in the format specified by the GFF spec
    
        my ( $sec, $min, $hr, $mday, $mon, $year ) = localtime;
        $year += 1900;    # correct the year
        $mon++;           # correct the month
        
        my $date = sprintf "%4d-%02d-%02d", $year, $mon, $mday;
        
        my $region      = $self->seq_region_name;
        my $start       = $self->start;
        my $end         = $self->end;
        my $assembly    = $self->coord_system->version;
        
        my $hdr =
            "##gff-version 3\n"
          . "##file-date $date\n"
          . "##genome-build ensembl $assembly\n";
        
        $hdr .= "##sequence-region $region $start $end\n" unless $args{no_sequence_region};
        
        return $hdr;
    }
    
    sub gvf_header {
        my $self = shift;
       
        my %args = @_;

        my $hdr = $self->gff_header(@_);
        
        my $mca = $self->adaptor->db->get_MetaContainerAdaptor;
        my $schema_version = $mca->get_schema_version;
        my $species_name = $mca->get_scientific_name;
        $species_name =~ s/ /_/g;
        my $url = 'http://e'.$schema_version.'.ensembl.org/'.$species_name;
        
        $hdr .= "##gvf-version 1.05\n";
        $hdr .= "##feature-ontology http://song.cvs.sourceforge.net/viewvc/song/ontology/so.obo?revision=1.283\n";
        $hdr .= "##data-source Source=ensembl;version=$schema_version;url=$url\n";
        $hdr .= "##file-version $schema_version\n";
        
        if (my $individual = $args{individual}) {
            $hdr .= "##individual-id ".$individual->to_gvf."\n";
        }

        if (my $population = $args{population}) {
            $hdr .= "##attribute-method ".$population->to_gvf."\n";
        }

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

        return undef unless defined $gff;

        # default optional columns, and check that all required columns are present

        $gff->{score}  = '.' unless defined $gff->{score};
        $gff->{strand} = '.' unless defined $gff->{strand};
        $gff->{phase}  = '.' unless defined $gff->{phase};

        for my $req (qw(source type start end)) {
            die "'$req' attribute required for GFF" unless $gff->{$req};
        }

        # order as per GFF3 spec: http://www.sequenceontology.org/gff3.shtml 
 
        my $gff_str = join( "\t",
            $gff->{seqid}, $gff->{source}, $gff->{type}, $gff->{start},
            $gff->{end}, $gff->{score}, $gff->{strand}, $gff->{phase},
        );

        if ($extra_attrs) {

            # combine the extra attributes with any existing ones (duplicate keys will get squashed,
            # so attributes specified in the extra_attrs hash will override existing ones)

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
        
        my $rebase      = $args{rebase}; # use absolute or slice-relative coordinates
        
        my $gff_seqid   = $args{gff_seqid} || $self->slice->seq_region_name;
        my $gff_source  = $args{gff_source} || $self->_gff_source;

        my $seqid = $rebase ? $gff_seqid.'_'.$self->slice->start.'-'.$self->slice->end : $gff_seqid;
        my $start = $rebase ? $self->start : $self->seq_region_start;
        my $end = $rebase ? $self->end : $self->seq_region_end;
        
        # GFF3 does not allow start > end, and mandates that for zero-length features (e.g. insertions) 
        # start = end and the implied insertion site is to the right of the specified base, so we use the
        # smaller of the two values
        
        if ($start > $end) {
            $start = $end;
        }
        
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
    
    use Bio::EnsEMBL::Utils::Sequence qw(expand);
    
    sub to_gvf {
        my $self = shift;
        return $self->to_gff(@_);
    }
    
    sub _gff_hash {
   
        my $self = shift;
        
        my $gff = $self->SUPER::_gff_hash(@_);
        
        my %args = @_;
        
        my $include_consequences = $args{include_consequences};

        $gff->{source} = $self->source;

        $gff->{type} = $self->class_SO_term;

        my $source = $self->source;

        $source .= '_'.$self->source_version if defined $self->source_version;

        $gff->{attributes}->{Dbxref} = "$source:".$self->variation_name;

        $gff->{attributes}->{ID} = $self->dbID;

        $gff->{attributes}->{Variant_effect} = [];

        # the Variant_seq attribute requires a comma separated list of alleles

        my @alleles = split '/', $self->allele_string;
        my $ref_seq = shift @alleles; # shift off the reference allele
       
        $gff->{attributes}->{Variant_seq} = join ',', @alleles;

        my $index = 0;

        # expand tandem repeat alleles, because TranscriptVariationAlleles use the expanded sequence
        
        map { expand(\$_) } @alleles; 

        # if you expand e.g. (T)0 you get an empty string, which we treat as a deletion, so default to '-'

        my %allele_index = map { ($_ || '-') => $index++ } @alleles;
        
        # the reference sequence should be set to '~' if the sequence is longer than 50 nucleotides

        $ref_seq = '~' if CORE::length($ref_seq) > 50;
        $gff->{attributes}->{Reference_seq} = $ref_seq;
        
        # Hack for HGMD mutations
       
        if ($self->allele_string eq 'HGMD_MUTATION') {
            $gff->{attributes}->{Reference_seq} = '~';
            $gff->{attributes}->{Variant_seq}   = '~';
            $allele_index{$self->allele_string} = 0;
        }
        
        if ($include_consequences) {
            for my $tv (@{ $self->get_all_TranscriptVariations }) {

                unless ($tv->get_all_alternate_TranscriptVariationAlleles) {
                    warn $self->variation_name." has no alternate alleles?";
                    next;
                }

                for my $tva (@{ $tv->get_all_alternate_TranscriptVariationAlleles }) {

                    my $allele_idx = $allele_index{$tva->variation_feature_seq};
                    
                    if (defined $allele_idx) {
                        for my $oc (@{ $tva->get_all_OverlapConsequences }) {

                            my $effect_string = join ' ', (
                                $oc->SO_term, 
                                $allele_idx,
                                $oc->feature_SO_term,
                                $tv->transcript_stable_id,
                            );

                            push @{ $gff->{attributes}->{Variant_effect} }, $effect_string;
                        }
                    }
                    else {
                        warn "No allele_index entry for allele: ".$tva->variation_feature_seq.
                            " of ".$self->variation_name."?\n";
                    }
                }
            }
        }
        
        return $gff;
    }
}

{
    package Bio::EnsEMBL::Variation::StructuralVariationFeature;
    
    sub _gff_hash {
        
        my $self = shift;
        
        my $gff = $self->SUPER::_gff_hash(@_);
        
        $gff->{attributes}->{ID} = $self->variation_name;
        
        $gff->{source} = $self->source;
        
        my $sv = $self->structural_variation;
        
        $gff->{attributes}->{Dbxref} = $self->source . ':' . $self->variation_name;
        
        $gff->{attributes}->{study_accession} = $sv->study->name if $sv->study->name;

        $gff->{type} = $self->class_SO_term;
        
        #$gff->{attributes}->{Reference_seq} = $self->end > $self->start+50 ? '~' : $self->get_reference_sequence;
     
        if ( (defined $self->inner_start) && (defined $self->outer_start) && ($self->inner_start != $self->outer_start) ) {
            $gff->{attributes}->{Start_range} = join ',', $self->outer_start, $self->inner_start;
        }

        if ( (defined $self->inner_end) && (defined $self->outer_end) && ($self->inner_end != $self->outer_end) ) {
            $gff->{attributes}->{End_range} = join ',', $self->inner_end, $self->outer_end;
        }

        return $gff;
    }
    
    sub to_gvf {
        my $self = shift;
        return $self->to_gff(@_);
    }
}

{
    package Bio::EnsEMBL::Variation::Individual;
        
    sub _gff_hash {
        
        my $self = shift;

        my $gff;

        $gff->{Gender} = $self->gender;

        $gff->{Display_name} = $self->name;

        $gff->{ensembl_description} = $self->description;

        $gff->{Type} = $self->type_description;

        $gff->{Population} = join ',', map { $_->name } @{ $self->get_all_Populations };
       
        return $gff;
    }

    sub to_gvf {
        my $self = shift;
       
        my $attrs = $self->_gff_hash(@_);

        # get rid of any empty attributes 
        map { delete $attrs->{$_} unless $attrs->{$_} } keys %$attrs;
        
        return join ';', map { $_.'='.$attrs->{$_} } keys %$attrs; 
    }
    
}

{
    package Bio::EnsEMBL::Variation::Population;
        
    sub _gff_hash {
        
        my $self = shift;

        my $gff;

        $gff->{Attribute} = 'Variant_freq';

        $gff->{population} = $self->name;
        
        $gff->{population_size} = $self->size;

        $gff->{Comment} = $self->description;
      
        return $gff;
    }

    sub to_gvf {
        my $self = shift;
       
        my $attrs = $self->_gff_hash(@_);

        # get rid of any empty attributes 
        map { delete $attrs->{$_} unless $attrs->{$_} } keys %$attrs;
        
        return join ';', map { $_.'='.$attrs->{$_} } keys %$attrs; 
    }
    
}

1;

