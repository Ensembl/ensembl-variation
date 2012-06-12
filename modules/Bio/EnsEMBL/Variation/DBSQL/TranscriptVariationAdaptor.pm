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

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor

=head1 SYNOPSIS
    my $reg = 'Bio::EnsEMBL::Registry';
  
    $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
    my $tva = $reg->get_adaptor('human','variation','TranscriptVariation');
  
    my $ta = $reg->get_adaptor('human','core','Transcript');
    my $va = $reg->get_adaptor('human','variation','Variation');
    my $vfa = $reg->get_adaptor('human','variation','VariationFeature');

    # fetch all TranscriptVariations related to a Transcript
    my $tran = $ta->fetch_by_stable_id('ENST00000380152');

    for my $tv (@{ $tva->fetch_all_by_Transcripts([$tran]) }) {
        print $tv->consequence_type, "\n";
        print $tv->cdna_start, '-', $tv->cdna_end, "\n";
    }
    
    # fetch all TranscriptVariations related to a VariationFeature
    my $vf = $vfa->fetch_all_by_Variation($va->fetch_by_name('rs669'))->[0];

    for my $tv (@{ $tva->fetch_all_by_VariationFeatures([$vf]) }) {
        print $tv->transcript->stable_id, "\n";
        print $tv->translation_start, '-', $tv->translation_end, "\n";
    }
    
=head1 DESCRIPTION

This adaptor allows you to fetch TranscriptVariation objects either by the Transcripts
the associated VariationFeature falls in, or by VariationFeature directly. Storing
TranscriptVariation objects in a variation schema database is also supported. In the
database there will a separate row for each alternative allele of a TranscriptVariation, 
but the methods here will fetch all alleles associated with the TranscriptVariation
at once.

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor;

use Bio::EnsEMBL::Variation::TranscriptVariation;
use Bio::EnsEMBL::Variation::TranscriptVariationAllele;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(MAX_DISTANCE_FROM_TRANSCRIPT);
use Bio::EnsEMBL::Variation::Utils::Constants qw(%OVERLAP_CONSEQUENCES);
use Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix;

use base qw(Bio::EnsEMBL::Variation::DBSQL::VariationFeatureOverlapAdaptor);

=head2 store

  Arg [1]    : Bio::EnsEMBL::Variation::TranscriptVariation $tv
  Description: Store the TranscriptVariation in the database
  Status     : At risk

=cut

sub store {
    my ($self, $tv) = @_;
    
    my $dbh = $self->dbc->db_handle;
    
    my $sth = $dbh->prepare_cached(q{
        INSERT DELAYED INTO transcript_variation (
            variation_feature_id,
            feature_stable_id,
            allele_string,
            somatic,
            consequence_types,
            cds_start,
            cds_end,
            cdna_start,
            cdna_end,
            translation_start,
            translation_end,
            distance_to_transcript,
            codon_allele_string,
            pep_allele_string,
            hgvs_genomic,
            hgvs_transcript,
            hgvs_protein,
            polyphen_prediction,
            polyphen_score,
            sift_prediction,
            sift_score
        ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
    });

    for my $allele (@{ $tv->get_all_alternate_TranscriptVariationAlleles }) {
        
        $sth->execute(
            $tv->variation_feature->dbID,
            $tv->feature->stable_id,
            $allele->allele_string,
            $tv->variation_feature->is_somatic,
            (join ',', map { $_->SO_term } @{ $allele->get_all_OverlapConsequences }),
            $tv->cds_start, 
            $tv->cds_end,
            $tv->cdna_start,
            $tv->cdna_end,
            $tv->translation_start,
            $tv->translation_end,
            $tv->distance_to_transcript,
            $allele->codon_allele_string,
            $allele->pep_allele_string,
            $allele->hgvs_genomic,
            $allele->hgvs_transcript,
            $allele->hgvs_protein,
            $allele->polyphen_prediction,
            $allele->polyphen_score,
            $allele->sift_prediction,
            $allele->sift_score
        );
    }
}

sub fetch_all_by_Transcripts_SO_terms {
    my ($self, $transcripts, $terms) = @_;
    my $constraint = $self->_get_consequence_constraint(@$terms);
    return $self->fetch_all_by_Transcripts_with_constraint($transcripts, $constraint.' AND somatic = 0');
}

sub fetch_all_somatic_by_Transcripts_SO_terms {
    my ($self, $transcripts, $terms) = @_;
    my $constraint = $self->_get_consequence_constraint(@$terms);
    return $self->fetch_all_by_Transcripts_with_constraint($transcripts, $constraint.' AND somatic = 1');
}

=head2 fetch_all_by_Transcripts

  Arg [1]    : listref of Bio::EnsEMBL::Transcripts
  Description: Fetch all germline TranscriptVariations associated with the
               given list of Transcripts
  Returntype : listref of Bio::EnsEMBL::Variation::TranscriptVariations
  Status     : At risk

=cut

sub fetch_all_by_Transcripts {
    my ($self, $transcripts) = @_;
    return $self->fetch_all_by_Transcripts_with_constraint($transcripts, 'somatic = 0');
}

=head2 fetch_all_somatic_by_Transcripts

  Arg [1]    : listref of Bio::EnsEMBL::Transcripts
  Description: Fetch all somatic TranscriptVariations associated with the
               given list of Transcripts
  Returntype : listref of Bio::EnsEMBL::Variation::TranscriptVariations
  Status     : At risk

=cut

sub fetch_all_somatic_by_Transcripts {
    my ($self, $transcripts) = @_;
    return $self->fetch_all_by_Transcripts_with_constraint($transcripts, 'somatic = 1');
}

=head2 fetch_all_by_Transcripts_with_constraint

  Arg [1]    : listref of Bio::EnsEMBL::Transcripts
  Arg [2]    : extra SQL constraint for the query
  Description: Fetch all TranscriptVariations associated with the
               given list of Transcripts
  Returntype : listref of Bio::EnsEMBL::Variation::TranscriptVariations
  Status     : At risk

=cut

sub fetch_all_by_Transcripts_with_constraint {
    my ($self, $transcripts, $constraint) = @_;
    
    return $self->SUPER::fetch_all_by_Features_with_constraint($transcripts, $constraint);
}

sub _objs_from_sth {
    my ($self, $sth) = @_;

    #warn $sth->sql;

    my (
        $transcript_variation_id,
        $variation_feature_id, 
        $feature_stable_id, 
        $allele_string,
        $consequence_types,
        $cds_start,
        $cds_end,
        $cdna_start,
        $cdna_end,
        $translation_start,
        $translation_end,
        $distance_to_transcript,
        $codon_allele_string,
        $pep_allele_string,
        $hgvs_genomic,
        $hgvs_transcript,
        $hgvs_protein,
        $polyphen_prediction,
        $polyphen_score,
        $sift_prediction,
        $sift_score,
    );
    
    $sth->bind_columns(
        \$transcript_variation_id,
        \$variation_feature_id, 
        \$feature_stable_id, 
        \$allele_string,
        \$consequence_types,
        \$cds_start,
        \$cds_end,
        \$cdna_start,
        \$cdna_end,
        \$translation_start,
        \$translation_end,
        \$distance_to_transcript,
        \$codon_allele_string,
        \$pep_allele_string,
        \$hgvs_genomic,
        \$hgvs_transcript,
        \$hgvs_protein,
        \$polyphen_prediction,
        \$polyphen_score,
        \$sift_prediction,
        \$sift_score,
    );
    
    my %tvs;
    
    while ($sth->fetch) {
        
        my ($ref_allele, $alt_allele)   = split /\//, $allele_string;
        my ($ref_codon, $alt_codon)     = split /\//, $codon_allele_string || '';
        my ($ref_pep, $alt_pep)         = split /\//, $pep_allele_string || '';
        
        # for HGMD mutations etc. just set the alt allele to the ref allele
        $alt_allele ||= $ref_allele;
        
        # for synonymous mutations the peptides are the same and 
        # there is no / in the string
        $alt_pep ||= $ref_pep;
        
        # for TranscriptVariations with multiple alternative alleles
        # there will be multiple rows in the database, so we construct
        # the TV object and the reference allele object when we see 
        # the first row, but then only add extra allele objects when 
        # we see further rows, we track existing TVs in the %tvs hash, 
        # keyed by variation_feature_id and feature_stable_id

        my $key = $variation_feature_id.'_'.$feature_stable_id;
        
        my $tv = $tvs{$key};
        
        unless ($tv) {
            $tv = Bio::EnsEMBL::Variation::TranscriptVariation->new_fast({
                _variation_feature_id   => $variation_feature_id,
                _feature_stable_id      => $feature_stable_id,
                cds_start               => $cds_start,
                cds_end                 => $cds_end,
                cdna_start              => $cdna_start,
                cdna_end                => $cdna_end,
                translation_start       => $translation_start,
                translation_end         => $translation_end,
                distance_to_transcript  => $distance_to_transcript,
                adaptor                 => $self,
            });
            
            $tvs{$key} = $tv;
            
            my $ref_allele = Bio::EnsEMBL::Variation::TranscriptVariationAllele->new_fast({
                is_reference                => 1,
                variation_feature_seq       => $ref_allele,
                transcript_variation        => $tv, 
                codon                       => $ref_codon,
                peptide                     => $ref_pep, 
                dbID                        => $transcript_variation_id,
            });

            $tv->add_TranscriptVariationAllele($ref_allele);
        }
       
        #my $overlap_consequences = $self->_transcript_variation_consequences_for_set_number($consequence_types);

        my $overlap_consequences = [ map { $OVERLAP_CONSEQUENCES{$_} } split /,/, $consequence_types ];
        
        my $allele = Bio::EnsEMBL::Variation::TranscriptVariationAllele->new_fast({
            is_reference                => 0,
            variation_feature_seq       => $alt_allele,
            transcript_variation        => $tv, 
            codon                       => $alt_codon,
            peptide                     => $alt_pep,
            hgvs_genomic                => $hgvs_genomic,
            hgvs_transcript             => $hgvs_transcript,
            hgvs_protein                => $hgvs_protein,
            overlap_consequences        => $overlap_consequences, 
            polyphen_prediction         => $polyphen_prediction,
            polyphen_score              => $polyphen_score,
            sift_prediction             => $sift_prediction, 
            sift_score                  => $sift_score, 
            dbID                        => $transcript_variation_id,
        });
        
        $tv->add_TranscriptVariationAllele($allele);
    }
    
    return [values %tvs];
}

sub _tables {
    return (
        ['transcript_variation', 'tv']
    );
}

sub _columns {
    return qw(
        transcript_variation_id 
        variation_feature_id 
        feature_stable_id 
        allele_string 
        consequence_types 
        cds_start 
        cds_end 
        cdna_start 
        cdna_end 
        translation_start 
        translation_end 
        distance_to_transcript 
        codon_allele_string 
        pep_allele_string 
        hgvs_genomic 
        hgvs_transcript 
        hgvs_protein 
        polyphen_prediction 
        polyphen_score 
        sift_prediction
        sift_score
    );
}

#sub _get_prediction_matrix {
#    my ($self, $analysis, $transcript_stable_id) = @_;
#    
#    # look in the protein function prediction table to see if there is 
#    # a prediction string for this transcript
#    
#    return undef unless ($analysis eq 'polyphen' || $analysis eq 'sift');
#     
#    my $dbh = $self->dbc->db_handle;
#    
#    my $col = $analysis.'_predictions';
#
#    my $sth = $dbh->prepare_cached(qq{
#        SELECT  $col
#        FROM    protein_function_predictions
#        WHERE   transcript_stable_id = ?
#    });
#    
#    $sth->execute($transcript_stable_id);
#    
#    my ($raw_matrix) = $sth->fetchrow_array;
#   
#    $sth->finish;
#
#    return undef unless $raw_matrix;
#
#    my $matrix = Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix->new(
#        -analysis   => $analysis,
#        -matrix     => $raw_matrix,
#    );
#    
#    return $matrix;
#}

1;
