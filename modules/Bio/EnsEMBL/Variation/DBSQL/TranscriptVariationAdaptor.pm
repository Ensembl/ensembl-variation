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


use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor;

use Bio::EnsEMBL::Variation::TranscriptVariation;
use Bio::EnsEMBL::Variation::TranscriptVariationAllele;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(Bio::EnsEMBL::Variation::DBSQL::VariationFeatureOverlapAdaptor);

sub store {
    my ($self, $tv) = @_;
    
    my $dbh = $self->dbc->db_handle;
    
    my $sth = $dbh->prepare_cached(q{
        INSERT INTO transcript_variation (
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
            codon_allele_string,
            pep_allele_string,
            hgvs_genomic,
            hgvs_coding,
            hgvs_protein,
            polyphen_prediction,
            sift_prediction
        ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
    });

    for my $allele (@{ $tv->alt_alleles }) {
        
        $sth->execute(
            $tv->variation_feature->dbID,
            $tv->feature->stable_id,
            $allele->allele_string,
            $tv->variation_feature->is_somatic,
            (join ',', map { $_->SO_term } @{ $allele->consequence_types }),
            $tv->cds_start, 
            $tv->cds_end,
            $tv->cdna_start,
            $tv->cdna_end,
            $tv->translation_start,
            $tv->translation_end,
            $allele->codon_allele_string,
            $allele->pep_allele_string,
            $allele->hgvs_genomic,
            $allele->hgvs_coding,
            $allele->hgvs_protein,
            $allele->polyphen_prediction,
            $allele->sift_prediction
        );
    }
}

sub fetch_all_by_Transcripts {
    my $self = shift;
    return $self->SUPER::fetch_all_by_Features(@_);
}

sub fetch_all_somatic_by_Transcripts {
    my $self = shift;
    return $self->SUPER::fetch_all_somatic_by_Features(@_);
}

sub fetch_all_by_Transcripts_with_constraint {
    my $self = shift;
    return $self->SUPER::fetch_all_by_Features_with_constraint(@_);
}

sub _objs_from_sth {
    my ($self, $sth) = @_;
    
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
        $codon_allele_string,
        $pep_allele_string,
        $hgvs_genomic,
        $hgvs_coding,
        $hgvs_protein,
        $polyphen_prediction,
        $sift_prediction,
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
        \$codon_allele_string,
        \$pep_allele_string,
        \$hgvs_genomic,
        \$hgvs_coding,
        \$hgvs_protein,
        \$polyphen_prediction,
        \$sift_prediction,
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
                adaptor                 => $self,
            });
            
            $tvs{$key} = $tv;
            
            my $ref_allele = Bio::EnsEMBL::Variation::TranscriptVariationAllele->new_fast({
                is_reference                => 1,
                variation_feature_seq       => $ref_allele,
                transcript_variation        => $tv, 
                codon                       => $ref_codon,
                peptide                     => $ref_pep, 
            });

            $tv->reference_allele($ref_allele);
        }
       
        #my $cons_types = $self->_transcript_variation_consequences_for_set_number($consequence_types);

        my $cons_types = [ map { $self->_overlap_consequence_for_SO_term($_) } 
            split /,/, $consequence_types ]; # / comment exists to satisfy eclipse!
        
        my $allele = Bio::EnsEMBL::Variation::TranscriptVariationAllele->new_fast({
            is_reference                => 0,
            variation_feature_seq       => $alt_allele,
            transcript_variation        => $tv, 
            codon                       => $alt_codon,
            peptide                     => $alt_pep,
            hgvs_genomic                => $hgvs_genomic,
            hgvs_coding                 => $hgvs_coding,
            hgvs_protein                => $hgvs_protein,
            consequence_types           => $cons_types, 
            polyphen_prediction         => $polyphen_prediction,
            sift_prediction             => $sift_prediction, 
        });
        
        $tv->alt_alleles($allele);
    }
    
    return [values %tvs];
}

sub _tables {
    return (
        ['transcript_variation']
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
        codon_allele_string 
        pep_allele_string 
        hgvs_genomic 
        hgvs_coding 
        hgvs_protein 
        polyphen_prediction 
        sift_prediction
    );
}

# looks in the (sift|polyphen)_prediction table to see if there is 
# a prediction for the given non-synonymous TranscriptVariationAllele
sub _get_nsSNP_prediction {
    my ($self, $program, $tva) = @_;
    
    return undef unless ($program eq 'polyphen' || $program eq 'sift');
    
    # we can only deal with single amino acid substitutions
    
    return undef unless ($tva->peptide && length($tva->peptide) == 1);
    
    my $dbh = $self->dbc->db_handle;
    
    my $sth = $dbh->prepare_cached(qq{
        SELECT  pred.prediction
        FROM    ${program}_prediction pred, protein_position pp, protein_info pi
        WHERE   pred.protein_position_id = pp.protein_position_id
        AND     pp.protein_info_id = pi.protein_info_id
        AND     pi.transcript_stable_id = ?
        AND     pp.position = ?
        AND     pred.amino_acid = ?
   });
    
    $sth->execute(
        $tva->transcript->stable_id,
        $tva->transcript_variation->translation_start,
        $tva->peptide,
    );
    
    my ($prediction) = $sth->fetchrow_array;
   
    $sth->finish;

    return $prediction;
}

1;
