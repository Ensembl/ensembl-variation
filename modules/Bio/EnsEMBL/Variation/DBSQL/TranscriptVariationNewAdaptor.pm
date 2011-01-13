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

package Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationNewAdaptor;

use Bio::EnsEMBL::Variation::TranscriptVariationNew;
use Bio::EnsEMBL::Variation::TranscriptVariationAllele;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp expand);

use base qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);

sub store {
    my ($self, $tv) = @_;
    
    my $dbh = $self->dbc->db_handle;
    
    my $sth = $dbh->prepare_cached(q{
        INSERT INTO transcript_variation_allele (
            variation_feature_id,
            feature_stable_id,
            allele_string,
            is_somatic,
            consequence_types,
            cds_start,
            cds_end,
            cdna_start,
            cdna_end,
            pep_start,
            pep_end,
            codon_allele_string,
            peptide_allele_string,
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
            (join ",", map { $_->SO_id } @{ $allele->consequence_types }),
            $tv->cds_start, 
            $tv->cds_end,
            $tv->cdna_start,
            $tv->cdna_end,
            $tv->pep_start,
            $tv->pep_end,
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

# converts VariationFeatureOverlap objects into TranscriptVariation 
# objects, fetching any relevant attributes from the variation_feature_overlap_attrib 
# and variation_feature_overlap_allele_attrib tables
sub _vfos_to_tvs {
    
    my ($self, $vfos) = @_;
    
    my $dbh = $self->dbc->db_handle;
    
    my %vfos_by_id;
    my %vfoas_by_id;
    
    # rebless the VariationFeatureOverlap and VariationFeatureOverlapAllele objects 
    # into TranscriptVariation and TranscriptVariationAllele objects
    
    for my $vfo (@$vfos) {
        
        bless $vfo, 'Bio::EnsEMBL::Variation::TranscriptVariationNew';
        
        $vfos_by_id{$vfo->dbID} = $vfo;
        
        for my $vfoa (@{ $vfo->alleles }) {
            
            bless $vfoa, 'Bio::EnsEMBL::Variation::TranscriptVariationAllele';
            
            $vfoas_by_id{$vfoa->dbID} = $vfoa;
        }
    }

    # fetch the TranscriptVariation specific attributes
    
    if (my $vfo_id_str = join ',', keys %vfos_by_id) {
        
        my $tv_attrib_sth = $dbh->prepare(qq{
            SELECT  variation_feature_overlap_id, attrib_type_id, value 
            FROM    variation_feature_overlap_attrib 
            WHERE   variation_feature_overlap_id IN ( $vfo_id_str )
        });
        
        $tv_attrib_sth->execute;
        
        my ($vfo_id, $vfo_attrib_id, $vfo_value);
        
        $tv_attrib_sth->bind_columns(\$vfo_id, \$vfo_attrib_id, \$vfo_value);
        
        while ($tv_attrib_sth->fetch) {
            my $attrib = $self->_attrib_code_for_id($vfo_attrib_id);
            $vfos_by_id{$vfo_id}->$attrib($vfo_value);
        }
    }
    
    # and the TranscriptVariationAllele specific attributes
    
    if (my $vfoa_id_str = join ',', keys %vfoas_by_id) {
    
        my $tva_attrib_sth = $dbh->prepare(qq{
            SELECT  variation_feature_overlap_allele_id, attrib_type_id, value 
            FROM    variation_feature_overlap_allele_attrib 
            WHERE   variation_feature_overlap_allele_id IN ( $vfoa_id_str ) 
        });
        
        my ($vfoa_id, $vfoa_attrib_id, $vfoa_value);
        
        $tva_attrib_sth->execute;
        
        $tva_attrib_sth->bind_columns(\$vfoa_id, \$vfoa_attrib_id, \$vfoa_value);
        
        while ($tva_attrib_sth->fetch) {
            my $attrib = $self->_attrib_code_for_id($vfoa_attrib_id);
            
            my $allele = $vfoas_by_id{$vfoa_id};
                    
            if ($attrib eq 'codon') {
                $allele->codon($vfoa_value);
            }
            elsif ($attrib eq 'peptide') {
                $allele->peptide($vfoa_value);
            }
            elsif ($attrib eq 'SO_id') {
                $allele->consequences($self->_overlap_consequences->{$vfoa_value});
            }
            elsif ($attrib =~ m/hgvs_/) {
                $allele->$attrib($vfoa_value);
            }
            else {
                warn "Unexpected attribute: $attrib";
            }
        }
    }
    
    return $vfos;
}

sub fetch_by_dbID {
    my ($self, $dbID) = @_;
    
    my $constraint;
    
    if ($dbID =~ /,/) {
        my ($vf_id, $feat_stable_id) = split /,/, $dbID;
        $constraint = "variation_feature_id = $vf_id AND tva.feature_stable_id = $feat_stable_id";
    }
    else {
        $constraint = "tva.transcript_variation_allele_id = $dbID";
    }
    
    return $self->generic_fetch($constraint);
}

sub fetch_all_by_Transcripts {
    my ($self, $transcripts) = @_;
    return $self->fetch_all_by_Transcripts_with_constraint($transcripts,'tva.is_somatic = 0');
}

sub fetch_all_somatic_by_Transcripts {
    my ($self, $transcripts) = @_;
    return $self->fetch_all_by_Transcripts_with_constraint($transcripts,'tva.is_somatic = 1');
}

sub fetch_all_by_Transcripts_with_constraint {
    
    my ($self, $transcripts, $constraint) = @_;
    
    my $dbh = $self->dbc->db_handle;
   
    my %trans_by_id = map { $_->stable_id => $_ } @$transcripts;
    
    my $id_str = join',', map {"'".$_."'"} keys %trans_by_id;
    
    my $full_constraint = "feature_stable_id in ( $id_str )";
    $full_constraint .= " AND $constraint" if $constraint;
    
    my $tvs = $self->generic_fetch($full_constraint);
    
    for my $tv (@$tvs) {
        if ($tv->{_feature_stable_id}) {
            my $tran_id = delete $tv->{_feature_stable_id};
            $tv->{feature} = $trans_by_id{$tran_id};
        }
    }
    
    return $tvs;
}

sub fetch_all_by_VariationFeatures {
    
    my ($self, $transcripts) = @_;
    
    my $dbh = $self->dbc->db_handle;
   
    my %trans_by_id = map { $_->stable_id => $_ } @$transcripts;
    
    my $id_str = join',', map {"'".$_."'"} keys %trans_by_id;
    
    my $full_constraint = "feature_stable_id in ( $id_str )";
    
    my $tvs = $self->generic_fetch($full_constraint);
    
    for my $tv (@$tvs) {
        if ($tv->{_feature_stable_id}) {
            my $tran_id = delete $tv->{_feature_stable_id};
        }
    }
    
    return $tvs;
}

sub _objs_from_sth {
    my ($self, $sth) = @_;
    
    my (
        $transcript_variation_allele_id,
        $variation_feature_id, 
        $feature_stable_id, 
        $allele_string,
        $consequence_types,
        $cds_start,
        $cds_end,
        $cdna_start,
        $cdna_end,
        $pep_start,
        $pep_end,
        $codon_allele_string,
        $peptide_allele_string,
        $hgvs_genomic,
        $hgvs_coding,
        $hgvs_protein,
        $polyphen_prediction,
        $sift_prediction,
    );
    
    $sth->bind_columns(
        \$transcript_variation_allele_id,
        \$variation_feature_id, 
        \$feature_stable_id, 
        \$allele_string,
        \$consequence_types,
        \$cds_start,
        \$cds_end,
        \$cdna_start,
        \$cdna_end,
        \$pep_start,
        \$pep_end,
        \$codon_allele_string,
        \$peptide_allele_string,
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
        my ($ref_pep, $alt_pep)         = split /\//, $peptide_allele_string || '';
        
        # for HGMD mutations etc. just set the alt allele to the ref allele
        $alt_allele ||= $ref_allele;
        
        # for synonymous mutations the peptides are the same and 
        # there is no / in the string
        $alt_pep ||= $ref_pep;
        
        my $key = $variation_feature_id.'_'.$feature_stable_id;
        
        my $tv = $tvs{$key};
        
        unless ($tv) {
            $tv = Bio::EnsEMBL::Variation::TranscriptVariationNew->new_fast({
                dbID                    => $variation_feature_id.','.$feature_stable_id,
                _variation_feature_id   => $variation_feature_id,
                _feature_stable_id      => $feature_stable_id,
                cds_start               => $cds_start,
                cds_end                 => $cds_end,
                cdna_start              => $cdna_start,
                cdna_end                => $cdna_end,
                pep_start               => $pep_start,
                pep_end                 => $pep_end,
                adaptor                 => $self,
            });
            
            $tvs{$key} = $tv;
            
            my $ref_allele = Bio::EnsEMBL::Variation::TranscriptVariationAllele->new_fast({
                is_reference                => 1,
                variation_feature_allele    => $ref_allele,
                codon                       => $ref_codon,
                peptide                     => $ref_pep, 
            });
        }
        
        my @cons_types = map { $self->_overlap_consequences->{$_} } split /,/, $consequence_types; # / comment exists to satisfy eclipse!
        my $allele = Bio::EnsEMBL::Variation::TranscriptVariationAllele->new_fast({
            is_reference                => 0,
            variation_feature_allele    => $alt_allele,
            codon                       => $alt_codon,
            peptide                     => $alt_pep,
            hgvs_genomic                => $hgvs_genomic,
            hgvs_coding                 => $hgvs_coding,
            hgvs_protein                => $hgvs_protein,
            consequence_types           => \@cons_types, 
            polyphen_prediction         => $polyphen_prediction,
            sift_prediction             => $sift_prediction, 
        });
        
        $tv->alt_alleles($allele);
    }
    
    return [values %tvs];
}

sub _tables {
    return (
        ['transcript_variation_allele', 'tva']
    );
}

sub _columns {
    return qw(
        tva.transcript_variation_allele_id 
        tva.variation_feature_id 
        tva.feature_stable_id 
        tva.allele_string 
        tva.consequence_types 
        tva.cds_start 
        tva.cds_end 
        tva.cdna_start 
        tva.cdna_end 
        tva.pep_start 
        tva.pep_end 
        tva.codon_allele_string 
        tva.peptide_allele_string 
        tva.hgvs_genomic 
        tva.hgvs_coding 
        tva.hgvs_protein 
        tva.polyphen_prediction 
        tva.sift_prediction
    );
}

# returns a hashref of OverlapConsequences keyed by SO id
sub _overlap_consequences {
    my ($self) = @_;
    
    unless ($self->{_overlap_consequences}) {
        
        my $oca = $self->db->get_OverlapConsequenceAdaptor;

        $self->{_overlap_consequences} = { 
            map { $_->SO_id => $_ } @{ $oca->fetch_all } 
        };
    }
    
    return $self->{_overlap_consequences};
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
        FROM    ${program}_prediction pred, protein_position pp
        WHERE   pred.protein_position_id = pp.protein_position_id
        AND     pp.transcript_stable_id = ?
        AND     pp.transcript_version = ?
        AND     pp.position = ?
        AND     pred.amino_acid = ?
    });
    
    $sth->execute(
        $tva->transcript->stable_id,
        $tva->transcript->version,
        $tva->transcript_variation->pep_start,
        $tva->peptide,
    );
    
    my ($prediction) = $sth->fetchrow_array;
    
    return $prediction;
}

1;
