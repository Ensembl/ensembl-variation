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

sub store_normalised {
    my ($self, $tv) = @_;
    
    # must call the superclass first to set the dbID
    
    $self->SUPER::store($tv);
    
    my $dbh = $self->dbc->db_handle;
    
    if (defined $tv->cds_start) {
        
        my $vfo_attrib_sth = $dbh->prepare_cached(q{
           INSERT INTO variation_feature_overlap_attrib (variation_feature_overlap_id, attrib_type_id, value) 
           VALUES (?,?,?)
        });
        
        my @ids;
        my @types;
        my @values;
        
        for my $attrib (qw(cds_start cds_end cdna_start cdna_end pep_start pep_end)) {
            if (defined $tv->$attrib) {
                push @ids, $tv->dbID;
                push @types, $self->_attrib_id_for_code($attrib);
                push @values, $tv->$attrib;
            }
        }
        
        $vfo_attrib_sth->execute_array({}, \@ids, \@types, \@values);
    }
    
    for my $tva (@{ $tv->alt_alleles }) {
       
        my $vfoa_attrib_sth = $dbh->prepare_cached(q{
           INSERT INTO variation_feature_overlap_allele_attrib (variation_feature_overlap_allele_id, attrib_type_id, value) 
           VALUES (?,?,?)
        });
        
        my @ids;
        my @types;
        my @values;
        
        if ($tva->codon) { 
            for my $attrib (qw(codon peptide)) {
                push @ids, $tva->dbID;
                push @types, $self->_attrib_id_for_code($attrib);
                push @values, $tva->$attrib;
            }
        }
        
        for my $cons (@{ $tva->consequences }) {
            push @ids, $tva->dbID;
            push @types, $self->_attrib_id_for_code('SO_id');
            push @values, $cons->SO_id;
        }
     
        #ÊHGVS notation
        for my $ref (qw( genomic coding protein rna mitochondrial )) {
            my $attrib = qq{hgvs_$ref};
            if (my $val = $tva->$attrib) {
                push @ids, $tva->dbID;
                push @types, $self->_attrib_id_for_code($attrib);
                push @values, $val;
            }
        }
           
        $vfoa_attrib_sth->execute_array({}, \@ids, \@types, \@values);
    }
}

sub _get_nsSNP_prediction {
    my ($self, $program, $tva) = @_;
    
    return undef unless ($program eq 'polyphen' || $program eq 'sift');
    
    return undef unless ($tva->peptide && length($tva->peptide) == 1);
    
    my $dbh = $self->dbc->db_handle;
    
    my $sth = $dbh->prepare_cached(qq{
        SELECT  p.prediction
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
        $tva->pep_start,
        $tva->peptide,
    );
    
    my ($prediction) = $sth->fetchrow_array;
    
    return $prediction;
}

sub _attrib_codes {
    return qw(
        cds_start
        cds_end
        cdna_start
        cdna_end
        pep_start
        pep_end
        codon
        peptide
        SO_id
        hgvs_genomic
        hgvs_coding
        hgvs_protein
        hgvs_rna
        hgvs_mitochondrial
    );
}

sub fetch_all_by_Transcripts2 {
    
    # another version of fetch_all_by_Transcripts that tries to get all
    # data needed in a single SQL query, but it's currently slower than
    # the simpler way calling the superclass and fetching the alleles
    # separately
    
    my ($self, $transcripts) = @_;
      
    my $dbh = $self->dbc->db_handle;
    
    my $id_str = join ',', map { "'".$_->stable_id."'" } @$transcripts;
    
    my $sth = $dbh->prepare(qq{
        SELECT vfo.variation_feature_overlap_id, vfo.variation_feature_id, vfo.feature_stable_id, 
            vfo.feature_type_id, vfo_att.code, vfoat.value, vfoa.variation_feature_overlap_allele_id, 
            vfoa.allele, vfoa_att.code, vfoaa.value
        FROM variation_feature_overlap vfo 
        LEFT JOIN variation_feature_overlap_allele vfoa         ON vfo.variation_feature_overlap_id = vfoa.variation_feature_overlap_id 
        LEFT JOIN variation_feature_overlap_attrib vfoat        ON vfo.variation_feature_overlap_id = vfoat.variation_feature_overlap_id 
        LEFT JOIN variation_feature_overlap_allele_attrib vfoaa ON vfoa.variation_feature_overlap_allele_id = vfoaa.variation_feature_overlap_allele_id
        LEFT JOIN attrib_type vfo_att            ON vfo_att.attrib_type_id = vfoat.attrib_type_id 
        LEFT JOIN attrib_type vfoa_att           ON vfoa_att.attrib_type_id = vfoaa.attrib_type_id
        JOIN variation_feature vf                ON vf.variation_feature_id = vfo.variation_feature_id
        JOIN source s                            ON s.source_id = vf.source_id
        WHERE vfo.feature_stable_id IN ( $id_str )
        AND s.somatic = 0
    });
    
    $sth->execute;
    
    my ($vfo_id, $vf_id, $tran_stable_id, $feat_type_id, $vfo_attrib, $vfo_value, 
        $vfoa_id, $allele, $vfoa_attrib, $vfoa_value);
    
    $sth->bind_columns(\$vfo_id, \$vf_id, \$tran_stable_id, \$feat_type_id, \$vfo_attrib, \$vfo_value, 
        \$vfoa_id, \$allele, \$vfoa_attrib, \$vfoa_value);
    
    my %vfos_by_id;
    my %vfoas_by_id;
    
    while ($sth->fetch) {
        
        my $vfo = $vfos_by_id{$vfo_id};
        
        unless ($vfo) {
            $vfo = Bio::EnsEMBL::Variation::TranscriptVariationNew->new_fast({
                dbID                    => $vfo_id,
                _variation_feature_id   => $vf_id,
                _feature_stable_id      => $tran_stable_id,
                feature_type_id         => $feat_type_id,
            });
            
            $vfos_by_id{$vfo_id} = $vfo
        }
        
        $vfo->$vfo_attrib($vfo_value) if $vfo_attrib;
        
        if ($vfoa_id) {
            my $vfoa = $vfoas_by_id{$vfoa_id};
            
            unless ($vfoa) {
                 $vfoa = Bio::EnsEMBL::Variation::TranscriptVariationAllele->new_fast({
                    variation_feature_overlap   => $vfo,
                    dbID                        => $vfoa_id,
                    seq                         => $allele,
                });
                
                $vfo->alt_alleles($vfoa);
                
                $vfoas_by_id{$vfoa_id} = $vfoa;
            }
            
            if ($vfoa_attrib) {
                if ($vfoa_attrib =~ /codon|peptide/) {
                    $vfoa->$vfoa_attrib($vfoa_value);
                }
                elsif ($vfoa_attrib eq 'SO_id') {
                    $vfoa->consequences($self->_overlap_consequences->{$vfoa_value});
                }
                elsif ($vfoa_attrib =~ m/hgvs_/) {
                    $allele->$vfoa_attrib($vfoa_value);
                }
                else {
                    warn "Unexpected attribute: $vfoa_attrib";
                }
            }
        }
    }
    
    return [values %vfos_by_id];
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

sub fetch_all_by_VariationFeatures_normalised {  
    my ($self, $vfs) = @_;
    return $self->_vfos_to_tvs($self->SUPER::fetch_all_by_VariationFeatures($vfs));
}

sub fetch_all_somatic_by_VariationFeatures_normalised {   
    my ($self, $vfs) = @_;
    return $self->_vfos_to_tvs($self->SUPER::fetch_all_somatic_by_VariationFeatures($vfs));
}

sub fetch_all_by_Transcripts_normalised {   
    my ($self, $transcripts) = @_;
    return $self->_vfos_to_tvs($self->SUPER::fetch_all_by_Features($transcripts));
}

sub fetch_all_somatic_by_Transcripts_normalised {
    my ($self, $transcripts) = @_;
    return $self->_vfos_to_tvs($self->SUPER::fetch_all_somatic_by_Features($transcripts));
}

sub fetch_by_dbID_normalised {
    my $self = shift;
    return $self->_vfos_to_tvs([$self->SUPER::fetch_by_dbID(@_)])->[0];
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
        
        my @cons_types = map { $self->_overlap_consequences->{$_} } split /,/, $consequence_types;
        
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

1;
