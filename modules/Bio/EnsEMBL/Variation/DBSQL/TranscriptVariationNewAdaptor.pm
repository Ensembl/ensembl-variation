
use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationNewAdaptor;

use Bio::EnsEMBL::Variation::TranscriptVariationNew;
use Bio::EnsEMBL::Variation::TranscriptVariationAllele;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureOverlapAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp expand);

use base qw(Bio::EnsEMBL::Variation::DBSQL::VariationFeatureOverlapAdaptor);

our $AUTOLOAD;

sub store {
    my ($self, $tv) = @_;
    
    # must call the superclass first to set the dbID
    
    $self->SUPER::store($tv);
    
    my $dbh = $self->dbc->db_handle;
    
    if ($tv->cds_start) {
        
        my $vfo_attrib_sth = $dbh->prepare_cached(q{
           INSERT INTO vf_overlap_attrib (vf_overlap_id, attrib_type_id, value) 
           VALUES (?,?,?)
        });
        
        my @ids;
        my @types;
        my @values;
        
        for my $attrib (qw(cds_start cds_end cdna_start cdna_end pep_start pep_end)) {
            push @ids, $tv->dbID;
            push @types, $self->_attrib_id_for_code($attrib);
            push @values, $tv->$attrib;
        }
        
        $vfo_attrib_sth->execute_array({}, \@ids, \@types, \@values);
    }
    
    for my $tva (@{ $tv->alleles }) {
       
        my $vfoa_attrib_sth = $dbh->prepare_cached(q{
           INSERT INTO vf_overlap_allele_attrib (vf_overlap_allele_id, attrib_type_id, value) 
           VALUES (?,?,?)
        });
        
        my @ids;
        my @types;
        my @values;
        
        if ($tva->codon) { 
            for my $attrib (qw(codon amino_acid)) {
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
        
        $vfoa_attrib_sth->execute_array({}, \@ids, \@types, \@values);
    }
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
        amino_acid
        SO_id
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
        SELECT vfo.vf_overlap_id, vfo.variation_feature_id, vfo.feature_stable_id, 
            vfo.feature_type_id, vfo_att.code, vfoat.value, vfoa.vf_overlap_allele_id, 
            vfoa.allele, vfoa_att.code, vfoaa.value
        FROM vf_overlap vfo 
        LEFT JOIN vf_overlap_allele vfoa         ON vfo.vf_overlap_id = vfoa.vf_overlap_id 
        LEFT JOIN vf_overlap_attrib vfoat        ON vfo.vf_overlap_id = vfoat.vf_overlap_id 
        LEFT JOIN vf_overlap_allele_attrib vfoaa ON vfoa.vf_overlap_allele_id = vfoaa.vf_overlap_allele_id
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
                dbID                => $vfo_id,
                _vf_id              => $vf_id,
                _feature_stable_id  => $tran_stable_id,
                feature_type_id     => $feat_type_id,
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
                if ($vfoa_attrib =~ /codon|amino_acid/) {
                    $vfoa->$vfoa_attrib($vfoa_value);
                }
                elsif ($vfoa_attrib eq 'SO_id') {
                    $vfoa->consequences($self->_overlap_consequences->{$vfoa_value});
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
# objects, fetching any relevant attributes from the vf_overlap_attrib 
# and vf_overlap_allele_attrib tables
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
        
        for my $vfoa (@{ $vfo->alt_alleles }) {
            
            bless $vfoa, 'Bio::EnsEMBL::Variation::TranscriptVariationAllele';
            
            $vfoas_by_id{$vfoa->dbID} = $vfoa;
        }
    }

    # fetch the TranscriptVariation specific attributes
    
    if (my $vfo_id_str = join ',', keys %vfos_by_id) {
        
        my $tv_attrib_sth = $dbh->prepare(qq{
            SELECT  vf_overlap_id, attrib_type_id, value 
            FROM    vf_overlap_attrib 
            WHERE   vf_overlap_id IN ( $vfo_id_str )
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
            SELECT  vf_overlap_allele_id, attrib_type_id, value 
            FROM    vf_overlap_allele_attrib 
            WHERE   vf_overlap_allele_id IN ( $vfoa_id_str ) 
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
            elsif ($attrib eq 'amino_acid') {
                $allele->amino_acid($vfoa_value);
            }
            elsif ($attrib eq 'SO_id') {
                $allele->consequences($self->_overlap_consequences->{$vfoa_value});
            }
            else {
                warn "Unexpected attribute: $attrib";
            }
        }
    }
    
    return $vfos;
}

sub fetch_all_by_VariationFeatures {  
    my ($self, $vfs) = @_;
    return $self->_vfos_to_tvs($self->SUPER::fetch_all_by_VariationFeatures($vfs));
}

sub fetch_all_somatic_by_VariationFeatures {   
    my ($self, $vfs) = @_;
    return $self->_vfos_to_tvs($self->SUPER::fetch_all_somatic_by_VariationFeatures($vfs));
}

sub fetch_all_by_Transcripts {   
    my ($self, $transcripts) = @_;
    return $self->_vfos_to_tvs($self->SUPER::fetch_all_by_Features($transcripts));
}

sub fetch_all_somatic_by_Transcripts {
    my ($self, $transcripts) = @_;
    return $self->_vfos_to_tvs($self->SUPER::fetch_all_somatic_by_Features($transcripts));
}

sub fetch_by_dbID {
    my $self = shift;
    return $self->_vfos_to_tvs($self->SUPER::fetch_by_dbID(@_));
}

sub AUTOLOAD {
    my $self = shift;
    my $method = $AUTOLOAD;
    $method =~ s/.*://;
    
    if ($self->SUPER::can($method)) {
        $method = 'SUPER::'.$method;
        return $self->_vfos_to_tvs($self->$method(@_));
    }
}

1;
