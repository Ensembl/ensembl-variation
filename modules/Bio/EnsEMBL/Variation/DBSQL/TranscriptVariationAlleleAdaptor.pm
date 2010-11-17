
use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationNewAdaptor;

use Bio::EnsEMBL::Variation::DBSQL::BaseVariationAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp expand);

use Bio::EnsEMBL::Variation::TranscriptVariationNew;

use base qw(Bio::EnsEMBL::Variation::DBSQL::VariationFeatureOverlapAdaptor);

sub store {
    my ($self, $tv) = @_;
    
    # must call the superclass first to set the dbID
    
    $self->SUPER::store($tv);
    
    if ($tv->feature->biotype eq 'protein_coding') {
        
        my $vf_overlap_id = $tv->dbID;
        
        my $sth = $self->prepare(qq{
           INSERT INTO vf_overlap_attrib (vf_overlap_id, attrib_type_id, value) 
           VALUES 
           ( $vf_overlap_id, (SELECT attrib_type_id FROM attrib_type WHERE code = 'cds_start'), ? ),
           ( $vf_overlap_id, (SELECT attrib_type_id FROM attrib_type WHERE code = 'cds_end'), ? ),
           ( $vf_overlap_id, (SELECT attrib_type_id FROM attrib_type WHERE code = 'cdna_start'), ? ),
           ( $vf_overlap_id, (SELECT attrib_type_id FROM attrib_type WHERE code = 'cdna_end'), ? ),
           ( $vf_overlap_id, (SELECT attrib_type_id FROM attrib_type WHERE code = 'translation_start'), ? ),
           ( $vf_overlap_id, (SELECT attrib_type_id FROM attrib_type WHERE code = 'translation_end'), ? ),
        });
        
        $sth->execute(
            $tv->cds_start,
            $tv->cds_end,
            $tv->cdna_start,
            $tv->cdna_end,
            $tv->translation_start,
            $tv->translation_end,
        );
        
    }
}

sub _tables {
    return (
        ['transcript_variation','tv'],
        ['variation_feature', 'vf'],
        ['source', 's']
    );
}

sub _default_where_clause {
  return 'tv.variation_feature_id = vf.variation_feature_id AND vf.source_id = s.source_id';
}

sub _columns {
    return qw (tv.transcript_variation_id tv.transcript_stable_id
	       tv.variation_feature_id tv.cdna_start tv.cdna_end tv.cds_start tv.cds_end
	       tv.translation_start tv.translation_end
	       tv.peptide_allele_string tv.consequence_type
	       );
}
1;
