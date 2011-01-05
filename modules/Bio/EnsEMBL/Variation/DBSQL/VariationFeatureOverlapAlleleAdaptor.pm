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

package Bio::EnsEMBL::Variation::DBSQL::VariationFeatureOverlapAlleleAdaptor;

use Bio::EnsEMBL::Variation::DBSQL::BaseVariationAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp expand);

use Bio::EnsEMBL::Variation::TranscriptVariationNew;

use base qw(Bio::EnsEMBL::Variation::DBSQL::BaseVariationAdaptor);

sub store {
    my ($self, $vfo) = @_;
    
    my $sth = $self->prepare(qq{
       INSERT INTO vf_overlap (variation_feature_id, feature_type_id, feature_stable_id)
       VALUES (?,?,?) 
    });
    
    $sth->execute(
        $vfo->variation_feature->dbID, 
        $vfo->feature_type_id,
        $vfo->feature->stable_id
    );
    
    $sth->finish;
    
    my $db_id = $sth->{mysql_insertid};
    
    $vfo->dbID($db_id);
}

sub _tables {
    return (
        ['vf_overlap','vfo'],
        ['variation_feature', 'vf'],
        ['source', 's']
    );
}

sub _default_where_clause {
  return 'vfo.variation_feature_id = vf.variation_feature_id AND vf.source_id = s.source_id';
}

sub _columns {
    return qw (tv.transcript_variation_id tv.transcript_stable_id
	       tv.variation_feature_id tv.cdna_start tv.cdna_end tv.cds_start tv.cds_end
	       tv.translation_start tv.translation_end
	       tv.peptide_allele_string tv.consequence_type
	       );
}
1;
