
use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::VariationFeatureOverlapAdaptor;

use Bio::EnsEMBL::Variation::VariationFeatureOverlap;
use Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele;
use Bio::EnsEMBL::Variation::DBSQL::OverlapConsequenceAdaptor;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);

sub new {
    my $class = shift;
    
    my $self = $class->SUPER::new(@_);
    
    return $self;
}

sub store {
    my ($self, $vfo) = @_;
    
    my $dbh = $self->dbc->db_handle;
    
    my $vfo_sth = $dbh->prepare_cached(qq{
       INSERT INTO vf_overlap (variation_feature_id, feature_type_id, feature_stable_id)
       VALUES (?,?,?) 
    });
    
    $vfo_sth->execute(
        $vfo->variation_feature->dbID, 
        $vfo->feature_type_id,
        $vfo->feature->stable_id
    );
    
    my $vfo_db_id = $vfo_sth->{mysql_insertid};
    
    $vfo_sth->finish;
    
    $vfo->dbID($vfo_db_id);
    
    for my $vfoa (@{ $vfo->alleles }) {
        
        my $vfoa_sth = $dbh->prepare_cached(qq{
            INSERT INTO vf_overlap_allele (vf_overlap_id, allele, is_reference)
            VALUES (?,?,?)
        });
        
        $vfoa_sth->execute($vfo_db_id, $vfoa->seq, $vfoa->is_reference);
        
        my $vfoa_db_id = $vfoa_sth->{mysql_insertid};
        
        $vfoa_sth->finish;
        
        $vfoa->dbID($vfoa_db_id);
    }
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

sub _fetch_attrib_types {
    my ($self, $thing) = @_;
    
    my $dbh = $self->dbc->db_handle;
     
    my $attrib_str = join ',', map {"'$_'"} $self->_attrib_codes;
    
    my $sth = $dbh->prepare_cached(qq{
        SELECT  code, attrib_type_id 
        FROM    attrib_type 
        WHERE   code IN ( $attrib_str )
    });
    
    $sth->execute;
    
    while (my $row = $sth->fetchrow_hashref) {
        $self->{_attrib_types_by_code}->{$row->{code}} = $row->{attrib_type_id};
        $self->{_attrib_types_by_id}->{$row->{attrib_type_id}} = $row->{code};
    }
}

sub _attrib_id_for_code {
    my ($self, $code) = @_;
    
    $self->_fetch_attrib_types unless $self->{_attrib_types_by_code};
    
    return $self->{_attrib_types_by_code}->{$code};
}

sub _attrib_code_for_id {
    my ($self, $id) = @_;
    
    $self->_fetch_attrib_types unless $self->{_attrib_types_by_id};
    
    return $self->{_attrib_types_by_id}->{$id};
}

# fetch the associated  VariationFeatureOverlapAlleles for a list of 
# VariationFeatureOverlap objects fetched by one of the various fetch_all* methods
sub _fetch_alleles {
    
    my ($self, $vfos) = @_;
    
    my %vfos_by_id = map { $_->dbID => $_ } @$vfos;
    
    my $dbh = $self->dbc->db_handle;
    
    my $vfo_id_str = join ',', keys %vfos_by_id;
    
    my $allele_sth = $dbh->prepare(qq{
        SELECT  vf_overlap_id, vf_overlap_allele_id, allele, is_reference
        FROM    vf_overlap_allele 
        WHERE   vf_overlap_id IN ( $vfo_id_str )
    });
    
    $allele_sth->execute;
    
    my ($vfo_id, $vfoa_id, $seq, $is_ref);
    
    $allele_sth->bind_columns(\$vfo_id, \$vfoa_id, \$seq, \$is_ref);
    
    while ($allele_sth->fetch) {
        
        my $vfo = $vfos_by_id{$vfo_id};
        
        my $vfoa = Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele->new_fast({
            variation_feature_overlap   => $vfo,
            dbID                        => $vfoa_id,
            seq                         => $seq,
            is_reference                => $is_ref,
        });
        
        if ($is_ref) {
            $vfo->reference_allele($vfoa);
        }
        else {
            $vfo->alt_alleles($vfoa);
        }
    }
    
    # return the list of vfos for your method chaining convenience
    
    return $vfos;
}

sub fetch_all_by_VariationFeatures_with_constraint {
    
    my ($self, $vfs, $constraint) = @_;
    
    my %vfs_by_id;
    
    map { $vfs_by_id{$_->dbID} = $_ } @$vfs;
    
    my $id_str = join ',', map {"'$_'"} keys %vfs_by_id;
    
    my $full_constraint = "vfo.variation_feature_id in ( $id_str )";
    
    $full_constraint .= " AND $constraint" if $constraint;
    
    my $vfos = $self->generic_fetch($full_constraint);
    
    # we might as well resolve the vf_ids to the actual
    # objects now as we have them anyway (otherwise they would be
    # lazy-loaded)
    
    for my $vfo (@$vfos) {
        if ($vfo->{_vf_id}) {
            my $vf_id = delete $vfo->{_vf_id};
            $vfo->variation_feature($vfs_by_id{$vf_id});
        }
    }
    
    # fetch the associated VariationFeatureOverlapAlleles and return
    
    return $self->_fetch_alleles($vfos);
}

sub fetch_all_by_Features_with_constraint {
    
    my ($self, $features, $constraint) = @_;
    
    my %trans_by_id;
    
    map { $trans_by_id{$_->stable_id} = $_ } @$features;
    
    my $id_str = join ',', map {"'$_'"} keys %trans_by_id;
    
    my $full_constraint = "vfo.feature_stable_id in ( $id_str )";
    
    $full_constraint .= " AND $constraint" if $constraint;
    
    my $vfos = $self->generic_fetch($full_constraint);
    
    # we might as well resolve the feature_stable_ids to the actual 
    # objects now as we have them anyway (otherwise they would be
    # lazy-loaded)
    
    for my $vfo (@$vfos) {
        if ($vfo->{_feature_stable_id}) {
            my $tran_id = delete $vfo->{_feature_stable_id};
            $vfo->feature($trans_by_id{$tran_id});
        }
    }
    
    # fetch the associated VariationFeatureOverlapAlleles and return
    
    return $self->_fetch_alleles($vfos);
}

sub fetch_all_by_Features {
    my ($self, $features) = @_;
    
    my $constraint = 's.somatic = 0';
    
    return $self->fetch_all_by_Features_with_constraint($features, $constraint);
}

sub fetch_all_somatic_by_Features {
    my ($self, $features) = @_;
    
    my $constraint = 's.somatic = 1';
    
    return $self->fetch_all_by_Features_with_constraint($features, $constraint);
}

sub fetch_all_by_VariationFeatures {
    my ($self, $vfs) = @_;
    
    my $constraint = 's.somatic = 0';
    
    return $self->fetch_all_by_VariationFeatures_with_constraint($vfs, $constraint);
}

sub fetch_all_somatic_by_VariationFeatures {
    my ($self, $vfs) = @_;
    
    my $constraint = 's.somatic = 1';
    
    return $self->fetch_all_by_VariationFeatures_with_constraint($vfs, $constraint);
}

sub _objs_from_sth {
    my ($self, $sth) = @_;
    
    my ($vfo_id, $vf_id, $feat_stable_id, $feat_type_id);
    
    $sth->bind_columns(\$vfo_id, \$vf_id, \$feat_stable_id, \$feat_type_id);
    
    my @results;
    
    while ($sth->fetch) {
        my $vfo = Bio::EnsEMBL::Variation::VariationFeatureOverlap->new_fast({
           dbID                 => $vfo_id,
           _vf_id               => $vf_id,
           _feature_stable_id   => $feat_stable_id,
           feature_type_id      => $feat_type_id,
        });
        
        push @results, $vfo;
    }
    
    return \@results;
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
    return qw(
        vfo.vf_overlap_id
        vfo.variation_feature_id 
        vfo.feature_stable_id 
        vfo.feature_type_id
	);
}

1;
