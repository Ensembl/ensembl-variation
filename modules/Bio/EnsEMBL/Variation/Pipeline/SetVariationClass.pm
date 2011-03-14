package Bio::EnsEMBL::Variation::Pipeline::SetVariationClass;

use strict;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Hive::Process;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(SO_variation_class);

use base ('Bio::EnsEMBL::Hive::Process');

sub run {

    my $self = shift;

    my $reg_file = $self->param('ensembl_registry')
        or die "ensembl_registry is a required parameter";
    
    my $species = $self->param('species')
        or die "species is a required parameter";

    my $var_id_start = $self->param('variation_id_start')
        or die "variation_id_start is required";

    my $var_id_stop  = $self->param('variation_id_stop')
        or die "variation_id_stop is required";

    my $temp_var_table = $self->param('temp_var_table')
        or die "temp_var_table is required";
    
    my $temp_var_feat_table = $self->param('temp_var_feat_table')
        or die "temp_var_feat_table is required";

    my $reg = 'Bio::EnsEMBL::Registry';
    
    $reg->load_all($reg_file, 0, 1);
  
    my $var_dba = $reg->get_DBAdaptor($species, 'variation')
        or die "failed to get variation DBA for $species";

    my $aa = $var_dba->get_AttributeAdaptor;

    my $dbh = $var_dba->dbc->db_handle;

    # fetch the failed_descriptions to avoid a join

    my $fds_sth = $dbh->prepare(qq{
        SELECT  failed_description_id, description
        FROM    failed_description
    });

    $fds_sth->execute;

    my %fds;

    while (my ($fd_id, $desc) = $fds_sth->fetchrow_array) {
        $fds{$fd_id} = $desc;
    }

    my $all_sth = $dbh->prepare_cached(qq{
        SELECT  v.variation_id, vf.variation_feature_id, vf.allele_string, fv.failed_description_id
        FROM    (variation v LEFT JOIN variation_feature vf ON v.variation_id = vf.variation_id) 
                LEFT JOIN failed_variation fv ON v.variation_id = fv.variation_id
        WHERE   v.variation_id >= ?
        AND     v.variation_id <= ?
    });
   
    my $vf_insert_sth = $dbh->prepare(qq{
        INSERT IGNORE INTO $temp_var_feat_table (variation_feature_id, class_attrib_id)
        VALUES (?,?)
    });
    
    my $v_insert_sth = $dbh->prepare(qq{
        INSERT IGNORE INTO $temp_var_table (variation_id, class_attrib_id)
        VALUES (?,?)
    });
    
    $all_sth->execute($var_id_start, $var_id_stop);

    my @unmapped_v_ids;

    while (my ($v_id, $vf_id, $allele_string, $fd_id) = $all_sth->fetchrow_array) {

        unless ($vf_id) {
            # this variation doesn't have a corresponding variation_feature
            push @unmapped_v_ids, $v_id;
            next;
        }

        my $ref_correct = 1;
        
        # check to see if this variation_feature is known not to match the reference allele,
        # as this tells us if we can call insertions or deletions, or have to resort to indel

        if (defined $fd_id) {
            my $fail_reason = $fds{$fd_id};

            if ($fail_reason eq 'None of the variant alleles match the reference allele') {
                $ref_correct = 0;
            }
        }

        my $so_term = SO_variation_class($allele_string, $ref_correct);

        my $attrib_id = $aa->attrib_id_for_type_value('SO_term', $so_term);

        die "No attrib_id for $so_term" unless defined $attrib_id;

        $vf_insert_sth->execute($vf_id, $attrib_id);
        
        $v_insert_sth->execute($v_id, $attrib_id);
    }

    # now we need to fetch the alleles for any variations that are not mapped
    # and work out their class
    
    if (@unmapped_v_ids) {
        
        my $id_str = join ',', @unmapped_v_ids;

        my $unmapped_sth = $dbh->prepare(qq{
            SELECT  variation_id, allele
            FROM    allele
            WHERE   variation_id IN ($id_str)
            GROUP BY allele
        });

        $unmapped_sth->execute or die "Failed to fetch unmapped variation alleles for variation ids: $id_str";

        my $unmapped_alleles;

        while (my ($v_id, $allele) = $unmapped_sth->fetchrow_array) {
            push @{ $unmapped_alleles->{$v_id} ||= [] }, $allele; 
        }

        for my $v_id (keys %$unmapped_alleles) {

            my $allele_string = join '/', @{ $unmapped_alleles->{$v_id} };

            # we don't know what the reference is here

            my $ref_correct = 0;

            my $so_term = SO_variation_class($allele_string, $ref_correct);

            my $attrib_id = $aa->attrib_id_for_type_value('SO_term', $so_term);

            die "No attrib_id for $so_term" unless defined $attrib_id;

            $v_insert_sth->execute($v_id, $attrib_id);
        }
    }
}

1;

