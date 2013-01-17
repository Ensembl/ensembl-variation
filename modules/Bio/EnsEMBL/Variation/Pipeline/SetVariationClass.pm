=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/legal/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut
package Bio::EnsEMBL::Variation::Pipeline::SetVariationClass;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

use Bio::EnsEMBL::Variation::Utils::Sequence qw(SO_variation_class);


sub run {

    my $self = shift;
    
    my $var_id_start = $self->required_param('variation_id_start');

    my $var_id_stop  = $self->required_param('variation_id_stop');

    my $temp_var_table = $self->param('temp_var_table');
    
    my $temp_var_feat_table = $self->param('temp_var_feat_table');
  
    my $var_dba = $self->get_species_adaptor('variation');

    my $aa = $var_dba->get_AttributeAdaptor;

    my $dbc = $var_dba->dbc;

    # fetch the failed_descriptions to avoid a join

    my $fds_sth = $dbc->prepare(qq{
        SELECT  failed_description_id, description
        FROM    failed_description
    });

    $fds_sth->execute;

    my %fds;

    while (my ($fd_id, $desc) = $fds_sth->fetchrow_array) {
        $fds{$fd_id} = $desc;
    }
    
    $fds_sth->finish();

    my $all_sth = $dbc->prepare(qq{
        SELECT  v.variation_id, vf.variation_feature_id, vf.allele_string, fv.failed_description_id
        FROM    (variation v LEFT JOIN variation_feature vf ON v.variation_id = vf.variation_id) 
                LEFT JOIN failed_variation fv ON v.variation_id = fv.variation_id, source s
        WHERE   v.variation_id >= ?
        AND     v.variation_id <= ?
        AND     v.source_id = s.source_id
        AND     s.name != 'HGMD-PUBLIC'
    });

    my $vf_insert_sth;
    my $v_insert_sth;

    if (defined $temp_var_feat_table) {
        $vf_insert_sth = $dbc->prepare(qq{
            INSERT IGNORE INTO $temp_var_feat_table (class_attrib_id, variation_feature_id)
            VALUES (?,?)
        });
    }
    else {
        $vf_insert_sth = $dbc->prepare(qq{
            UPDATE variation_feature SET class_attrib_id = ? WHERE variation_feature_id = ?
        });
    }

    if (defined $temp_var_table) {
        $v_insert_sth = $dbc->prepare(qq{
            INSERT IGNORE INTO $temp_var_table (class_attrib_id, variation_id)
            VALUES (?,?)
        });
    }
    else {
        $v_insert_sth = $dbc->prepare(qq{
            UPDATE variation SET class_attrib_id = ? WHERE variation_id = ?
        });
    }
    
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

        $vf_insert_sth->execute($attrib_id, $vf_id);
        
        $v_insert_sth->execute($attrib_id, $v_id);
    }
    
    $all_sth->finish();

    # now we need to fetch the alleles for any variations that are not mapped
    # and work out their class
    
    if (@unmapped_v_ids) {
        
        my $id_str = join ',', @unmapped_v_ids;

        my $unmapped_sth = $dbc->prepare(qq{
            SELECT  a.variation_id, ac.allele
            FROM    allele a, allele_code ac
            WHERE   a.variation_id IN ($id_str)
            AND     a.allele_code_id = ac.allele_code_id
            GROUP BY ac.allele
        });

        $unmapped_sth->execute or die "Failed to fetch unmapped variation alleles for variation ids: $id_str";

        my $unmapped_alleles;

        while (my ($v_id, $allele) = $unmapped_sth->fetchrow_array) {
            push @{ $unmapped_alleles->{$v_id} ||= [] }, $allele; 
        }
        
        $unmapped_sth->finish();

        for my $v_id (keys %$unmapped_alleles) {

            my $allele_string = join '/', @{ $unmapped_alleles->{$v_id} };

            # we don't know what the reference is here

            my $ref_correct = 0;

            my $so_term = SO_variation_class($allele_string, $ref_correct);

            my $attrib_id = $aa->attrib_id_for_type_value('SO_term', $so_term);

            die "No attrib_id for $so_term" unless defined $attrib_id;

            $v_insert_sth->execute($attrib_id, $v_id);
        }
    }
    
    $vf_insert_sth->finish();
    $v_insert_sth->finish();
}

1;

