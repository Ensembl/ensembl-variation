=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut
package Bio::EnsEMBL::Variation::Pipeline::SetVariationClass;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

use Bio::EnsEMBL::Variation::Utils::Sequence qw(SO_variation_class);
use Bio::EnsEMBL::Variation::Utils::Constants qw(:SO_class_terms);


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

        my $so_term = $self->assign_SO_variation_class($allele_string, $ref_correct);

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
            SELECT  DISTINCT a.variation_id, ac.allele
            FROM    allele a, allele_code ac
            WHERE   a.variation_id IN ($id_str)
            AND     a.allele_code_id = ac.allele_code_id
            ORDER BY a.variation_id, ac.allele
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

            my $so_term = $self->assign_SO_variation_class($allele_string, $ref_correct);

            my $attrib_id = $aa->attrib_id_for_type_value('SO_term', $so_term);

            die "No attrib_id for $so_term" unless defined $attrib_id;

            $v_insert_sth->execute($attrib_id, $v_id);
        }
    }
    
    $vf_insert_sth->finish();
    $v_insert_sth->finish();
}

sub assign_SO_variation_class {
    my $self = shift;
    my $allele_string = shift;
    my $ref_correct = shift;
    my $identify_marker_e = $self->param('identify_marker_e');
#    my $identify_marker_eg = $self->param('identify_marker_eg');
    my $so_term = SO_variation_class($allele_string, $ref_correct);

    my $sequence_alteration = SO_TERM_SEQUENCE_ALTERATION;
    my $genetic_marker = SO_TERM_GENETIC_MARKER;

    if (($so_term eq $sequence_alteration) && $identify_marker_e) {
        if ($allele_string =~ /^\([^\(\)\/]+\)$/) {
            # (D19S912) 
            $allele_string =~ s/\(|\)//g;
            my $cdba = $self->get_species_adaptor('core');
            my $ma = $cdba->get_MarkerAdaptor;
            my $marker_list = $ma->fetch_all_by_synonym($allele_string);
            if (scalar @$marker_list > 0) {
                $so_term = $genetic_marker;
            }
        } else {
            return $so_term;
        }
    }
    return $so_term;
}


1;

