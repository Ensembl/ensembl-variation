=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::Config;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess');

use FileHandle;
use JSON;

sub fetch_input {
	my $self = shift;
}

sub run {
	my $self = shift;
    $self->write_config_file();
}

sub write_output {
	my $self = shift;
    $self->dataflow_output_id({'config_file' => $self->param('config_file')}, 1);
    return 1;
}

sub write_config_file {
    my $self = shift;
    # structural_variation svs
    # somatic
    # incl consequences: protein info: sift, polyphen
    # evidence, clinical_significance, ancestral_allele, minor_allele_freq, validation_status
    # populations:
    # individuals:
    # sets: phenotypes, clinically_associated
    my $config = {};
    $self->variation_data_survey($config);

    my $populations = $self->param('populations');
    my $individuals = $self->param('individuals');

    my $species_variation_data = $self->get_all_species();

    foreach my $species (keys %$species_variation_data) {
        my $species_config = {
            failed => ['failed'],
            generic => ['evidence', 'validation_status'],
            incl_consequences => ['incl_consequences', 'protein_coding_details', 'evidence'],
        };
        foreach my $attribute (qw/ancestral_allele global_maf clinical_significance/) {
            if ($config->{$species}->{$attribute}) {
                push @{$species_config->{generic}}, $attribute;
                push @{$species_config->{incl_consequences}}, $attribute;
            }
        } 
        
        if ($config->{$species}->{sift}) {
            push @{$species_config->{incl_consequences}}, 'sift';
        }
        if ($config->{$species}->{svs}) {
            $species_config->{structural_variations} = ['structural_variations'];
            if ($config->{$species}->{clinical_significance_svs}) {
                push @{$species_config->{structural_variations}}, 'clinical_significance';
            }
        }
#        if ($config->{$species}->{phenotypes}) {
#            $species_config->{sets}->{phenotype_associated} = $species_config->{generic};
#        }
        if (defined $populations->{$species}) {
            $species_config->{populations} = $populations->{$species};
        }
        if (defined $individuals->{$species}) {
            $species_config->{individuals} = $individuals->{$species};
        }
        if ($species eq 'Homo_sapiens') {
            $species_config->{sets}->{clinically_associated} = ['evidence', 'ancestral_allele', 'clinical_significance', 'global_maf'];
            $species_config->{sets}->{phenotype_associated} =  ['evidence', 'ancestral_allele', 'clinical_significance', 'global_maf'];
            $species_config->{incl_consequences} =  ['sift', 'polyphen', 'incl_consequences', 'protein_coding_details', 'evidence', 'ancestral_allele', 'clinical_significance', 'global_maf'];
            $species_config->{somatic_incl_consequences} =  ['somatic', 'sift', 'polyphen', 'incl_consequences', 'protein_coding_details', 'evidence', 'ancestral_allele', 'clinical_significance', 'global_maf'];
            $species_config->{somatic} = ['somatic', 'evidence', 'ancestral_allele', 'clinical_significance', 'global_maf'];
            $species_config->{generic} = ['evidence', 'ancestral_allele', 'clinical_significance', 'global_maf', 'variation_id', 'allele_string'];
        }
        $config->{$species} = $species_config;
    }
    my $pipeline_dir = $self->param('pipeline_dir');
    my $config_file = "$pipeline_dir/data_dumps_config_human.json"; 
    my $fh = FileHandle->new($config_file, 'w');
    my $json = JSON->new->allow_nonref;
    print $fh $json->encode($config);
    $fh->close();
    $self->param('config_file', $config_file);
}

sub variation_data_survey {
    my $self = shift;
    my $config = shift;
    my $registry = 'Bio::EnsEMBL::Registry';
    $registry->load_all($self->param('registry_file'));
    my $vdbas = $registry->get_all_DBAdaptors(-group => 'variation');

    my $queries = {
        sift => 'select count(*) from protein_function_predictions;',
        ancestral_allele => 'select variation_id from variation where ancestral_allele is not null limit 1',
        global_maf => 'select variation_id from variation where minor_allele is not null limit 1',
        clinical_significance => 'select variation_id from variation where clinical_significance is not null limit 1',
        clinical_significance_svs => 'select structural_variation_id from structural_variation where clinical_significance is not null limit 1',
#        phenotypes => 'select count(*) from phenotype_feature',
        svs => 'select count(*) from structural_variation',
    };

    foreach my $data_type (keys %$queries) {
        my $sub_set_species = query_database($vdbas, $queries->{$data_type});
        foreach my $species (keys %$sub_set_species) {
            $config->{$species}->{$data_type} = 1;
        }
    }
}

sub query_database {
    my $vdbas = shift;
    my $query = shift;
    my $species_names = {};
    foreach my $vdba (@$vdbas) {
        my $species_name = $vdba->species();
        my $dbh = $vdba->dbc->db_handle;
        my $sth = $dbh->prepare($query);
        $sth->execute();
        while (my @row = $sth->fetchrow_array) {
            my $count = $row[0];
            if ($count > 0) {
                $species_names->{$species_name} = 1;
            }
        }
        $sth->finish();
    }
    return $species_names;
}

1;
