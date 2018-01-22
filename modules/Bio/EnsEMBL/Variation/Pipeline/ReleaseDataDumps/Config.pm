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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::Config;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess');
use Bio::EnsEMBL::Registry;
use FileHandle;
use JSON;

sub run {
	my $self = shift;
    $self->write_config_file();
}


sub write_config_file {
    my $self = shift;
    my $species = $self->param('species');
    my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');
    my $vdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'variation');
    my @input;
    my $params = {};
    my $division = $self->division($dba);


    # structural_variation svs
    # somatic
    # incl consequences: protein info: sift, polyphen
    # evidence, clinical_significance, ancestral_allele, minor_allele_freq, validation_status
    # populations:
    # individuals:
    # sets: phenotypes, clinically_associated
    my $config = {};

    $self->variation_data_survey($config,$species,$vdba);

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
    if ($species eq 'Homo_sapiens') {
        $species_config->{sets}->{clinically_associated} = ['evidence', 'ancestral_allele', 'clinical_significance', 'global_maf'];
        $species_config->{sets}->{phenotype_associated} =  ['evidence', 'ancestral_allele', 'clinical_significance', 'global_maf'];
        $species_config->{incl_consequences} =  ['sift', 'polyphen', 'incl_consequences', 'protein_coding_details', 'evidence', 'ancestral_allele', 'clinical_significance', 'global_maf'];
        $species_config->{somatic_incl_consequences} =  ['somatic', 'sift', 'polyphen', 'incl_consequences', 'protein_coding_details', 'evidence', 'ancestral_allele', 'clinical_significance', 'global_maf'];
        $species_config->{somatic} = ['somatic', 'evidence', 'ancestral_allele', 'clinical_significance', 'global_maf'];
        $species_config->{generic} = ['evidence', 'ancestral_allele', 'clinical_significance', 'global_maf', 'variation_id', 'allele_string'];
    }
    $config->{$species} = $species_config;
    
    $params->{species} = $species;
    $params->{config}  = $config->{$species};
    if ($division ne '') {
      $params->{species_division} = $division;
    }
    push @input, $params;
    
    $self->dataflow_output_id(\@input, 2);
    $self->dataflow_output_id(\@input, 1);
}

sub variation_data_survey {
    my ($self,$config,$species,$vdba)=@_;
    my $vdbc = $vdba->dbc();
    my $queries = {
        sift => 'select count(*) from protein_function_predictions;',
        ancestral_allele => 'select variation_id from variation where ancestral_allele is not null limit 1;',
        global_maf => 'select variation_id from variation where minor_allele is not null limit 1;',
        clinical_significance => 'select variation_id from variation where clinical_significance is not null limit 1;',
        clinical_significance_svs => 'select structural_variation_id from structural_variation where clinical_significance is not null limit 1;',
        svs => 'select count(*) from structural_variation;',
    };

    foreach my $data_type (keys %$queries) {
        my $sth = $vdbc->prepare($queries->{$data_type});
        $sth->execute();
        while (my @row = $sth->fetchrow_array) {
            my $count = $row[0];
            if ($count > 0) {
              $config->{$species}->{$data_type} = 1;
            }
        }
        $sth->finish();
    }
    $vdbc->disconnect_if_idle();
}


sub division {
    my ($self, $dba) = @_;
    my ($division) = @{$dba->get_MetaContainer()->list_value_by_key('species.division')};
    return if ! $division;
    $division =~ s/^Ensembl//;

return lc($division);
}

1;
