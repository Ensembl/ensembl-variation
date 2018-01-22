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

package Bio::EnsEMBL::Variation::Pipeline::VariantQC::UpdatePopulationGenotype;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);


=head1 NAME

  Bio::EnsEMBL::Variation::Pipeline::VariantQC::UpdatePopulationGenotype

=head1 DESCRIPTION

This module migrates the population genotype table from old to new schema version
It is run as a seperate independant process after the main variant QC 
Reports back to hive db when processes complete in case of job time outs (the hive
job may die, but the database process continue to completion 

=cut

sub run {
    
   my $self = shift;
   
   
   
   my $var_dba      = $self->get_species_adaptor('variation');


  ## populate temp table with genotypes from sample_genotype_multiple_bp
  $var_dba->dbc->do(qq[INSERT IGNORE INTO genotype_code_tmp(allele_1, allele_2, phased) 
                      SELECT distinct allele_1, allele_2, 0 FROM sample_genotype_multiple_bp]);

  ## populate temp table with genotypes from flipped population_genotype table 
  $var_dba->dbc->do(qq[INSERT IGNORE INTO genotype_code_tmp(allele_1, allele_2, phased) 
                      SELECT distinct allele_1, allele_2, 0 FROM MTMP_population_genotype_working]);

  ## add any missing allele codes present in genotype tables only
  $var_dba->dbc->do(qq[INSERT IGNORE INTO allele_code(allele) SELECT allele_1 FROM genotype_code_tmp]);
  $var_dba->dbc->do(qq[INSERT IGNORE INTO allele_code(allele) SELECT allele_2 FROM genotype_code_tmp]);

  ## populate genotype code with both alleles
  $var_dba->dbc->do(qq[INSERT INTO genotype_code (genotype_code_id, allele_code_id, haplotype_id, phased )
                       SELECT t.genotype_code_id, ac.allele_code_id, 1, phased
                       FROM genotype_code_tmp t, allele_code ac
                       WHERE t.allele_1 = ac.allele ]);

  $var_dba->dbc->do(qq[INSERT INTO genotype_code (genotype_code_id, allele_code_id, haplotype_id, phased )
                      SELECT t.genotype_code_id, ac.allele_code_id, 2, phased
                      FROM genotype_code_tmp t, allele_code ac
                      WHERE t.allele_2 = ac.allele ]);
  $self->warning( 'UpdatePopulationGenotype: genotype_code updated');

  $var_dba->dbc->do(qq[ALTER TABLE genotype_code ORDER BY genotype_code_id, haplotype_id ASC]);
  $self->warning( 'UpdatePopulationGenotype: genotype_code sorted');

  ## Create coded genotypes from Mart table in which minus-strand single-mapping variants have been flipped
  $var_dba->dbc->do(qq{ ALTER TABLE population_genotype_working DISABLE KEYS});

  $var_dba->dbc->do(qq[insert into population_genotype_working 
                      select pg.population_genotype_id, pg.variation_id, pg.subsnp_id,  
                      gc.genotype_code_id, pg.frequency, pg.population_id, pg.count
                      from MTMP_population_genotype_working pg, genotype_code_tmp gc 
                      where pg.allele_1 = gc.allele_1 and pg.allele_2 = gc.allele_2 
                      and gc.phased =0 ]);
   $self->warning( 'UpdatePopulationGenotype: population_genotype_working populated');

   $var_dba->dbc->do(qq{ ALTER TABLE population_genotype_working ENABLE KEYS});

}



1;
