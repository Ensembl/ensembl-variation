=head1 LICENSE

 Copyright (c) 1999-2012 The European Bioinformatics Institute and
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

package Bio::EnsEMBL::Variation::Pipeline::VariantQC::UpdatePopulationGenotype;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);


=head1 NAME

  Bio::EnsEMBL::Variation::Pipeline::VariantQC::UpdatePopulationGenotype

=head1 DESCRIPTION

This module migrates the population genotype table from old to new schema version
It is run as a seperate independant process after the main variant QC 

=cut

sub run {
    
   my $self = shift;
   
   
   
   my $var_dba      = $self->get_species_adaptor('variation');


  ## populate temp table with genotypes from individual_genotype_multiple_bp
  $var_dba->dbc->do(qq[INSERT IGNORE INTO genotype_code_tmp(allele_1, allele_2) 
                      SELECT distinct allele_1, allele_2 FROM individual_genotype_multiple_bp]);

  ## populate temp table with genotypes from flipped population_genotype table 
  $var_dba->dbc->do(qq[INSERT IGNORE INTO genotype_code_tmp(allele_1, allele_2) 
                      SELECT distinct allele_1, allele_2 FROM MTMP_population_genotype_working]);

  ## add any missing allele codes present in genotype tables only
  $var_dba->dbc->do(qq[INSERT IGNORE INTO allele_code(allele) SELECT allele_1 FROM genotype_code_tmp]);
  $var_dba->dbc->do(qq[INSERT IGNORE INTO allele_code(allele) SELECT allele_2 FROM genotype_code_tmp]);

  ## populate genotype code with both alleles
  $var_dba->dbc->do(qq[INSERT INTO genotype_code
                       SELECT t.genotype_code_id, ac.allele_code_id, 1
                       FROM genotype_code_tmp t, allele_code ac
                       WHERE t.allele_1 = ac.allele ]);

  $var_dba->dbc->do(qq[INSERT INTO genotype_code
                      SELECT t.genotype_code_id, ac.allele_code_id, 2
                      FROM genotype_code_tmp t, allele_code ac
                      WHERE t.allele_2 = ac.allele ]);

  $var_dba->dbc->do(qq[ALTER TABLE genotype_code ORDER BY genotype_code_id, haplotype_id ASC]);

  ## Create coded genotypes from Mart table in which minus-strand single-mapping variants have been flipped
   $var_dba->dbc->do(qq{ ALTER TABLE population_genotype_working DISABLE KEYS});

  $var_dba->dbc->do(qq[insert into population_genotype_working 
                      select pg.population_genotype_id, pg.variation_id, pg.subsnp_id,  
                      gc.genotype_code_id, pg.frequency, pg.sample_id, pg.count
                      from MTMP_population_genotype_working pg, genotype_code_tmp gc 
                      where pg.allele_1 = gc.allele_1 and pg.allele_2 = gc.allele_2 ]);
   
 $var_dba->dbc->do(qq{ ALTER TABLE population_genotype_working ENABLE KEYS});
}



1;
