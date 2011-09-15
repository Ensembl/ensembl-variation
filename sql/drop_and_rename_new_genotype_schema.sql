# script to drop and rename tables after new genotype schema conversion

# drop tables
drop table allele;
drop table population_genotype;
drop table compressed_genotype_single_bp;
drop table genotype_code_tmp;

# rename tables
rename table allele_proxy to allele;
rename table population_genotype_proxy to population_genotype;