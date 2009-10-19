# added subsnp_id in the  variation_synonym, allele, population_genotype, tmp_individual_genotype_single_bp, individual_genotype_multiple_bp table
##################
alter table variation_synonym add subsnp_id int(15) unsigned after variation_id, add key subsnp_idx(subsnp_id);
alter table allele add subsnp_id int(15) unsigned after variation_id, add key subsnp_idx(subsnp_id);
alter table population_genotype add subsnp_id int(15) unsigned after variation_id, add key subsnp_idx(subsnp_id);
alter table tmp_individual_genotype_single_bp add subsnp_id int(15) unsigned after variation_id, add key subsnp_idx(subsnp_id);
alter table individual_genotype_multiple_bp add subsnp_id int(15) unsigned after variation_id, add key subsnp_idx(subsnp_id);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_55_56_b.sql|add subsnp_id');
