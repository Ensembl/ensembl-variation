## This patch changes all fields that hold alleles to be able to hold longer strings (up to 25000 characters)

# ancestral allele in variation
ALTER TABLE variation CHANGE ancestral_allele ancestral_allele varchar(25000);

# allele_string in variation_feature - this is necessarily longer as it contains more than one allele
ALTER TABLE variation_feature CHANGE allele_string allele_string varchar(50000);

# allele_1 and allele_2 in population_genotype
ALTER TABLE population_genotype CHANGE allele_1 allele_1 varchar(25000);
ALTER TABLE population_genotype CHANGE allele_2 allele_2 varchar(25000);

# allele in allele
ALTER TABLE allele CHANGE allele allele varchar(25000);

# allele_1 and allele_2 in individual_genotype_multiple_bp
ALTER TABLE individual_genotype_multiple_bp CHANGE allele_1 allele_1 varchar(25000);
ALTER TABLE individual_genotype_multiple_bp CHANGE allele_2 allele_2 varchar(25000);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_61_62_j.sql|change fields that hold alleles to large varchars');
