# Foreign key relationships in the Ensembl schema (see table.sql);
#
# This file is intended as a reference since some of the relationships
# are not obvious.
#
# Note that these constraints are not actually used by Ensembl for 
# performance reasons, and referential integrity is enforced at the
# application level. Also MySQL currently does not support foreign
# key constraints on MyISAM tables.

ALTER TABLE allele ADD FOREIGN KEY (subsnp_id) REFERENCES subsnp_handle(subsnp_id);
ALTER TABLE allele ADD FOREIGN KEY (variation_id) REFERENCES variation(variation_id);
ALTER TABLE allele ADD FOREIGN KEY (sample_id) REFERENCES population(sample_id);

ALTER TABLE allele_group ADD FOREIGN KEY (source_id) REFERENCES source(source_id);
ALTER TABLE allele_group ADD FOREIGN KEY (variation_group_id) REFERENCES variation_group(variation_group_id);
ALTER TABLE allele_group ADD FOREIGN KEY (sample_id) REFERENCES population(sample_id);

ALTER TABLE allele_group_allele ADD FOREIGN KEY (variation_id) REFERENCES variation(variation_id);
ALTER TABLE allele_group_allele ADD FOREIGN KEY (allele_group_id) REFERENCES allele_group(allele_group_id);

ALTER TABLE compressed_genotype_single_bp ADD FOREIGN KEY (sample_id) REFERENCES individual(sample_id);
ALTER TABLE compressed_genotype_single_bp ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER TABLE failed_variation ADD FOREIGN KEY (variation_id) REFERENCES variation(variation_id);
ALTER TABLE failed_variation ADD FOREIGN KEY (failed_description_id) REFERENCES failed_description(failed_description_id);

ALTER TABLE flanking_sequence ADD FOREIGN KEY (variation_id) REFERENCES variation(variation_id);
ALTER TABLE flanking_sequence ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER TABLE httag ADD FOREIGN KEY (source_id) REFERENCES source(source_id);
ALTER TABLE httag ADD FOREIGN KEY (variation_group_id) REFERENCES variation_group(variation_group_id);

ALTER TABLE individual ADD FOREIGN KEY (sample_id) REFERENCES sample(sample_id);
ALTER TABLE individual ADD FOREIGN KEY (father_individual_sample_id) REFERENCES sample(sample_id);
ALTER TABLE individual ADD FOREIGN KEY (mother_individual_sample_id) REFERENCES sample(sample_id);
ALTER TABLE individual ADD FOREIGN KEY (individual_type_id) REFERENCES individual_type(individual_type_id);

ALTER TABLE individual_genotype_multiple_bp ADD FOREIGN KEY (subsnp_id) REFERENCES subsnp_handle(subsnp_id);
ALTER TABLE individual_genotype_multiple_bp ADD FOREIGN KEY (variation_id) REFERENCES variation(variation_id);
ALTER TABLE individual_genotype_multiple_bp ADD FOREIGN KEY (sample_id) REFERENCES individual(sample_id);

ALTER TABLE individual_population ADD FOREIGN KEY (individual_sample_id) REFERENCES individual(sample_id);
ALTER TABLE individual_population ADD FOREIGN KEY (population_sample_id) REFERENCES population(sample_id);

ALTER TABLE population ADD FOREIGN KEY (sample_id) REFERENCES sample(sample_id);

ALTER TABLE population_genotype ADD FOREIGN KEY (subsnp_id) REFERENCES subsnp_handle(subsnp_id);
ALTER TABLE population_genotype ADD FOREIGN KEY (variation_id) REFERENCES variation(variation_id);
ALTER TABLE population_genotype ADD FOREIGN KEY (sample_id) REFERENCES population(sample_id);

ALTER TABLE population_structure ADD FOREIGN KEY (super_population_sample_id) REFERENCES population(sample_id);
ALTER TABLE population_structure ADD FOREIGN KEY (sub_population_sample_id) REFERENCES population(sample_id);

ALTER TABLE read_coverage ADD FOREIGN KEY (sample_id) REFERENCES individual(sample_id);
ALTER TABLE read_coverage ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER TABLE sample_synonym ADD FOREIGN KEY (source_id) REFERENCES source(source_id);
ALTER TABLE sample_synonym ADD FOREIGN KEY (sample_id) REFERENCES sample(sample_id);

ALTER TABLE structural_variation ADD FOREIGN KEY (source_id) REFERENCES source(source_id);
ALTER TABLE structural_variation ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER TABLE tagged_variation_feature ADD FOREIGN KEY (variation_feature_id) REFERENCES variation_feature(variation_feature_id);
ALTER TABLE tagged_variation_feature ADD FOREIGN KEY (sample_id) REFERENCES sample(sample_id);

ALTER TABLE transcript_variation ADD FOREIGN KEY (variation_feature_id) REFERENCES variation_feature(variation_feature_id);

ALTER TABLE variation ADD FOREIGN KEY (source_id) REFERENCES source(source_id);

ALTER TABLE variation_annotation ADD FOREIGN KEY (source_id) REFERENCES source(source_id);
ALTER TABLE variation_annotation ADD FOREIGN KEY (variation_id) REFERENCES variation(variation_id);
ALTER TABLE variation_annotation ADD FOREIGN KEY (phenotype_id) REFERENCES phenotype(phenotype_id);

ALTER TABLE variation_feature ADD FOREIGN KEY (source_id) REFERENCES source(source_id);
ALTER TABLE variation_feature ADD FOREIGN KEY (variation_id) REFERENCES variation(variation_id);
ALTER TABLE variation_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER TABLE variation_group ADD FOREIGN KEY (source_id) REFERENCES source(source_id);

ALTER TABLE variation_group_feature ADD FOREIGN KEY (variation_group_id) REFERENCES variation_group(variation_group_id);
ALTER TABLE variation_group_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER TABLE variation_group_variation ADD FOREIGN KEY (variation_group_id) REFERENCES variation_group(variation_group_id);
ALTER TABLE variation_group_variation ADD FOREIGN KEY (variation_id) REFERENCES variation(variation_id);

ALTER TABLE variation_set_structure ADD FOREIGN KEY (variation_set_super) REFERENCES variation_set(variation_set_id);
ALTER TABLE variation_set_structure ADD FOREIGN KEY (variation_set_sub) REFERENCES variation_set(variation_set_id);

ALTER TABLE variation_set_variation ADD FOREIGN KEY (variation_id) REFERENCES variation(variation_id);
ALTER TABLE variation_set_variation ADD FOREIGN KEY (variation_set_id) REFERENCES variation_set(variation_set_id);

ALTER TABLE variation_synonym ADD FOREIGN KEY (subsnp_id) REFERENCES subsnp_handle(subsnp_id);
ALTER TABLE variation_synonym ADD FOREIGN KEY (source_id) REFERENCES source(source_id);
ALTER TABLE variation_synonym ADD FOREIGN KEY (variation_id) REFERENCES variation(variation_id);
