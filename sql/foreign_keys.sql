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
ALTER TABLE allele ADD FOREIGN KEY (allele_code_id) REFERENCES allele_code(allele_code_id);
ALTER TABLE allele ADD FOREIGN KEY (population_id) REFERENCES population(population_id);
ALTER TABLE allele ADD FOREIGN KEY (frequency_submitter_handle) REFERENCES submitter_handle(handle_id);

ALTER TABLE associate_study ADD FOREIGN KEY (study1_id) REFERENCES study(study_id);
ALTER TABLE associate_study ADD FOREIGN KEY (study2_id) REFERENCES study(study_id);

ALTER TABLE attrib ADD FOREIGN KEY (attrib_type_id) REFERENCES attrib_type(attrib_type_id);

ALTER TABLE attrib_set ADD FOREIGN KEY (attrib_id) REFERENCES attrib(attrib_id);

ALTER TABLE compressed_genotype_region ADD FOREIGN KEY (individual_id) REFERENCES individual(individual_id);
ALTER TABLE compressed_genotype_region ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER TABLE compressed_genotype_var ADD FOREIGN KEY (variation_id) REFERENCES variation(variation_id);
ALTER TABLE compressed_genotype_var ADD FOREIGN KEY (subsnp_id) REFERENCES subsnp_handle(subsnp_id);

ALTER TABLE failed_allele ADD FOREIGN KEY (allele_id) REFERENCES allele(allele_id);
ALTER TABLE failed_allele ADD FOREIGN KEY (failed_description_id) REFERENCES failed_description(failed_description_id);

ALTER TABLE failed_variation ADD FOREIGN KEY (variation_id) REFERENCES variation(variation_id);
ALTER TABLE failed_variation ADD FOREIGN KEY (failed_description_id) REFERENCES failed_description(failed_description_id);

ALTER TABLE failed_structural_variation ADD FOREIGN KEY (structural_variation_id) REFERENCES structural_variation(structural_variation_id);
ALTER TABLE failed_structural_variation ADD FOREIGN KEY (failed_description_id) REFERENCES failed_description(failed_description_id);

ALTER TABLE genotype_code ADD FOREIGN KEY (allele_code_id) REFERENCES allele_code(allele_code_id);

ALTER TABLE individual ADD FOREIGN KEY (father_individual_id) REFERENCES individual(individual_id);
ALTER TABLE individual ADD FOREIGN KEY (mother_individual_id) REFERENCES individual(individual_id);
ALTER TABLE individual ADD FOREIGN KEY (individual_type_id) REFERENCES individual_type(individual_type_id);

ALTER TABLE individual_genotype_multiple_bp ADD FOREIGN KEY (subsnp_id) REFERENCES subsnp_handle(subsnp_id);
ALTER TABLE individual_genotype_multiple_bp ADD FOREIGN KEY (variation_id) REFERENCES variation(variation_id);
ALTER TABLE individual_genotype_multiple_bp ADD FOREIGN KEY (individual_id) REFERENCES individual(individual_id);

ALTER TABLE individual_population ADD FOREIGN KEY (individual_id) REFERENCES individual(individual_id);
ALTER TABLE individual_population ADD FOREIGN KEY (population_id) REFERENCES population(population_id);

ALTER TABLE population_genotype ADD FOREIGN KEY (subsnp_id) REFERENCES subsnp_handle(subsnp_id);
ALTER TABLE population_genotype ADD FOREIGN KEY (variation_id) REFERENCES variation(variation_id);
ALTER TABLE population_genotype ADD FOREIGN KEY (genotype_code_id) REFERENCES genotype_code(genotype_code_id);
ALTER TABLE population_genotype ADD FOREIGN KEY (population_id) REFERENCES population(population_id);

ALTER TABLE population_structure ADD FOREIGN KEY (super_population_id) REFERENCES population(population_id);
ALTER TABLE population_structure ADD FOREIGN KEY (sub_population_id) REFERENCES population(population_id);

ALTER TABLE protein_function_predictions ADD FOREIGN KEY (translation_md5_id) REFERENCES translation_md5(translation_md5_id);
ALTER TABLE protein_function_predictions ADD FOREIGN KEY (analysis_attrib_id) REFERENCES attrib(attrib_id);

ALTER TABLE read_coverage ADD FOREIGN KEY (individual_id) REFERENCES individual(individual_id);
ALTER TABLE read_coverage ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER TABLE individual_synonym ADD FOREIGN KEY (source_id) REFERENCES source(source_id);
ALTER TABLE individual_synonym ADD FOREIGN KEY (individual_id) REFERENCES individual(individual_id);

ALTER TABLE population_synonym ADD FOREIGN KEY (source_id) REFERENCES source(source_id);
ALTER TABLE population_synonym ADD FOREIGN KEY (population_id) REFERENCES population(population_id);

ALTER TABLE seq_region ADD FOREIGN KEY (coord_system_id) REFERENCES coord_system(coord_system_id);

ALTER TABLE structural_variation ADD FOREIGN KEY (source_id) REFERENCES source(source_id);
ALTER TABLE structural_variation ADD FOREIGN KEY (study_id) REFERENCES study(study_id);
ALTER TABLE structural_variation ADD FOREIGN KEY (class_attrib_id) REFERENCES attrib(attrib_id);

ALTER TABLE structural_variation_association ADD FOREIGN KEY (structural_variation_id) REFERENCES structural_variation(structural_variation_id);
ALTER TABLE structural_variation_association ADD FOREIGN KEY (supporting_structural_variation_id) REFERENCES structural_variation(structural_variation_id);

ALTER TABLE structural_variation_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);
ALTER TABLE structural_variation_feature ADD FOREIGN KEY (structural_variation_id) REFERENCES structural_variation(structural_variation_id);
ALTER TABLE structural_variation_feature ADD FOREIGN KEY (source_id) REFERENCES source(source_id);
ALTER TABLE structural_variation_feature ADD FOREIGN KEY (study_id) REFERENCES study(study_id);
ALTER TABLE structural_variation_feature ADD FOREIGN KEY (class_attrib_id) REFERENCES attrib(attrib_id);

ALTER TABLE structural_variation_sample ADD FOREIGN KEY (individual_id) REFERENCES individual(individual_id);
ALTER TABLE structural_variation_sample ADD FOREIGN KEY (strain_id) REFERENCES individual(individual_id);
ALTER TABLE structural_variation_sample ADD FOREIGN KEY (structural_variation_id) REFERENCES structural_variation(structural_variation_id);


ALTER TABLE study ADD FOREIGN KEY (source_id) REFERENCES source(source_id);

ALTER TABLE study_variation ADD FOREIGN KEY (study_id) REFERENCES study(study_id);
ALTER TABLE study_variation ADD FOREIGN KEY (variation_id) REFERENCES variation(variation_id);

ALTER TABLE tagged_variation_feature ADD FOREIGN KEY (variation_feature_id) REFERENCES variation_feature(variation_feature_id);
ALTER TABLE tagged_variation_feature ADD FOREIGN KEY (population_id) REFERENCES population(population_id);

ALTER TABLE transcript_variation ADD FOREIGN KEY (variation_feature_id) REFERENCES variation_feature(variation_feature_id);

ALTER TABLE motif_feature_variation ADD FOREIGN KEY (variation_feature_id) REFERENCES variation_feature(variation_feature_id);

ALTER TABLE regulatory_feature_variation ADD FOREIGN KEY (variation_feature_id) REFERENCES variation_feature(variation_feature_id);

ALTER TABLE variation ADD FOREIGN KEY (source_id) REFERENCES source(source_id);
ALTER TABLE variation ADD FOREIGN KEY (class_attrib_id) REFERENCES attrib(attrib_id);

ALTER TABLE phenotype_citation ADD FOREIGN KEY (phenotype_id) REFERENCES phenotype(phenotype_id);
ALTER TABLE phenotype_citation ADD FOREIGN KEY (publication_id) REFERENCES publication(publication_id);

ALTER TABLE phenotype_feature ADD FOREIGN KEY (source_id) REFERENCES study(study_id);
ALTER TABLE phenotype_feature ADD FOREIGN KEY (study_id) REFERENCES study(study_id);
ALTER TABLE phenotype_feature ADD FOREIGN KEY (object_id) REFERENCES variation(name);
ALTER TABLE phenotype_feature ADD FOREIGN KEY (object_id) REFERENCES structural_variation(variation_name);
ALTER TABLE phenotype_feature ADD FOREIGN KEY (phenotype_id) REFERENCES phenotype(phenotype_id);
ALTER TABLE phenotype_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER TABLE phenotype_feature_attrib ADD FOREIGN KEY (phenotype_feature_id) REFERENCES phenotype_feature(phenotype_feature_id);
ALTER TABLE phenotype_feature_attrib ADD FOREIGN KEY (attrib_type_id) REFERENCES attrib_type(attrib_type_id);

ALTER TABLE variation_citation ADD FOREIGN KEY (variation_id) REFERENCES variation(variation_id);
ALTER TABLE variation_citation ADD FOREIGN KEY (publication_id) REFERENCES publication(publication_id);

ALTER TABLE variation_feature ADD FOREIGN KEY (source_id) REFERENCES source(source_id);
ALTER TABLE variation_feature ADD FOREIGN KEY (variation_id) REFERENCES variation(variation_id);
ALTER TABLE variation_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);
ALTER TABLE variation_feature ADD FOREIGN KEY (class_attrib_id) REFERENCES attrib(attrib_id);

ALTER TABLE variation_set ADD FOREIGN KEY (short_name_attrib_id) REFERENCES attrib(attrib_id);

ALTER TABLE variation_set_structural_variation ADD FOREIGN KEY (structural_variation_id) REFERENCES structural_variation(structural_variation_id);
ALTER TABLE variation_set_structural_variation ADD FOREIGN KEY (variation_set_id) REFERENCES variation_set(variation_set_id);

ALTER TABLE variation_set_structure ADD FOREIGN KEY (variation_set_super) REFERENCES variation_set(variation_set_id);
ALTER TABLE variation_set_structure ADD FOREIGN KEY (variation_set_sub) REFERENCES variation_set(variation_set_id);

ALTER TABLE variation_set_variation ADD FOREIGN KEY (variation_id) REFERENCES variation(variation_id);
ALTER TABLE variation_set_variation ADD FOREIGN KEY (variation_set_id) REFERENCES variation_set(variation_set_id);

ALTER TABLE variation_synonym ADD FOREIGN KEY (subsnp_id) REFERENCES subsnp_handle(subsnp_id);
ALTER TABLE variation_synonym ADD FOREIGN KEY (source_id) REFERENCES source(source_id);
ALTER TABLE variation_synonym ADD FOREIGN KEY (variation_id) REFERENCES variation(variation_id);
