-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2018] EMBL-European Bioinformatics Institute
-- 
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
-- 
--      http://www.apache.org/licenses/LICENSE-2.0
-- 
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.


## Add two more consequence terms (regulatory_region_variant, TF_binding_site_variant) to consequence_types set in variation_feature table

#Create backup table
CREATE TABLE variation_feature_backup like variation_feature;
INSERT INTO variation_feature_backup SELECT * FROM variation_feature;

# Create tmp table
CREATE TABLE variation_feature_tmp like variation_feature;
ALTER TABLE variation_feature_tmp CHANGE consequence_types consequence_types SET('intergenic_variant','splice_acceptor_variant','splice_donor_variant','stop_lost','coding_sequence_variant','missense_variant','stop_gained','synonymous_variant','frameshift_variant','nc_transcript_variant','non_coding_exon_variant','mature_miRNA_variant','NMD_transcript_variant','5_prime_UTR_variant','3_prime_UTR_variant','incomplete_terminal_codon_variant','intron_variant','splice_region_variant','downstream_gene_variant','upstream_gene_variant','initiator_codon_variant','stop_retained_variant','inframe_insertion','inframe_deletion','transcript_ablation','transcript_fusion','transcript_amplification','transcript_translocation','TFBS_ablation','TFBS_fusion','TFBS_amplification','TFBS_translocation','regulatory_region_ablation','regulatory_region_fusion','regulatory_region_amplification','regulatory_region_translocation','feature_elongation','feature_truncation','regulatory_region_variant','TF_binding_site_variant') NOT NULL DEFAULT 'intergenic_variant';

# Insert into tmp table
INSERT INTO variation_feature_tmp SELECT * FROM variation_feature;

DROP TABLE variation_feature;
RENAME TABLE variation_feature_tmp TO variation_feature;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_69_70_e.sql|add consequence terms to consequence_types set in variation_feature');
