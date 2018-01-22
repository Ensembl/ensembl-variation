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


################
#
# Table structure for table compressed_genotype_single_bp
#
################

CREATE TABLE compressed_genotype_single_bp(
  sample_id int not null,
  seq_region_id int not null,
  seq_region_start int not null,
  seq_region_end int not null,
  seq_region_strand tinyint not null,
  genotypes blob,

  key pos_idx(seq_region_id,seq_region_start)
);

################
#
#  Remove the LD table. LD values are calculated on the fly
#
################

DROP TABLE pairwise_ld;

################
#
# 2 new fields for a variation: ancestral_allele and molecular_type
#
################


ALTER TABLE variation add ancestral_allele text;

ALTER TABLE variation_synonym add moltype varchar(50);

################
#
# And finally, it will be necessary to import the genotype data and compress it
# for doing it, you must run the
# ensembl-variation/scripts/import/compress_genotypes.pl
#
################

SELECT 'In order to populate the Compressed_gentoype_single_bp table, you must
run the script ensembl-variation/scripts/import/compress_genotypes.pl. Then, you
can delete individual_genotype_single_bp. Be aware that the process of
compressing the genotype table can take 7/8 hours in human' AS '';
