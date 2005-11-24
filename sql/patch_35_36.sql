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
