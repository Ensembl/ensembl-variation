# add MAF columns to variation_feature

# create tmp table and add columns
CREATE TABLE tmp_vf LIKE variation_feature;
ALTER TABLE tmp_vf ADD COLUMN `minor_allele` char(1) DEFAULT NULL;
ALTER TABLE tmp_vf ADD COLUMN `minor_allele_freq` float DEFAULT NULL;
ALTER TABLE tmp_vf ADD COLUMN `minor_allele_count` int(10) unsigned DEFAULT NULL;

# copy data with join to variation
INSERT INTO tmp_vf
SELECT vf.*, v.minor_allele, v.minor_allele_freq, v.minor_allele_count
FROM variation v, variation_feature vf
WHERE v.variation_id = vf.variation_id;

# drop and rename
DROP TABLE variation_feature;
RENAME TABLE tmp_vf TO variation_feature;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_67_68_e.sql|add MAF columns to variation_feature');