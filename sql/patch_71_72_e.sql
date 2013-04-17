# replace variation.clinical_significance_attrib_id with variation.clinical_significance set to allow different statuses for different traits

alter table variation add column clinical_significance SET('drug-response','histocompatibility','non-pathogenic','other','pathogenic','probable-non-pathogenic','probable-pathogenic''unknown','untested');


alter table variation drop column clinical_significance_attrib_id;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_71_72_e.sql|change variation clinical_significance column');
