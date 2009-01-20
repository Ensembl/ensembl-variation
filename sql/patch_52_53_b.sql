# patch_52_53_b.sql
#
# title: change variation_annotation table
#
# description:
# Table containing annotation associated with the variation
# such as GWAS

alter table variation_annotation 
add column study varchar(30) default NULL after source_id,
add column associated_gene varchar(255) default NULL,
add column associated_variant_risk_allele varchar(255) default NULL,
add column variation_names varchar(255) default NULL,
add column risk_allele_freq_in_controls varchar(30) default NULL,
add column p_value varchar(20) default NULL;




INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_52_53_b.sql|add columns in variation_annotation table');
