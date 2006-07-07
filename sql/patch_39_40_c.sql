#to add index and patch column length
alter table allele_group modify allele_group_id int(10) unsigned not null auto_increment;

alter table transcript_variation add key consequence_type_idx(consequence_type);

INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_39_40_c.sql|key_consequence_type');
