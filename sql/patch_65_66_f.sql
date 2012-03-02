# update the tagged_variation_feature table
##################
alter table tagged_variation_feature drop primary key;

alter table tagged_variation_feature
add column tagged_variation_feature_id int(10) unsigned default null
after variation_feature_id;

# create indices
create index tag_idx on tagged_variation_feature(variation_feature_id);
create index tagged_idx on tagged_variation_feature(tagged_variation_feature_id);
create index sample_idx on tagged_variation_feature(sample_id);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_65_66_f.sql|change tagged_variation_feature to store relationship between tag and tagged');
