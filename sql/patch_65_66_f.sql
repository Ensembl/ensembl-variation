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
