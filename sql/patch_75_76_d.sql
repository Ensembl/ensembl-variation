# create read_coverage table
CREATE TABLE IF NOT EXISTS read_coverage (
  seq_region_id int(10) unsigned not null,
  seq_region_start int not null,
  seq_region_end int not null,
  level tinyint not null,
  individual_id int(10) unsigned not null,
  
  key seq_region_idx(seq_region_id,seq_region_start)   
);

#patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_75_76_d.sql|Restore read_coverage table');
