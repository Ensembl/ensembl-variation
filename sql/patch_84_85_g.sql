-- allow the column description to store more text in the source table
ALTER TABLE source CHANGE description description varchar(400);

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_84_85_g.sql|allow the column description to store more text in the source table');
