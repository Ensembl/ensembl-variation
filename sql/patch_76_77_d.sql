# update deprecated SO terms in attrib table

UPDATE attrib SET value='non_coding_transcript_variant' WHERE value='nc_transcript_variant';
UPDATE attrib SET value='non_coding_transcript_exon_variant' WHERE value='non_coding_exon_variant';


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_76_77_d.sql|update SO terms in attrib table');
