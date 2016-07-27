-- add qualifier & index to phenotype_onology_accession

ALTER TABLE phenotype_ontology_accession add index accession_idx(accession);

ALTER TABLE phenotype_ontology_accession CHANGE column linked_by_attrib mapped_by_attrib set('437','438','439','440','441','442','443','444') ;

ALTER TABLE phenotype_ontology_accession add column ( mapping_type enum('is','involves') default NULL ) ;

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_85_86_b.sql|add qualifier & index to phenotype_onology_accession');
