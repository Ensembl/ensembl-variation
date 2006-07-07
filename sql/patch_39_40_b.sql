#this patch is suposed to change the order in the set column

ALTER TABLE variation_feature ADD column new_consequence_type SET ('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST','FRAMESHIFT_CODING','NON_SYNONYMOUS_CODING','SPLICE_SITE','SYNONYMOUS_CODING','REGULATORY_REGION','5PRIME_UTR','3PRIME_UTR','INTRONIC','UPSTREAM','DOWNSTREAM','INTERGENIC') DEFAULT 'INTERGENIC' not null;
UPDATE variation_feature SET new_consequence_type = consequence_type;
ALTER TABLE variation_feature DROP COLUMN consequence_type, CHANGE new_consequence_type consequence_type  SET ('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST','FRAMESHIFT_CODING','NON_SYNONYMOUS_CODING','SPLICE_SITE','SYNONYMOUS_CODING','REGULATORY_REGION','5PRIME_UTR','3PRIME_UTR','INTRONIC','UPSTREAM','DOWNSTREAM','INTERGENIC') DEFAULT 'INTERGENIC' not null;

#do the same in the trancript_variation table

ALTER TABLE transcript_variation ADD column new_consequence_type SET ('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST','FRAMESHIFT_CODING','NON_SYNONYMOUS_CODING','SPLICE_SITE','SYNONYMOUS_CODING','REGULATORY_REGION','5PRIME_UTR','3PRIME_UTR','INTRONIC','UPSTREAM','DOWNSTREAM') not null;
UPDATE transcript_variation SET new_consequence_type = consequence_type;
ALTER TABLE transcript_variation DROP COLUMN consequence_type, CHANGE new_consequence_type consequence_type  SET ('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST','FRAMESHIFT_CODING','NON_SYNONYMOUS_CODING','SPLICE_SITE','SYNONYMOUS_CODING','REGULATORY_REGION','5PRIME_UTR','3PRIME_UTR','INTRONIC','UPSTREAM','DOWNSTREAM') not null;

INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_39_40_b.sql|modify_set_column');
