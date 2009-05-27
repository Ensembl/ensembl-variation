# update the consequence_type entry in the variation_feature table
# add a type 'NO_CONSEQUENCE' for illumina cnv entries
##################
alter table variation_feature change consequence_type consequence_type set('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST','COMPLEX_INDEL','FRAMESHIFT_CODING','NON_SYNONYMOUS_CODING','SPLICE_SITE','SYNONYMOUS_CODING','REGULATORY_REGION','WITHIN_MATURE_miRNA','5PRIME_UTR','3PRIME_UTR','INTRONIC','UPSTREAM','DOWNSTREAM','WITHIN_NON_CODING_GENE','NO_CONSEQUENCE','INTERGENIC') default 'INTERGENIC';;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_54_55_b.sql|add NO_CONSEQUENCE for cnv');
