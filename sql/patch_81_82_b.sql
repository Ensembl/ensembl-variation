# update the description in the failed_description table
UPDATE failed_description SET description='Variant maps to more than 3 different locations' WHERE description='Variation maps to more than 3 different locations';
UPDATE failed_description SET description='Variant has more than 3 different alleles' WHERE description='Variation has more than 3 different alleles';
UPDATE failed_description SET description='Variant does not map to the genome' WHERE description='Variation does not map to the genome';
UPDATE failed_description SET description='Variant has no genotypes' WHERE description='Variation has no genotypes';
UPDATE failed_description SET description='Variant has no associated sequence' WHERE description='Variation has no associated sequence';
UPDATE failed_description SET description='Variant submission has been withdrawn by the 1000 genomes project due to high false positive rate' WHERE description='Variation submission has been withdrawn by the 1000 genomes project due to high false positive rate';
UPDATE failed_description SET description='Variant has more than 3 different submitted alleles' WHERE description='Variation has more than 3 different submitted alleles';
UPDATE failed_description SET description='Variant can not be re-mapped to the current assembly' WHERE description='Variation can not be re-mapped to the current assembly';
UPDATE failed_description SET description='Variant maps to more than one genomic location' WHERE description='Variation maps to more than one genomic location';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_81_82_b.sql|update the description in the failed_description table');

