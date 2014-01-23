#
# Document describing how to create a Variation database and import it into a local copy
# v1.0
# 

First of all, you should be familiar with the database schema. I would recommend to have a look at the schema under ensembl-variation/schema/database_schema.pdf, so you familiarise with the tables, columns and relations.

Then, you will need to create and empty database. You will have to run the ensembl-variation/sql/tables.sql in your database to create a ensembl-variation flavour database.

After that, you can start importing your data in the database. I would have a look at the script, ensembl-variation/scripts/import/import_Sanger_database.pl, this script copies the data from a new variation database into an existing one, and will tell you which tables you need to put information. In summary:

Sample table -- you need to create a new Sample for each Population or Individual you have in your data. If you have genotypes for 5 different individuals, you will need to create 5 entries for the Individual, and 1 entry for the Population they belong to (assuming they all belong to the same Population). The necessary field is the name.
Source table -- this is a sinlge entry that refers to who called the SNPs. In that case, it will probably be your institute. The necessary field is the name
Population table -- for the Sample you have previously created that was the population, simply add in this table the sample_id
Individual table -- for the Sample you have created that are individuals, add the sample_id here plus the individual_type_id (1 for fully inbred organism, 2 for partly_inbred, 3 for outbred and 4 for mutant).
IndividualPopulation table-- this table should contain the relation between population and individual entries. Following the above example, if you have 5 individuals belonging to the same population, you would need to add 5 entries in the table, with the relation between individual_sample_id and population_sample_id
Variation table -- here you import all Variations you have called. Source_id will probably be the same for all of them. You need to assign a name (you can make up one if you don't have it).
Allele table--containing the alleles in the Populations for the different Variations. You will need a variation_id, the 2 (or more) alleles as text, and the sample_id for the population they belong. You don't have to put frequency if you don't have this value.
FlankingSequence table--if you are using the same coordinate system we have our core database, you should just fill in the values in variation_id, to refer to the variation contained in the flanking sequence, and the columns up_seq_region_start, up_seq_region_end, down_seq_region_start, down_seq_region_end, seq_region_id and seq_region_strand. The first 4 are coordinates for the upstream sequence and the downstream sequences. Seq_region_id contains the region in which the coordinates are (you should get this value from our core database), and the seq_region_strand should contain either 1 or -1, indicating if they are in the positive or negative strand.
VariationFeature table--here you store the position of the variations. Seq_region_id, seq_region_start, seq_region_end and seq_region_strand define the variation in your region, variation_id should be the same you used in the variation table, allele_string contains a non-repetitive concatenation of the alleles, having the first as the reference (so, A/T would mean that the A is the reference allele). variation_name and source_id are the same you assigned in the variation table. map_weight contains the number of times this variation maps to the genome. If you have just 1 mapping for each variation (as it should be, but dbSNP contains weird things...), just put 1 here.
tmp_individual_genotype_single_bp table--notice that this table doesn't exist in the tables.sql file, because it is a temporary table used to initially store the genotypes with 1 base pair, but it is later compressed in the compressed_genotype_single_bp table. You can create it writing in your mysql: CREATE TABLE tmp_individual_genotype_single_bp LIKE individual_genotype_multiple_bp. Just notice again that this should only store 1 base genotypes, so if you have any indel or similar, they should go into the multiple_bp table. The format here is quite simple, you put the variation_id, the 2 alleles in your genotype (it doesn't matter the order) and the individual they have been genotyped (sample_id refers to the individual one). Once this table is populated, you will need to compress the data using the compress_genotypes.pl script. You need to have all the ensembl API (ensembl and ensembl-variation) in order to run it. If you have it, running like perl compress_genotypes.pl -species 'human' -tmpdir '/my/tmp/dir' -tmpfile 'myfile.txt' should compress the data and load it in the compressed_genotype table.

This is what would be the "core" information necessary to have a useful (and usable) API and website, even you might want to add more things (like transcript_variation information or read_coverage data).

In any case, once you have your database, you can use the import_Sanger_database.pl script to import this database in your local copy. This script should take care of all the id headaches.

If you have any problems/questions creating the database, don't hesitate to contact us at http://lists.ensembl.org/mailman/listinfo/dev
