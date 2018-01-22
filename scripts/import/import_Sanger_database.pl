#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.



=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut


use strict;
use warnings;

use Getopt::Long;
use ImportUtils qw(debug load dumpSQL);
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);
use DBH;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use FindBin qw( $Bin );
use DBI qw(:sql_types);
use Data::Dumper;
my $species;

GetOptions('tmpdir=s'  => \$ImportUtils::TMP_DIR,
	   'tmpfile=s' => \$ImportUtils::TMP_FILE,
	   'species=s' => \$species
	   );

warn("Make sure you have a updated ensembl.registry file!\n");
die "you must specify the species:" if (!$species);
my $registry_file ||= $Bin . "/ensembl.registry";

Bio::EnsEMBL::Registry->load_all( $registry_file );

my $dbVar = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
my $dbSanger = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'sanger');

my $TMP_DIR  = $ImportUtils::TMP_DIR;
my $LOCAL_TMP_DIR  = '/tmp';
my $TMP_FILE = $ImportUtils::TMP_FILE;

my $old_new_variation_id = {}; #reference to a hash with the old_variation_id => new_variation_id
my $old_new_variation_feature_id = {}; #reference to a hash with the old_Variation_feature_id => new_variation_feature_id
my $old_new_source_id = {}; # reference to a hash with the old_source_id => new_source_id
my $old_new_sample_id = {}; #reference to a hash with the old_sample_id => new_sample_id
my $sanger_sample = {}; #reference to a hash containing the samples that are in the variation database

my $last_source_id = get_last_table_id($dbVar,"source"); #last source_id used in the database
my $last_sample_id = get_last_table_id($dbVar,"sample"); #last sample_id used in the database
my $last_variation_id = get_last_table_id($dbVar,"variation"); #last variation_id used in the database
my $last_variation_feature_id = get_last_table_id($dbVar,"variation_feature"); #last variation_feature_id used in the database

import_Sample_table($dbSanger, $dbVar, $old_new_sample_id, $last_sample_id, $sanger_sample);
import_Source_table($dbSanger, $dbVar, $old_new_source_id, $last_source_id);
import_Population_table($dbSanger, $dbVar, $old_new_sample_id, $sanger_sample);
#import_Individual_table($dbSanger,$dbVar,$old_new_sample_id);
#import_Individual_Population_table($dbSanger,$dbVar,$old_new_sample_id);
#import_Meta_table($dbSanger,$dbVar);
#import_Meta_Coord_table($dbSanger,$dbVar);
import_Variation_table($dbSanger,$dbVar,$old_new_variation_id,$last_variation_id, $old_new_source_id);
import_Allele_table($dbSanger,$dbVar,$old_new_variation_id, $old_new_sample_id);
import_Flanking_sequence_table($dbSanger,$dbVar,$old_new_variation_id);
#import_Variation_synonym_table($dbSanger, $dbVar, $old_new_variation_id, $old_new_source_id);
import_Variation_feature_table($dbSanger,$dbVar,$old_new_variation_feature_id,$last_variation_feature_id, $old_new_variation_id, $old_new_source_id);
import_Transcript_variation_table($dbSanger,$dbVar,$old_new_variation_feature_id);
import_Read_coverage_table($dbSanger,$dbVar, $old_new_sample_id);
import_Tmp_individual_genotype_single_bp_table($dbSanger,$dbVar,$old_new_variation_id,$old_new_sample_id);



debug("Finished reading data, ready to start loading in the database....");


sub import_Sample_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_sample_id = shift;
    my $last_sample_id = shift;
    my $sanger_sample = shift;

    debug("Loading Sample table");
    my ($sample_id, $name, $size, $description);
    my $new_sample_id;
    my $sth = $dbSanger->dbc->prepare(qq{SELECT sample_id, name, size, description from sample});
    $sth->execute();
    $sth->bind_columns(\$sample_id, \$name, \$size, \$description);
    while ($sth->fetch){
	#need to check if the strain is already in the Variation database, as a population or individual
	$new_sample_id = &get_sample_variation_database($dbVar, $name);
	if ($new_sample_id == 0){
	    #get the new id for the sample in the variation table
	    $new_sample_id = $last_sample_id + 1;
	    $last_sample_id++;	    
	    #and copy to the variation database
	    write_file($new_sample_id,$name,$size,$description);
	}
	else{
	    #this sample is already in the variation database
	    $sanger_sample->{$sample_id}++;
	}
        #and store the relation with the old one
	$old_new_sample_id->{$sample_id} = $new_sample_id;
    }   
    $sth->finish;
    #and finally import the table
    my $call = "lsrcp $LOCAL_TMP_DIR/$TMP_FILE ecs4a:$TMP_DIR/$TMP_FILE";
    system($call);
    load($dbVar->dbc,"sample", "sample_id", "name", "size", "description");
    $call = "$LOCAL_TMP_DIR/$TMP_FILE";
    unlink ($call);   
#    copy_file("sample.txt");
}

#check wether the sample_id is already present in the Variation database as a individual or a population
sub get_sample_variation_database{
    my $dbVariation = shift;
    my $sanger_sample_name = shift;

    my $variation_sample_id = 0;

    my $pop_adaptor = $dbVariation->get_PopulationAdaptor();
    my $ind_adaptor = $dbVariation->get_IndividualAdaptor();

    my $population = $pop_adaptor->fetch_by_name($sanger_sample_name);
    if (defined($population)){
	$variation_sample_id = $population->dbID();
    }
    else{
	my $individual = $ind_adaptor->fetch_by_name($sanger_sample_name);
	if (defined $individual){
	    $variation_sample_id = $individual->dbID();
	}
    }
    return $variation_sample_id;
}

sub import_Source_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_source_id = shift;
    my $last_source_id = shift;

    debug("Loading Source table");
    my ($source_id, $name, $version);
    my $new_source_id;
    my $sth = $dbSanger->dbc->prepare(qq{SELECT source_id, name, version from source});
    $sth->execute();
    $sth->bind_columns(\$source_id, \$name, \$version);
    while ($sth->fetch){
	#get the new id for the source in the variation table
	$new_source_id = $last_source_id + 1;
	#and store the relation with the old one
	$old_new_source_id->{$source_id} = $new_source_id;
	$last_source_id++;
	write_file($new_source_id,$name,$version);
    }   
    $sth->finish;
    #and finally import the table
    my $call = "lsrcp $LOCAL_TMP_DIR/$TMP_FILE ecs4a:$TMP_DIR/$TMP_FILE";
    system($call);
    load($dbVar->dbc,"source","source_id", "name", "version");
    $call = "$LOCAL_TMP_DIR/$TMP_FILE";
    unlink ($call);   
#    copy_file("source.txt");
}

sub import_Population_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_sample_id = shift;
    my $sanger_sample = shift;

    debug("Load Population table");
    my ($sample_id, $is_strain);
    my $sth = $dbSanger->dbc->prepare(qq{SELECT sample_id, is_strain from population});
    $sth->execute();
    $sth->bind_columns(\$sample_id, \$is_strain);
    while ($sth->fetch){
	if (!defined $sanger_sample->{$sample_id}){
	    #get the new id for the sample in the variation table
	    write_file($old_new_sample_id->{$sample_id}, $is_strain);
	}
    }   
    $sth->finish;
    #and finally import the table    
    my $call = "lsrcp $LOCAL_TMP_DIR/$TMP_FILE ecs4a:$TMP_DIR/$TMP_FILE";
    system($call);

    load($dbVar->dbc,"population","sample_id", "is_strain");
    $call = "$LOCAL_TMP_DIR/$TMP_FILE";
    unlink ($call);   
    #copy_file("population.txt");
}

sub import_Individual_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_sample_id = shift;

    debug("Load Individual table");
    my ($sample_id, $father_individual_sample_id, $mother_individual_sample_id, $gender);
    my $sth = $dbSanger->dbc->prepare(qq{SELECT sample_id, father_individual_sample_id, mother_individual_sample_id, gender from individual});
    $sth->execute();
    $sth->bind_columns(\$sample_id, \$father_individual_sample_id, \$mother_individual_sample_id, \$gender);
    while ($sth->fetch){
	#get the new id for the sample in the variation table
	write_file($old_new_sample_id->{$sample_id}, $father_individual_sample_id, $mother_individual_sample_id, $gender);
    }   
    $sth->finish;
    #and finally import the table    
    my $call = "lsrcp $LOCAL_TMP_DIR/$TMP_FILE ecs4a:$TMP_DIR/$TMP_FILE";
    system($call);

    load($dbVar->dbc,"individual","sample_id", "father_individual_sample_id","mother_individual_sample_id","gender");
    $call = "$LOCAL_TMP_DIR/$TMP_FILE";
    unlink ($call);   
    #copy_file("individual.txt");
}

sub import_Individual_Population_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_sample_id = shift;

    debug("Load Individual_Population table");
    my ($individual_sample_id, $population_sample_id);
    my $sth = $dbSanger->dbc->prepare(qq{SELECT individual_sample_id, population_sample_id from individual_population});
    $sth->execute();
    $sth->bind_columns(\$individual_sample_id, \$population_sample_id);
    while ($sth->fetch){
	#get the new id for the sample in the variation table
	write_file($old_new_sample_id->{$individual_sample_id}, $old_new_sample_id->{$population_sample_id});
    }   
    $sth->finish;
    #and finally import the table    
    my $call = "lsrcp $LOCAL_TMP_DIR/$TMP_FILE ecs4a:$TMP_DIR/$TMP_FILE";
    system($call);

    load($dbVar->dbc,"individual_population","individual_sample_id", "population_sample_id");
    $call = "$LOCAL_TMP_DIR/$TMP_FILE";
    unlink ($call);   
    #copy_file("individual_population.txt");
}

sub import_Meta_table{
    my $dbSanger = shift;
    my $dbVar = shift;

    debug("Load Meta table");

    dumpSQL($dbSanger->dbc,qq{SELECT meta_key,meta_value FROM meta});    
    load($dbVar->dbc,"meta","meta_key","meta_value");
}

sub import_Meta_Coord_table{
    my $dbSanger = shift;
    my $dbVar = shift;

    debug("Load MetaCoord table");

    dumpSQL($dbSanger->dbc,qq{SELECT table_name, coord_system_id, max_length FROM meta_coord});    
    load($dbVar->dbc,"meta_coord","table_name","coord_system_id","max_length");
}

sub import_Variation_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_variation_id = shift;
    my $last_variation_id = shift;
    my $old_new_source_id = shift;

    debug("Loading Variation table");
    my ($variation_id, $source_id, $name, $validation_status, $ancestral_allele);
    my $new_variation_id;
    my $sth = $dbSanger->dbc->prepare(qq{SELECT variation_id, source_id, name, validation_status, ancestral_allele from variation});
    $sth->{'mysql_use_result'} = 1;
    $sth->execute();
    $sth->bind_columns(\$variation_id, \$source_id, \$name, \$validation_status, \$ancestral_allele);
    while ($sth->fetch){
	#get the new id for the variation in the variation table
	$new_variation_id = $last_variation_id + 1;
	#and store the relation with the old one
	$old_new_variation_id->{$variation_id} = $new_variation_id;
	$last_variation_id++;
	write_file($new_variation_id,$old_new_source_id->{$source_id}, $name, $validation_status, $ancestral_allele);
    }   
    $sth->finish;
    #and finally import the table
    copy_file("variation.txt");
}

sub import_Allele_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_variation_id = shift;
    my $old_new_sample_id = shift;

    debug("Load Allele table");
    my ($variation_id, $sample_id, $allele, $frequency);
    my $sth = $dbSanger->dbc->prepare(qq{SELECT variation_id, sample_id, allele, frequency from allele});
    $sth->{'mysql_use_result'} = 1;
    $sth->execute();
    $sth->bind_columns(\$variation_id, \$sample_id, \$allele, \$frequency);
    while ($sth->fetch){
	#get the new id for the sample and variation in the variation table
	write_file($old_new_variation_id->{$variation_id}, $old_new_sample_id->{$sample_id}, $allele, $frequency);
    }   
    $sth->finish;
    #and finally import the table    
    copy_file("allele.txt");

}

sub import_Flanking_sequence_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_variation_id = shift;

    debug("Load Flanking Sequence table");
    my ($variation_id, $up_seq, $down_seq, $up_seq_region_start, $up_seq_region_end, $down_seq_region_start, $down_seq_region_end, $seq_region_id, $seq_region_strand);
    my $sth = $dbSanger->dbc->prepare(qq{SELECT variation_id, up_seq, down_seq, up_seq_region_start, up_seq_region_end, down_seq_region_start, down_seq_region_end, seq_region_id, seq_region_strand from flanking_sequence});
    $sth->{'mysql_use_result'} = 1;
    $sth->execute();
    $sth->bind_columns(\$variation_id, \$up_seq, \$down_seq, \$up_seq_region_start, \$up_seq_region_end, \$down_seq_region_start, \$down_seq_region_end,\$seq_region_id, \$seq_region_strand);
    while ($sth->fetch){
	#get the new id for the variation in the variation table
	write_file($old_new_variation_id->{$variation_id},$up_seq, $down_seq, $up_seq_region_start, $up_seq_region_end, $down_seq_region_start, $down_seq_region_end, $seq_region_id, $seq_region_strand);
    }   
    $sth->finish;
    #and finally import the table    
    copy_file("flanking_sequence.txt");

}

sub import_Variation_synonym_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_variation_id = shift;
    my $old_new_source_id = shift;

    debug("Load Variation synonym table");
    my ($variation_id, $source_id, $name, $moltype);
    my $sth = $dbSanger->dbc->prepare(qq{SELECT variation_id, source_id, name, moltype from variation_synonym});
    $sth->{'mysql_use_result'} = 1;
    $sth->execute();
    $sth->bind_columns(\$variation_id, \$source_id, \$name, \$moltype);
    while ($sth->fetch){
	#get the new id for the source and variation in the variation table
	write_file($old_new_variation_id->{$variation_id}, $old_new_source_id->{$source_id}, $name, $moltype);
    }   
    $sth->finish;
    #and finally import the table    
    copy_file("variation_synonym.txt");
}

sub import_Variation_feature_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_variation_feature_id = shift;
    my $last_variation_feature_id = shift;
    my $old_new_variation_id = shift;
    my $old_new_source_id = shift;

    debug("Loading VariationFeature table");
    my ($variation_feature_id, $variation_id, $seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand, $allele_string, $variation_name, $map_weight, $flags, $source_id, $validation_status, $consequence_type);
    my $new_variation_feature_id;
    my $sth = $dbSanger->dbc->prepare(qq{SELECT variation_feature_id, variation_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, allele_string, variation_name, map_weight, flags, source_id, validation_status, consequence_type from variation_feature});
    $sth->{'mysql_use_result'} = 1;
    $sth->execute();
    $sth->bind_columns(\$variation_feature_id, \$variation_id, \$seq_region_id, \$seq_region_start, \$seq_region_end, \$seq_region_strand, \$allele_string, \$variation_name, \$map_weight, \$flags, \$source_id, \$validation_status, \$consequence_type);
    while ($sth->fetch){
	#get the new id for the variation_feature in the variation table
	$new_variation_feature_id = $last_variation_feature_id + 1;
	#and store the relation with the old one
	$old_new_variation_feature_id->{$variation_feature_id} = $new_variation_feature_id;
	$last_variation_feature_id++;
	write_file($new_variation_feature_id,$old_new_variation_id->{$variation_id}, $seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand, $allele_string, $variation_name, $map_weight, $flags, $old_new_source_id->{$source_id}, $validation_status, $consequence_type);
    }   
    $sth->finish;
    #and finally import the table
    copy_file("variation_feature.txt");
}

sub import_Transcript_variation_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_variation_feature_id = shift;

    debug("Load TranscriptVariation table");
    my ($transcript_id, $variation_feature_id, $cdna_start, $cdna_end, $translation_start, $translation_end, $peptide_allele_string, $consequence_type);
    my $sth = $dbSanger->dbc->prepare(qq{SELECT transcript_id, variation_feature_id, cdna_start, cdna_end, translation_start, translation_end, peptide_allele_string, consequence_type from transcript_variation});
    $sth->{'mysql_use_result'} = 1;
    $sth->execute();
    $sth->bind_columns(\$transcript_id, \$variation_feature_id, \$cdna_start, \$cdna_end, \$translation_start, \$translation_end, \$peptide_allele_string, \$consequence_type);
    while ($sth->fetch){
	#get the new id for the variation_feature in the variation table
	write_file($transcript_id, $old_new_variation_feature_id->{$variation_feature_id}, $cdna_start, $cdna_end, $translation_start, $translation_end, $peptide_allele_string, $consequence_type);
    }   
    $sth->finish;
    #and finally import the table    
    copy_file("transcript_variation.txt");

}

sub import_Read_coverage_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_sample_id = shift;

    debug("Load ReadCoverage table");
    my ($seq_region_id, $seq_region_start, $seq_region_end, $level, $sample_id);
    my $sth = $dbSanger->dbc->prepare(qq{SELECT seq_region_id, seq_region_start, seq_region_end, level, sample_id from read_coverage});
    $sth->{'mysql_use_result'} = 1;
    $sth->execute();
    $sth->bind_columns(\$seq_region_id, \$seq_region_start, \$seq_region_end, \$level, \$sample_id);
    while ($sth->fetch){
	#get the new id for the sample in the variation table
	write_file($seq_region_id,$seq_region_start, $seq_region_end, $level, $old_new_sample_id->{$sample_id});
    }   
    $sth->finish;
    #and finally import the table    
    copy_file("read_coverage.txt");
}

sub import_Tmp_individual_genotype_single_bp_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_variation_id = shift;
    my $old_new_sample_id = shift;

    debug("Load tmp_individual_genotype_single_bp table");
    my ($variation_id, $sample_id, $allele_1, $allele_2);
    my $sth = $dbSanger->dbc->prepare(qq{SELECT variation_id, sample_id, allele_1, allele_2 from tmp_individual_genotype_single_bp});
    $sth->{'mysql_use_result'} = 1;
    $sth->execute();
    $sth->bind_columns(\$variation_id, \$sample_id, \$allele_1, \$allele_2);
    while ($sth->fetch){
	#get the new id for the sample and variation in the variation table
	write_file($old_new_variation_id->{$variation_id}, $old_new_sample_id->{$sample_id}, $allele_1, $allele_2);
    }   
    $sth->finish;
    #and finally import the table    
    copy_file("tmp_individual_genotype_single_bp.txt");

}

sub write_file{
    my @values = @_;

    open FH, ">>$LOCAL_TMP_DIR/$TMP_FILE" || die "Could not open file with information: $!\n";
    my @a = map {defined($_) ? $_ : '\N'} @values; #to replace undefined values by \N in the file
    print FH join("\t", @a), "\n";
    close FH || die "Could not close file with information: $!\n";
    
}

#function to return the last id used in the table tablename (the id must be called "tablename_id")
sub get_last_table_id{
    my $dbSanger = shift;
    my $tablename = shift;

    my $max_id;
    my $sth = $dbSanger->dbc->prepare(qq{SELECT MAX($tablename\_id) from $tablename});
    $sth->execute();
    $sth->bind_columns(\$max_id);
    $sth->fetch;
    $sth->finish();

    return $max_id if (defined $max_id);
    return 0 if (!defined $max_id);
}

#to cpy the file from bc node to ecs4
sub copy_file{
    my $filename = shift;

    my $call = "lsrcp $LOCAL_TMP_DIR/$TMP_FILE ecs4a:$TMP_DIR/$filename";
    system($call);
    #and remove the file
    $call = "$LOCAL_TMP_DIR/$TMP_FILE";
    unlink ($call);   
}
