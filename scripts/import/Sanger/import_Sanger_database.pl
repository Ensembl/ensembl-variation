#!/usr/bin/env perl
# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
my ($species,$redo, $add_new_name);

GetOptions('tmpdir=s'  => \$ImportUtils::TMP_DIR,
	   'tmpfile=s' => \$ImportUtils::TMP_FILE,
	   'species=s' => \$species,
	   'redo'      => \$redo, #if set -redo, will not load tables like sample/population/individual etc again
	   'add_new_name' => \$add_new_name #if add ENS name after current one
	   );

warn("Make sure you have a updated ensembl.registry file!\n");
die "you must specify the species:" if (!$species);
my $registry_file ||= $Bin . "/ensembl.registry";


Bio::EnsEMBL::Registry->load_all( $registry_file );

my $dbVar = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
my $dbSanger = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'sanger');

my $TMP_DIR  = $ImportUtils::TMP_DIR;
#my $LOCAL_TMP_DIR  = '/tmp';
my $TMP_FILE = $ImportUtils::TMP_FILE;
my $buffer = {};



my $old_new_variation_id = {}; #reference to a hash with the old_variation_id => new_variation_id
my $old_new_variation_feature_id = {}; #reference to a hash with the old_Variation_feature_id => new_variation_feature_id
my $old_new_source_id = {}; # reference to a hash with the old_source_id => new_source_id
my $old_new_sample_id = {}; #reference to a hash with the old_sample_id => new_sample_id
my $sanger_sample = {}; #reference to a hash containing the samples that are in the variation database
my $old_new_variation_name = {}; #reference to a hash with the old_variation_name => new_variation_name

my $last_source_id = get_last_table_id($dbVar,"source"); #last source_id used in the database
my $last_sample_id = get_last_table_id($dbVar,"sample"); #last sample_id used in the database
my $last_variation_id = get_last_table_id($dbVar,"variation"); #last variation_id used in the database
my $last_variation_feature_id = get_last_table_id($dbVar,"variation_feature"); #last variation_feature_id used in the database

load_existing_stable_id($dbSanger,$dbVar) if ($add_new_name);
import_Sample_table($dbSanger, $dbVar, $old_new_sample_id, $last_sample_id, $sanger_sample);
import_Source_table($dbSanger, $dbVar, $old_new_source_id, $last_source_id);
import_Population_table($dbSanger, $dbVar, $old_new_sample_id, $sanger_sample);
import_Individual_table($dbSanger,$dbVar,$old_new_sample_id);
import_Individual_Population_table($dbSanger,$dbVar,$old_new_sample_id);
import_Meta_table($dbSanger,$dbVar);
import_Meta_Coord_table($dbSanger,$dbVar);
import_Variation_table($dbSanger,$dbVar,$old_new_variation_id,$last_variation_id, $old_new_source_id, $old_new_variation_name);
import_Allele_table($dbSanger,$dbVar,$old_new_variation_id, $old_new_sample_id);
import_Flanking_sequence_table($dbSanger,$dbVar,$old_new_variation_id);
#import_Variation_synonym_table($dbSanger, $dbVar, $old_new_variation_id, $old_new_source_id);
#import_Failed_variation_table($dbSanger, $dbVar, $old_new_variation_id);
import_Variation_feature_table($dbSanger,$dbVar,$old_new_variation_feature_id,$last_variation_feature_id, $old_new_variation_id, $old_new_source_id, $old_new_variation_name);
#import_Transcript_variation_table($dbSanger,$dbVar,$old_new_variation_feature_id);
#import_Read_coverage_table($dbSanger,$dbVar, $old_new_sample_id);
#import_Tmp_individual_genotype_single_bp_table($dbSanger,$dbVar,$old_new_variation_id,$old_new_sample_id);

sub load_existing_stable_id {
  my $dbSanger = shift;
  my $dbVar = shift;

  my $dbvarname = $dbVar->dbc->dbname();

  #$dbSanger->dbc->do(qq{ALTER TABLE variation add column (exist_var_flag int default 0 not null,exist_var_name varchar(50))});
  $dbSanger->dbc->do(qq{DROP TABLE IF EXISTS var_stableid_old_new});
  $dbSanger->dbc->do(qq{CREATE TABLE var_stableid_old_new (exist_name varchar(50),new_name varchar(50),unique(new_name))});
  my $source_id_ref = $dbVar->dbc->db_handle->selectall_arrayref(qq{SELECT source_id from source where name like "ENSEMBL%"});
  my @source_ids = map {$_->[0]} @{$source_id_ref};
  #print "source_ids is @source_ids\n";

  foreach my $source_id (@source_ids) {
    debug("Processing source_id $source_id...");
    $dbSanger->dbc->do(qq{INSERT IGNORE INTO var_stableid_old_new
                          SELECT vs.name as exist_name, vf1.variation_name as new_name 
                          FROM variation_feature vf1, $dbvarname.variation_feature vf2, $dbvarname.variation_synonym vs 
                          WHERE vf1.seq_region_id=vf2.seq_region_id 
                          AND vf1.seq_region_start=vf2.seq_region_start 
                          AND vf1.seq_region_end=vf2.seq_region_end 
                          AND vs.source_id = $source_id 
                          AND vf2.variation_id=vs.variation_id});
    $dbSanger->dbc->do(qq{INSERT IGNORE INTO var_stableid_old_new
                          SELECT vf2.variation_name as exist_name, vf1.variation_name as new_name
                          FROM variation_feature vf1, $dbvarname.variation_feature vf2
                          WHERE vf1.seq_region_id=vf2.seq_region_id 
                          AND vf1.seq_region_start=vf2.seq_region_start 
                          AND vf1.seq_region_end=vf2.seq_region_end 
                          AND vf2.source_id = $source_id });
  }

  $dbSanger->dbc->do(qq{UPDATE variation v, var_stableid_old_new o
                        SET v.exist_var_name = o.exist_name,v.exist_var_flag = 1
                        WHERE v.name = o.new_name});

}

debug("Finished reading data, ready to start loading in the database....");


sub import_Sample_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_sample_id = shift;
    my $last_sample_id = shift;
    my $sanger_sample = shift;

    debug("Loading Sample table");
    my ($sample_id, $name, $size, $description);

    my $sth = $dbSanger->dbc->db_handle->prepare(qq{SELECT sample_id, name, size, description from sample});
    $sth->execute();
    $sth->bind_columns(\$sample_id, \$name, \$size, \$description);
    while ($sth->fetch){
      my $new_sample_id;
      #need to check if the strain is already in the Variation database, as a population or individual
      if ($description =~ /individual/i) {
	$new_sample_id = &get_sample_variation_database($dbVar, $name,"ind");
      }
      elsif ($description =~ /population/i) {
	$new_sample_id = &get_sample_variation_database($dbVar, $name,"pop");
	#print " sample_id is $sample_id and new_sample_id = $new_sample_id\n" if $new_sample_id;
      }
      if (! $new_sample_id){
	#get the new id for the sample in the variation table
	$new_sample_id = $last_sample_id + 1;
	$last_sample_id++;
	#and write to tmpdir/tmpfile
	write_file($new_sample_id,$name,$size,$description);
      }
      else{
	#this sample is already in the variation database
	$sanger_sample->{$sample_id}++;
      }
      #and store the relation with the old one
      #print "sample_id is $sample_id and new_sample_id is $new_sample_id\n";
      $old_new_sample_id->{$sample_id} = $new_sample_id;
    }
    #print Dumper($old_new_sample_id);
    $sth->finish;
    #and finally import the table
    load($dbVar->dbc,"sample", "sample_id", "name", "size", "description") if (!$redo);
    my $call = "$TMP_DIR/$TMP_FILE";
    unlink ($call);   
}

#check wether the sample_id is already present in the Variation database as a individual or a population
sub get_sample_variation_database{
    my $dbVariation = shift;
    my $sanger_sample_name = shift;
    my $sample_type = shift;
    my $population_sample_id = 0;
    my $individual_sample_id = 0;

    my $pop_adaptor = $dbVariation->get_PopulationAdaptor();
    my $ind_adaptor = $dbVariation->get_IndividualAdaptor();

    my $population = $pop_adaptor->fetch_by_name($sanger_sample_name);
    if (defined($population) and $sample_type eq "pop"){
	$population_sample_id = $population->dbID();
	return $population_sample_id;
    }
    else{
	my $individual = $ind_adaptor->fetch_all_by_name($sanger_sample_name);
	if (defined $individual and $$individual[0] and $sample_type eq "ind") {
	    $individual_sample_id = $$individual[0]->dbID();
	    return $individual_sample_id;
	}
    }
}

sub import_Source_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_source_id = shift;
    my $last_source_id = shift;

    debug("Loading Source table");
    my ($source_id, $name, $version, $description);
    my $new_source_id;
    my $sth = $dbSanger->dbc->db_handle->prepare(qq{SELECT source_id, name, version, description from source});
    $sth->execute();
    $sth->bind_columns(\$source_id, \$name, \$version, \$description);
    while ($sth->fetch){
      #need to check if the source is already in the Variation database
      $new_source_id = get_source_variation_database($dbVar,$name);
      if (! $new_source_id) {
	#get the new id for the source in the variation table
	$new_source_id = $last_source_id + 1;
	$last_source_id++;
	write_file($new_source_id,$name,$version, $description);
	#and finally import the table
	load($dbVar->dbc,"source","source_id", "name", "version", "description") if !$redo;
	my $call = "$TMP_DIR/$TMP_FILE";
	unlink ($call);  
      }
      #and store the relation with the old one
      $old_new_source_id->{$source_id} = $new_source_id;
    }   
    $sth->finish;
}

#check whether the source_id is already present in the Variation database
sub get_source_variation_database{
    my $dbVariation = shift;
    my $sanger_source_name = shift;
    my $source_id = 0;

    my $source_id_ref = $dbVariation->dbc->db_handle->selectall_arrayref(qq{SELECT source_id from source 
                    where name = "$sanger_source_name"});
    $source_id = $source_id_ref->[0][0] if $source_id_ref;
    return $source_id;
}

sub import_Population_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_sample_id = shift;
    my $sanger_sample = shift;

    debug("Load Population table");
    my ($sample_id);
    my $sth = $dbSanger->dbc->db_handle->prepare(qq{SELECT sample_id from population});
    $sth->execute();
    $sth->bind_columns(\$sample_id);
    while ($sth->fetch){
	if (!defined $sanger_sample->{$sample_id}){
	    #get the new id for the sample in the variation table
	    write_file($old_new_sample_id->{$sample_id});
	}
    }   
    $sth->finish;
    #and finally import the table    
    load($dbVar->dbc,"population","sample_id") if !$redo;
    my $call = "$TMP_DIR/$TMP_FILE";
    unlink ($call);   
}

sub import_Individual_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_sample_id = shift;

    debug("Load Individual table");
    my ($sample_id, $father_individual_sample_id, $mother_individual_sample_id, $gender, $individual_type_id);
    my $sth = $dbSanger->dbc->db_handle->prepare(qq{SELECT sample_id, father_individual_sample_id, mother_individual_sample_id, gender, individual_type_id from individual});
    $sth->execute();
    $sth->bind_columns(\$sample_id, \$father_individual_sample_id, \$mother_individual_sample_id, \$gender,\$individual_type_id);
    while ($sth->fetch){
      if (! defined $sanger_sample->{$sample_id}) {
	#get the new id for the sample in the variation table
	write_file($old_new_sample_id->{$sample_id}, $father_individual_sample_id, $mother_individual_sample_id, $gender, $individual_type_id);
      }
    }
    $sth->finish;
    #and finally import the table    
    load($dbVar->dbc,"individual","sample_id", "father_individual_sample_id","mother_individual_sample_id","gender","individual_type_id") if !$redo;
    my $call = "$TMP_DIR/$TMP_FILE";
    unlink ($call);   

}

sub import_Individual_Population_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_sample_id = shift;

    debug("Load Individual_Population table");
    my ($individual_sample_id, $population_sample_id);
    my $sth = $dbSanger->dbc->db_handle->prepare(qq{SELECT individual_sample_id, population_sample_id from individual_population});
    $sth->execute();
    $sth->bind_columns(\$individual_sample_id, \$population_sample_id);
    while ($sth->fetch){
	#get the new id for the sample in the variation table
	write_file($old_new_sample_id->{$individual_sample_id}, $old_new_sample_id->{$population_sample_id});
    }   
    $sth->finish;
    #and finally import the table    
    load($dbVar->dbc,"individual_population","individual_sample_id", "population_sample_id") if !$redo ;
    my $call = "$TMP_DIR/$TMP_FILE";
    unlink ($call);   
}

sub import_Meta_table{
    my $dbSanger = shift;
    my $dbVar = shift;

    debug("Load Meta table");

    dumpSQL($dbSanger->dbc,qq{SELECT meta_key,meta_value FROM meta});    
    load($dbVar->dbc,"meta","meta_key","meta_value") if !$redo;
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
  my $old_new_variation_name = shift;
  my $stable_id_num =0;
  my $species_pre_name_length;

  if ($species =~ /hum/i) {
    $species_pre_name_length = 7;
  }
  else {
    $species_pre_name_length = 10;
  }

  if ($add_new_name) {
    #get last ensembl_stable name
    my $stable_id_ref = $dbVar->dbc->db_handle->selectall_arrayref(qq{select if(t.number > r.number,t.number,r.number) from (select max(round(substring(v.name,$species_pre_name_length))) as number from variation v where v.name like 'ENS%') as t, (select max(round(substring(vs.name,$species_pre_name_length))) as number from variation_synonym vs where vs.name like 'ENS%') as r});

    $stable_id_num = $stable_id_ref->[0][0] if $stable_id_ref;
    $stable_id_num =~ s/^0+// if $stable_id_num > 0;
  }

  print "the biggest stable_id_num is $stable_id_num\n";
  debug("Loading Variation table");

  my $new_variation_id;
  my $new_stable_id_num = $stable_id_num;

  dumpSQL($dbSanger->dbc->db_handle,qq{SELECT * FROM variation});
  system("mv $TMP_DIR/$TMP_FILE $TMP_DIR/$TMP_FILE\_variation");
  open IN, "$TMP_DIR/$TMP_FILE\_variation" or die "Can't open imput file\n";
  while (<IN>){
    my ($variation_id, $source_id, $name, $validation_status, $ancestral_allele,$exist_var_flag,$exist_var_name) = split;
    #get the new id for the variation in the variation table
    $new_variation_id = $last_variation_id + 1;
    #and store the relation with the old one
    $old_new_variation_id->{$variation_id} = $new_variation_id;
    $last_variation_id++;
    my ($pre_name,$num,$new_name);
    #if ensembl stable_id is already exist, new stable_id = old_stable_id + 1
    if ($add_new_name) {
      if ($name =~ /(^ENS.*SNP)(\d+)/ and $stable_id_num != 0) {
	$pre_name = $1;
	$num = $2;
	if (!$exist_var_flag) {##this is original data, i.e not overlap with existing ensembl_stable_ids
	  $new_stable_id_num++ ;
	  $new_name = "$pre_name$new_stable_id_num";
	}
	else {
	  $new_name = $exist_var_name;####this is stable_ids that overlap with existing ensembl stable_ids
	}
      }
    }
    else {
      $new_name = $name;####this is stable_ids that already existing
    }
    $old_new_variation_name->{$name} = $new_name;
    #write_file($new_variation_id,$old_new_source_id->{$source_id}, $new_name, $validation_status, $ancestral_allele);

    print_buffered($buffer,"$TMP_DIR/variation.txt",join ("\t", $new_variation_id,$old_new_source_id->{$source_id}, $new_name, $validation_status, $ancestral_allele)."\n");
  }
  close IN;
  print_buffered($buffer);

  #copy_file("variation.txt") ;
  unlink "$TMP_DIR/$TMP_FILE\_variation";

}

sub import_Allele_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_variation_id = shift;
    my $old_new_sample_id = shift;
    my $count = 0;
    my (@a,@lines);

    debug("Load Allele table");

    dumpSQL($dbSanger->dbc->db_handle,qq{SELECT variation_id, allele from allele});
    system("mv $TMP_DIR/$TMP_FILE $TMP_DIR/$TMP_FILE\_allele");
    open IN, "$TMP_DIR/$TMP_FILE\_allele" or die "Can't open imput file\n";
    while (<IN>){
      my ($variation_id, $allele) = split;
      my $frequency ||='\N';
      my $sample_id ||='\N';
      #get the new id for the sample and variation in the variation table
      #write_file($old_new_variation_id->{$variation_id}, $allele, $frequency, $old_new_sample_id->{$sample_id});

      my $new_sample_id = $old_new_sample_id->{$sample_id};
      $new_sample_id |='\N';
      print_buffered($buffer,"$TMP_DIR/allele.txt", join ("\t",$old_new_variation_id->{$variation_id}, $allele, $frequency, $new_sample_id)."\n");
    }
    close IN;
    print_buffered($buffer);
    #and finally import the table    
#    copy_file("allele.txt");
    unlink "$TMP_DIR/$TMP_FILE\_allele";
}

sub import_Flanking_sequence_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_variation_id = shift;

    debug("Load Flanking Sequence table");

    dumpSQL($dbSanger->dbc->db_handle,qq{SELECT variation_id, up_seq, down_seq, up_seq_region_start, up_seq_region_end, down_seq_region_start, down_seq_region_end, seq_region_id, seq_region_strand from flanking_sequence});

    system("mv $TMP_DIR/$TMP_FILE $TMP_DIR/$TMP_FILE\_flank");
    open IN, "$TMP_DIR/$TMP_FILE\_flank" or die "Can't open imput file\n";
    while (<IN>){
      my ($variation_id, $up_seq, $down_seq, $up_seq_region_start, $up_seq_region_end, $down_seq_region_start, $down_seq_region_end, $seq_region_id, $seq_region_strand) = split;
      #get the new id for the variation in the variation table
      #write_file($old_new_variation_id->{$variation_id},$up_seq, $down_seq, $up_seq_region_start, $up_seq_region_end, $down_seq_region_start, $down_seq_region_end, $seq_region_id, $seq_region_strand);
      $up_seq ||='\N';
      $down_seq ||='\N';

      print_buffered($buffer,"$TMP_DIR/flanking_sequence.txt",join ("\t",$old_new_variation_id->{$variation_id},$up_seq, $down_seq, $up_seq_region_start, $up_seq_region_end, $down_seq_region_start, $down_seq_region_end, $seq_region_id, $seq_region_strand)."\n");
    }
    close IN;
    print_buffered($buffer);
    #and finally import the table    
#    copy_file("flanking_sequence.txt");
    unlink "$TMP_DIR/$TMP_FILE\_flank";

}

sub import_Variation_synonym_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_variation_id = shift;
    my $old_new_source_id = shift;

    debug("Load Variation synonym table");

    dumpSQL($dbSanger->dbc->db_handle,(qq{SELECT variation_id, source_id, name, moltype from variation_synonym}));
    system("mv $TMP_DIR/$TMP_FILE $TMP_DIR/$TMP_FILE\_synonym");
    open IN, "$TMP_DIR/$TMP_FILE\_synonym" or die "Can't open imput file\n";
    while (<IN>){
      my ($variation_id, $source_id, $name, $moltype) = split;
      #get the new id for the source and variation in the variation table
      #write_file($old_new_variation_id->{$variation_id}, $old_new_source_id->{$source_id}, $name, $moltype);
      print_buffered($buffer,"$TMP_DIR/variation_synonym.txt",join ("\t",$old_new_variation_id->{$variation_id}, $old_new_source_id->{$source_id}, $name, $moltype)."\n");
    }

    #and finally import the table    
    #copy_file("variation_synonym.txt");
    print_buffered($buffer);
    unlink "$TMP_DIR/$TMP_FILE\_synonym";
}

sub import_Failed_variation_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_variation_id = shift;

    debug("Load Failed variation table");

    dumpSQL($dbSanger->dbc->db_handle,(qq{SELECT variation_id, failed_description_id from failed_variation}));
    system("mv $TMP_DIR/$TMP_FILE $TMP_DIR/$TMP_FILE\_failed");
    open IN, "$TMP_DIR/$TMP_FILE\_failed" or die "Can't open imput file\n";
    while (<IN>){
      my ($variation_id, $failed_description_id) = split;
      #get the new id for the source and variation in the variation table
      #write_file($old_new_variation_id->{$variation_id}, $old_new_source_id->{$source_id}, $name, $moltype);
      print_buffered($buffer,"$TMP_DIR/failed_variation.txt",join ("\t",$old_new_variation_id->{$variation_id}, $failed_description_id)."\n");
    }

    #and finally import the table    
    #copy_file("failed_variation.txt");
    print_buffered($buffer);
    unlink "$TMP_DIR/$TMP_FILE\_failed";
}

sub import_Variation_feature_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_variation_feature_id = shift;
    my $last_variation_feature_id = shift;
    my $old_new_variation_id = shift;
    my $old_new_source_id = shift;
    my $old_new_variation_name = shift;

    debug("Loading VariationFeature table");

    dumpSQL($dbSanger->dbc->db_handle,qq{SELECT variation_feature_id, variation_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, allele_string, variation_name, map_weight, flags, source_id, validation_status, consequence_type from variation_feature});

    system("mv $TMP_DIR/$TMP_FILE $TMP_DIR/$TMP_FILE\_feature");
    open IN, "$TMP_DIR/$TMP_FILE\_feature" or die "Can't open imput file\n";
    while (<IN>){
      my ($variation_feature_id, $variation_id, $seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand, $allele_string, $variation_name, $map_weight, $flags, $source_id, $validation_status, $consequence_type) = split;

	#get the new id for the variation_feature in the variation table
	my $new_variation_feature_id = $last_variation_feature_id + 1;
	#and store the relation with the old one
	$old_new_variation_feature_id->{$variation_feature_id} = $new_variation_feature_id;
	$last_variation_feature_id++;
	my $new_variation_name = $old_new_variation_name->{$variation_name};
	$variation_name = $new_variation_name if $new_variation_name;
        $validation_status ||='\N';
	#write_file($new_variation_feature_id,$seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand, $old_new_variation_id->{$variation_id}, $allele_string, $variation_name, $map_weight, $flags, $old_new_source_id->{$source_id}, $validation_status, $consequence_type);

	print_buffered($buffer,"$TMP_DIR/variation_feature.txt",join ("\t",$new_variation_feature_id,$seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand, $old_new_variation_id->{$variation_id}, $allele_string, $variation_name, $map_weight, $flags, $old_new_source_id->{$source_id}, $validation_status, $consequence_type)."\n");

    }   
    close IN;
    print_buffered($buffer);
    #and finally import the table
#    copy_file("variation_feature.txt");
    unlink "$TMP_DIR/$TMP_FILE\_feature";
}

sub import_Transcript_variation_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_variation_feature_id = shift;

    debug("Load TranscriptVariation table");
    my ($transcript_id, $variation_feature_id, $cdna_start, $cdna_end, $translation_start, $translation_end, $peptide_allele_string, $consequence_type);
    my $sth = $dbSanger->dbc->db_handle->prepare(qq{SELECT transcript_id, variation_feature_id, cdna_start, cdna_end, translation_start, translation_end, peptide_allele_string, consequence_type from transcript_variation});
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

    dumpSQL($dbSanger->dbc->db_handle,qq{SELECT seq_region_id, seq_region_start, seq_region_end, level, sample_id from read_coverage});

    system("mv $TMP_DIR/$TMP_FILE $TMP_DIR/$TMP_FILE\_read");
    open IN, "$TMP_DIR/$TMP_FILE\_read" or die "Can't open imput file\n";
    while (<IN>){
      my ($seq_region_id, $seq_region_start, $seq_region_end, $level, $sample_id) = split;
      #get the new id for the sample in the variation table
      #write_file($seq_region_id,$seq_region_start, $seq_region_end, $level, $old_new_sample_id->{$sample_id});
      print_buffered($buffer,"$TMP_DIR/read_coverage.txt",join ("\t",$seq_region_id,$seq_region_start, $seq_region_end, $level, $old_new_sample_id->{$sample_id})."\n");
    }   

    print_buffered($buffer);
    #and finally import the table    
    #copy_file("read_coverage.txt");
    unlink "$TMP_DIR/$TMP_FILE\_read";
}

sub import_Tmp_individual_genotype_single_bp_table{
    my $dbSanger = shift;
    my $dbVar = shift;
    my $old_new_variation_id = shift;
    my $old_new_sample_id = shift;

    foreach my $table ("tmp_individual_genotype_single_bp","individual_genotype_multiple_bp") {
      debug("Load $table table");

      dumpSQL($dbSanger->dbc->db_handle,qq{SELECT variation_id, allele_1, allele_2, sample_id from $table});

      system("mv $TMP_DIR/$TMP_FILE $TMP_DIR/$TMP_FILE\_$table");
      open IN, "$TMP_DIR/$TMP_FILE\_$table" or die "Can't open imput file\n";
      while (<IN>){
	my ($variation_id, $allele_1, $allele_2, $sample_id) = split;
	#get the new id for the sample and variation in the variation table
	#if ($old_new_variation_id->{$variation_id} >2811182) {
	#write_file($old_new_variation_id->{$variation_id}, $allele_1, $allele_2, $old_new_sample_id->{$sample_id});

	print_buffered($buffer,"$TMP_DIR/$table\.txt",join ("\t",$old_new_variation_id->{$variation_id}, $allele_1, $allele_2, $old_new_sample_id->{$sample_id})."\n");
      }

      close IN;
      print_buffered($buffer);
      #and finally import the table    
      #copy_file("tmp_individual_genotype_single_bp.txt");
      unlink  "$TMP_DIR/$TMP_FILE\_$table";
    }
}

sub write_file{
    my @values = @_;

    open FH, ">>$TMP_DIR/$TMP_FILE" || die "Could not open file with information: $!\n";
    my @a = map {defined($_) ? $_ : '\N'} @values; #to replace undefined values by \N in the file
    print FH join("\t", @a), "\n";
    close FH || die "Could not close file with information: $!\n";
    
}

#function to return the last id used in the table tablename (the id must be called "tablename_id")
sub get_last_table_id{
    my $dbSanger = shift;
    my $tablename = shift;

    my $max_id;
    my $sth = $dbSanger->dbc->db_handle->prepare(qq{SELECT MAX($tablename\_id) from $tablename});
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

    my $call = "mv $TMP_DIR/$TMP_FILE $TMP_DIR/$filename";
    system($call);
}

sub print_buffered {
    my $buffer = shift;
    my $filename = shift;
    my $text = shift;
    #print "file name is $filename\n";
    local *FH;

    if( ! $filename ) {
	# flush the buffer
	foreach my $file (keys %{$buffer}){
	    open( FH, ">>$file" ) or die "Could not print to file $! \n";
	    print FH $buffer->{ $file };
	    close FH;
	}
	%{$buffer} = (); #flush buffer
    } else {
	$buffer->{ $filename } .= $text;
	if( length( $buffer->{ $filename } ) > 10_000 ) {
	    open( FH, ">>$filename" ) or die;
	    print FH $buffer->{ $filename };
	    close FH;
	    $buffer->{ $filename } = '';
	}
    }
}

#for three individuals for human, import one, then do second one
#for most tables, can use mysqlimport to load the data, but for allele and transcript_variation tables, need:
#mysqlimport -uensadmin -pensembl -hens-genomics2 -L -c "variation_id,allele,frequency,sample_id" yuan_rattus_norvegicus_variation_42_34l allele.txt
#mysqlimport -uensadmin -pensembl -hens-genomics2 -L -c "transcript_id,variation_feature_id,cdna_start,cdna_end,translation_start,translation_end,peptide_allele_string,consequence_type" yuan_rattus_norvegicus_variation_42_34l transcript_variation.txt

