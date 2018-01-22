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

use Bio::EnsEMBL::Registry;
use Bio::Index::Fastq;
use FindBin qw( $Bin );
use Getopt::Long;
use Data::Dumper;
use Fcntl ':flock';
use DBI;
use DBH;
use Time::HiRes qw(tv_interval gettimeofday);
use Bio::EnsEMBL::Utils::Cache;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);
use Bio::EnsEMBL::Utils::Sequence qw(expand reverse_comp);
use ImportUtils qw(debug load create_and_load dumpSQL);

our ($TMP_DIR, $TMP_FILE, $check, $SEQ_REGION_ID, $SEQ_REGION_NAME, $species, $VAR_DBNAME);



GetOptions('species=s'   => \$species,
           'check=s'   => \$check,#such as check1,just to make file name different
	   'tmpdir=s'  => \$ImportUtils::TMP_DIR,
	   'tmpfile=s' => \$ImportUtils::TMP_FILE,
	  );

my $registry_file ||= $Bin . "/ensembl.registry";

Bio::EnsEMBL::Registry->load_all( $registry_file );

my $cdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $vdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');

my $dbVar = $vdba->dbc;
my $dbCore = $cdba;
my $vdbname = $dbVar->dbname;

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

#find out seq_region_id for haplotype chromosomes    
my $sr_ref = $dbCore->dbc->db_handle->selectall_arrayref(qq{select sra.seq_region_id from seq_region_attrib sra, attrib_type at where sra.attrib_type_id=at.attrib_type_id and at.name="Non Reference"});

my %hap_seq_id = map {$_->[0],1} @$sr_ref;
my $hap_id_string = "(". join ",",keys %hap_seq_id, . ")";

check_allele_table($dbVar,$check);
check_genotype_table($dbVar,$check);
#check_variation_feature_table($dbVar);
#change_allele_table($dbVar);
#change_allele_table1($dbVar);
#change_allele_table2($dbVar);
#change_allele_table3($dbVar);
#delete_allele_gtype_tables($dbVar);
#change_allele_gtype_tables($dbVar);
#change_variation_feature_table($dbVar);


sub check_allele_table {
  
  my $dbVar = shift;
  my $check = shift;

  debug("Checking for allele table...One variation should have at least two alleles");
  open OUT, ">$TMP_DIR/$vdbname\_allele_$check";
  #open OUT1,">$TMP_DIR/$vdbname\_allele_check_one_allele";

  #do not use this command, it's too big
  #$dbVar->do(qq{create table varid_one_allele select variation_id, count(*) as count from (select variation_id,allele from allele group by variation_id,allele) as tmp_table group by variation_id having count =1});
  #my $varid_ref = $dbVar->db_handle->selectall_arrayref(qq{select count(*) from varid_one_allele});
  #my $count = $varid_ref->[0][0];
  #print "There is $count record in table varid_one_allele\n";

  debug("Checking for allele table...All alleles in allele_string of variation_feature table should match alleles in allele table for the variation"); 

  my $hap_line;
  if ($hap_id_string =~ /\(\d+\)/) {
    $hap_line = "AND seq_region_id not in $hap_id_string";
  }
  else {
    $hap_line ='';
  }

  #only find out variation with one mapping that is vf.map_weight=1 and seq_region_id not in (226063,226064,226065,226066)
  my $sth = $dbVar->prepare(qq{select a.*,vf.allele_string from allele a, variation_feature vf where a.variation_id = vf.variation_id and vf.map_weight=1 $hap_line}, {mysql_use_result=>1});
  
  my (%total_var,%matched_var,$variation_id,$allele,$allele_string,$allele_id,$subsnp_id,$frequency,$sample_id);
  
  $sth->execute();
  $sth->bind_columns(\$allele_id,\$variation_id,\$subsnp_id,\$allele,\$frequency,\$sample_id,\$allele_string);
  while($sth->fetch()) {
     #print "$variation_id,$allele,$allele_string\n";
     next if ($allele =~ /N|^\(.*\)$|\+/ig);
     next if ($allele =~ /\+/);
     #print "again $variation_id,$allele,$allele_string again\n";
     #$total_var{$variation_id}=1;
     if ($allele =~ /\(|\)/g) {
       $allele =~ /\((\S+)\)(\d+)/;
       my ($al,$nu);
       if ($1 and $2) {
         $al = $1;
         $nu = $2;
		 next if $al =~ /\(/;
         next if ($al and $nu and $allele_string =~ /$al/);
       }
     }
     if ($allele_string !~ /$allele/i) {
        $frequency = '\N' if ! defined $frequency;
        $sample_id = '\N' if ! defined $sample_id;
        print OUT "$allele_id\t$variation_id\t$subsnp_id\t$allele\t$frequency\t$sample_id\t$allele_string\n";
     }
     #$matched_var{$variation_id}++;
  }

  #create table allele_change to hold them
  system("mv $TMP_DIR/$vdbname\_allele_$check $TMP_DIR/$TMP_FILE");

  create_and_load($dbVar,"$vdbname\_allele_$check","allele_id i*","variation_id i*","subsnp_id i","allele","frequency","sample_id i","allele_string");
}

sub check_genotype_table {

  my $dbVar = shift;
  my $check = shift;
  
  debug("Checking for genotype tables...All genotypes in genotype table should match alleles in allele_string of variation_feature table for the variation, also check if we have four genotypes for the same variation");

  my $hap_line;
  if ($hap_id_string =~ /\(\d+\)/) {
    $hap_line = "AND vf.seq_region_id not in $hap_id_string";
  }
  else {
    $hap_line ='';
  }

  #foreach my $table ("yuan_hum_var_51.tmp_individual_genotype_single_bp","individual_genotype_multiple_bp","population_genotype") {  
  foreach my $table ("population_genotype") {  
    debug("checking table $table...");

    open OUT, ">$TMP_DIR/$vdbname\_$table\_$check";
    #only find out variation with one mapping that is vf.map_weight=1 and seq_region_id not in(226063,226064,226065,226066
    my $sth = $dbVar->prepare(qq{select gt.*, vf.allele_string, vf.seq_region_id from $table gt, variation_feature vf where gt.variation_id = vf.variation_id and vf.map_weight=1 $hap_line}, {mysql_use_result=>1});

    my (%rec_var,$population_genotype_id,$variation_id,$subsnp_id,$allele_1,$allele_2,$allele_string,$seq_region_id,$frequency,$sample_id);
    $sth->execute();
    if ($table =~ /pop/i) {
      $sth->bind_columns(\$population_genotype_id,\$variation_id,\$subsnp_id,\$allele_1,\$allele_2,\$frequency,\$sample_id,\$allele_string,\$seq_region_id);
    }
    else {
      $sth->bind_columns(\$variation_id,\$subsnp_id,\$allele_1,\$allele_2,\$sample_id,\$allele_string,\$seq_region_id);
    }
    LINE : while($sth->fetch()) {
      next if $allele_1 =~ /N/i;
      next if $allele_2 =~ /N/i;
      next if $allele_1 =~ /\+/;
      next if $allele_2 =~ /\+/;
      next if $allele_1 =~ /^\(.*\)$/;
      next if $allele_2 =~ /^\(.*\)$/;
      #next if (!$allele_1 or $allele_1 =~ /\+|N|^\(.*\)$/ig);
      #next if (!$allele_2 or $allele_2 =~ /\+|N|^\(.*\)$/ig);
      my ($pass_allele);
      foreach my $allele ($allele_1,$allele_2) {
        my ($al,$nu);
        if ($allele =~ /\(|\)/g) {
          $allele =~ /\((\S+)\)(\d+)/;
          if ($1 and $2) {
            $al = $1;
            $nu = $2;
          }
        }
        $pass_allele++ if ($al and $nu and $allele_string =~ /$al/);
        next LINE if ($al and $nu and $pass_allele == 2);
      }
      if ($allele_string and ($allele_string !~ /$allele_1/i or $allele_string !~ /$allele_2/i)) {
        if ($table =~ /pop/i) {
          print OUT "$population_genotype_id\t$variation_id\t$subsnp_id\t$allele_1\t$allele_2\t$frequency\t$sample_id\t$allele_string\t$seq_region_id\n";
        }
        else {
          print OUT "$variation_id\t$subsnp_id\t$allele_1\t$allele_2\t$sample_id\t$allele_string\t$seq_region_id\n";
        }
      }
    }
    
    #system("mv $TMP_DIR/$vdbname\_$table\_$check $TMP_DIR/$TMP_FILE");
    
    if ($table =~ /pop/i) {
      create_and_load($dbVar,"$vdbname\_$table\_$check","population_genotype_id i","variation_id i*","subsnp_id i","allele_1","allele_2","frequency","sample_id i","allele_string","seq_region_id");
    }
    else {
      create_and_load($dbVar,"$vdbname\_$table\_$check","variation_id i*","subsnp_id i","allele_1","allele_2","sample_id i","allele_string","seq_region_id");
    }
  }
}

sub check_variation_feature_table {

  my $dbVar = shift;
  open OUT, ">$TMP_DIR/variation_feature_check\_$vdbname";
  debug("Checking for variation_feature table...All alleles in allele_string of variation_feature table should less than four for the variation");  
    
  my $sth = $dbVar->prepare(qq{select variation_id, allele_string from variation_feature where length(allele_string) - length(REPLACE(allele_string,'/','')) > 2}, {mysql_use_result=>1});  
    
  my ($variation_id,$allele_string);  
  $sth->execute();  
  $sth->bind_columns(\$variation_id,\$allele_string);  
  while($sth->fetch()) {
    print OUT "variation_id\t$allele_string\n";
  }
}

sub change_allele_table {

  my $dbVar = shift;
  open OUT, ">$TMP_DIR/allele_change\_$vdbname";
    
  debug("Change allele table...");
  my $allele_table = "allele";
  my $sth = $dbVar->prepare(qq{select allele_id,allele from $allele_table where allele like "%)%("}, {mysql_use_result=>1});

  my ($allele_id,$allele);
  $sth->execute();
  $sth->bind_columns(\$allele_id,\$allele);
  while($sth->fetch()) {
    ###3)GCTGG(
    if ($allele =~ /(\d+)\)(\S+)\(/ ) {
      $allele = "(".reverse($2).")".reverse($1);
      print OUT "update allele set allele = \"$allele\" where allele_id=$allele_id;\n";
    }
    elsif($allele =~ /^(\).*\()$/) {
      $allele = reverse($1);
      print OUT "update allele set allele = \"$allele\" where allele_id=$allele_id;\n";
    }
  }
}

sub change_allele_table1 {  

  my $dbVar = shift;  
  open OUT, ">$TMP_DIR/allele_change1";  
  debug("Change allele table1...");  
  my $allele_table = "allele";  
  my $sth = $dbVar->prepare(qq{select allele_id,allele,vf.allele_string from $allele_table a,variation_feature vf, varid_one_allele v where allele like "(%)%" and a.variation_id=v.variation_id and v.variation_id=vf.variation_id}, {mysql_use_result=>1});

  my ($allele_id,$allele,$allele_string);
  $sth->execute();
  $sth->bind_columns(\$allele_id,\$allele,\$allele_string);
  while($sth->fetch()) {
  ###3)GCTGG(
    #my @alleles = split /\//, $allele_string;
    if ($allele_string !~ /$allele/) {
      $allele =~ /(\(\S+\))(\d+)/;
      $allele = "$1".reverse($2);
      print OUT "update allele set allele = \"$allele\" where allele_id=$allele_id;\n";
    }
  }
}

sub change_allele_table2 {  

  my $dbVar = shift;  
  open OUT, ">$TMP_DIR/allele_change3";  
  debug("Change allele table2...");  
  my $allele_table = "allele";  
  #my $sth = $dbVar->prepare(qq{select allele_id,variation_id,allele,allele_string from $allele_table}, {mysql_use_result=>1});
  my $sth = $dbVar->prepare(qq{select a.allele_id,v.variation_id,a.allele,vf.allele_string from $allele_table a,variation_feature vf, varid_one_allele v where allele like "(%)%" and a.variation_id=v.variation_id and v.variation_id=vf.variation_id}, {mysql_use_result=>1});
  
  my ($allele_id,$allele,$allele_string,$variation_id,%done);
  $sth->execute();
  $sth->bind_columns(\$allele_id,\$variation_id,\$allele,\$allele_string);
  while($sth->fetch()) {
  ###3)GCTGG(
    print "$allele_id,$allele,$allele_string\n";
    my @alleles = split /\//, $allele_string;
    if ($alleles[0] eq $allele) {
      if (! $done{$variation_id}{$allele}) {
        $done{$variation_id}{$allele}=1;
      }
      else {
        print OUT "update $allele_table set allele = \"$alleles[1]\" where allele_id=$allele_id;\n";
      }
    }
    elsif ($alleles[1] eq $allele) { 
      if (! $done{$variation_id}{$allele}) {        
        $done{$variation_id}{$allele}=1;      
      }      
      else {        
        print OUT "update $allele_table set allele = \"$alleles[0]\" where allele_id=$allele_id;\n";
      }
    }
  }
}

sub change_allele_table3 {  
  my $dbVar = shift;  
  open OUT, ">$TMP_DIR/allele_change3";  
  debug("Change allele table3...");  
  my $allele_table = "allele";  
  #my $sth = $dbVar->prepare(qq{select allele_id,variation_id,allele,allele_string fro+m $allele_table}, {mysql_use_result=>1});  
  my $sth = $dbVar->prepare(qq{select a.allele_id,v.variation_id,a.allele,vf.allele_string from $allele_table a,variation_feature vf, varid_one_allele v where allele like "(%)%" and a.variation_id=v.variation_id and v.variation_id=vf.variation_id}, {mysql_use_result=>1});

  my ($allele_id,$allele,$allele_string,$variation_id,%done);
  $sth->execute();
  $sth->bind_columns(\$allele_id,\$variation_id,\$allele,\$allele_string);
  while($sth->fetch()) {
  ###3)GCTGG(
    print "$allele_id,$allele,$allele_string\n";
    my @alleles = split /\//, $allele_string;
    if ($alleles[0] eq $allele) {
      print OUT "insert into allele (variation_id,allele) values($variation_id,\"$alleles[1]\");\n";
    }
    else {
      print OUT "insert into allele (variation_id,allele) values($variation_id,\"$alleles[0]\");\n";                 
    }
  }
}

sub delete_allele_gtype_tables {
  my $dbVar = shift;  
  my $buffer = {}; #hash containing all the files to be parallelized
  my ($table_name,$population_genotype_id,$variation_id,$ss_id,$allele_id,$allele,$allele_1,$allele_2,$frequency,$sample_id,$allele_string,$seq_region_id);

  my (%rec_seq_region_ids,%rec_table_name);
  if ($species =~ /hum|homo/i) {
    my $sth1 = $dbCore->dbc->prepare(qq{SELECT sr.seq_region_id, sr.name
                                       FROM seq_region sr, coord_system cs 
                                       WHERE sr.coord_system_id=cs.coord_system_id 
                                       AND cs.name='chromosome'
                                       AND cs.version='GRCh37'});
    $sth1->execute();
  
    while (my ($seq_reg_id,$name) = $sth1->fetchrow_array()) {
      $rec_seq_region_ids{$seq_reg_id}=$name;
    }
  }
  
  foreach my $table ("tmp_individual_genotype_single_bp","individual_genotype_multiple_bp","population_genotype","allele") {  #we only change things that have map_weight=1
  #foreach my $table ("population_genotype","tmp_individual_genotype_single_bp") { 
  #foreach my $table ("individual_genotype_multiple_bp") {
    my $allele_column = "allele_1";
    debug("processling table $table");    

    if ($table =~ /allele/) {
      $allele_column = "allele";
    }
      
    $dbVar->do(qq{DELETE FROM $table\_check WHERE $allele_column like "(%)"});
    
    my $sth = $dbVar->prepare(qq{select v.* from $table\_check v}, {mysql_use_result=>1});

    $sth->execute();
  
    if ($table =~ /pop/i) {
      $sth->bind_columns(\$population_genotype_id,\$variation_id,\$ss_id,\$allele_1,\$allele_2,\$frequency,\$sample_id,\$allele_string,\$seq_region_id);
    }
    elsif($table =~ /allele/i) {
      $sth->bind_columns(\$allele_id,\$variation_id,\$ss_id,\$allele,\$frequency,\$sample_id,\$allele_string);
    }
    else {
      $sth->bind_columns(\$variation_id,\$ss_id,\$allele_1,\$allele_2,\$sample_id,\$allele_string,\$seq_region_id);
    }

    LINE : while($sth->fetch()) {
      if ($table =~ /allele/i and $allele_id) {
	print_buffered($buffer, "$TMP_DIR/$table\_del",
		       "DELETE FROM $table WHERE allele_id = $allele_id;\n");
      }
      elsif ($table =~ /pop/ and $population_genotype_id) {
	print_buffered($buffer, "$TMP_DIR/$table\_del",
		       "DELETE FROM $table WHERE population_genotype_id = $population_genotype_id;\n");
      }
      elsif ($variation_id and $allele_1 and $allele_2 and $sample_id) { 
	print_buffered($buffer, "$TMP_DIR/$table\_del",
		       "DELETE FROM $table WHERE variation_id = $variation_id and allele_1 = \"$allele_1\" and allele_2 = \"$allele_2\" and sample_id = $sample_id;\n");
      }
    }
    print_buffered($buffer);
  }
}



sub change_allele_gtype_tables {
  my $dbVar = shift;  
  my $buffer = {}; #hash containing all the files to be parallelized
  my ($table_name,$population_genotype_id,$variation_id,$ss_id,$allele_id,$allele,$allele_1,$allele_2,$frequency,$sample_id,$allele_string,$seq_region_id);

  my (%rec_seq_region_ids,%rec_table_name);
  if ($species =~ /hum|homo/i) {
    my $sth1 = $dbCore->dbc->prepare(qq{SELECT sr.seq_region_id, sr.name
                                       FROM seq_region sr, coord_system cs 
                                       WHERE sr.coord_system_id=cs.coord_system_id 
                                       AND cs.name='chromosome'
                                       AND cs.version='GRCh37'});
    $sth1->execute();
  
    while (my ($seq_reg_id,$name) = $sth1->fetchrow_array()) {
      $rec_seq_region_ids{$seq_reg_id}=$name;
    }
  }
  
  foreach my $table ("tmp_individual_genotype_single_bp") {
  #foreach my $table ("tmp_individual_genotype_single_bp","individual_genotype_multiple_bp","population_genotype","allele") {#we only change things that have map_weight=1
  #foreach my $table ("ind_gtype_minus_strand_A_T_OR_C_G","pop_gtype_minus_strand_A_T_OR_C_G") { 
  foreach my $table ("population_genotype") {
    my $allele_column = "allele_1";
    debug("processling table $table");    

    open OUT, ">$TMP_DIR/$table\_change\_$vdbname" or die "can't open file $table\_changed : $!";
    open OUT1, ">$TMP_DIR/$table\_change\_$vdbname\_del" or die "can't open file $table\_del : $!";

    if ($table =~ /allele/) {
      $allele_column = "allele";
    }
    #delete allele = (INDETERMINATE)  
    $dbVar->do(qq{DELETE FROM $vdbname\_$table\_$check WHERE $allele_column like "(%)"});
    print "processing table $vdbname\_$table\_$check\n";
    my $sth = $dbVar->prepare(qq{select v.* from $vdbname\_$table\_$check v where v.allele_2 != "R" and v.allele_2 != "Y"}, {mysql_use_result=>1});

    $sth->execute();
  
    if ($table =~ /pop/i) {
      $sth->bind_columns(\$population_genotype_id,\$variation_id,\$ss_id,\$allele_1,\$allele_2,\$frequency,\$sample_id,\$allele_string,\$seq_region_id);
    }
    elsif($table =~ /allele/i) {
      $sth->bind_columns(\$allele_id,\$variation_id,\$ss_id,\$allele,\$frequency,\$sample_id,\$allele_string);
    }
    else {
      $sth->bind_columns(\$variation_id,\$ss_id,\$allele_1,\$allele_2,\$sample_id,\$allele_string,\$seq_region_id);
    }

    LINE : while($sth->fetch()) {
      my ($old_allele,$old_allele_1,$old_allele_2);
      next if ($allele_1 and $allele_1 =~ /\+/);
      next if ($allele_2 and $allele_2 =~ /\+/);
      next if ($allele_1 and $allele_1 =~ /N|^\(.*\)$/ig);
      next if ($allele_2 and $allele_2 =~ /N|^\(.*\)$/ig);
      next if ($allele_1 and $allele_1 =~/homozygous|indeterminate/ig or $allele_2 and $allele_2 =~/homozygous|indeterminate/ig);
      next if ($allele_1 and $allele_1 =~/SEQUENCE/ig or $allele_2 and $allele_2 =~/SEQUENCE/ig);
      next if ($allele and $allele =~ /\+/);
      next if ($allele and $allele =~ /N|^\(.*\)$/ig);
      my @alleles = split /\//, $allele_string;
			
      foreach my $a (@alleles) {
	if ($a =~ /\(|\)/g) {
          $a =~ /\((\S+)\)(\d+)/;
          if ($1 and $2) {
	    $a = $1;
	    reverse_comp(\$a);
	    $a = "($a)$2";
	  }
        }
	else {
	  reverse_comp(\$a);
	}
      }
      my $new_allele_string = join "/", @alleles;
      
      if ($table =~ /allele/) {
	$old_allele = $allele;
        $frequency = '\N' if ! defined $frequency;
        $sample_id = '\N' if ! defined $sample_id;
        #cases like (AC)3 as allele, and allele_string like (AC)3/(AC)4
        my ($al,$nu,$allele_new);
        if ($allele =~ /\(|\)/g) {
          $allele =~ /\((\S+)\)(\d+)/;
          if ($1 and $2) {
            $al = $1;
            $nu = $2;
          }
        }
        next if ($al and $nu and $allele_string =~ /$al/);
        if ($allele and $allele_string !~ /$allele/) {
	  if ($allele =~ /\(|\)/g) {
            $allele =~ /\((\S+)\)(\d+)/;
            if ($1 and $2) {
	      $allele = $1;
	      reverse_comp(\$allele);
	      $allele_new = "($allele)$2";
	    }
          }
	  else {
	    reverse_comp(\$allele);
            $allele_new = $allele;
	  } 
          if ($allele_string =~ /$allele/) {
	    print OUT "update $table set allele = \"$allele_new\" where allele_id=$allele_id;\n";
          }
	  else {
	    print OUT1 "DELETE FROM $table WHERE allele_id = $allele_id;\n";
          }
	}
      }
      else {#for genotype tables

	$old_allele_1 = $allele_1;
	$old_allele_2 = $allele_2;
	my ($new_allele_1,$new_allele_2);
        #cases like (AC)3 as allele_1, and allele_string like (AC)3/(AC)4
        my ($allele_pass);
	foreach my $allele ($allele_1,$allele_2) {
	  my ($al,$nu);
          if ($allele =~ /\(|\)/g) {
	    $allele =~ /\((\S+)\)(\d+)/;
	    if ($1 and $2) {
	      $al = $1;
	      $nu = $2;
	    }
	  } 
          #if $allele_string matching $al, it will keep unchanged, otherwise will go next if
	  $allele_pass++ if ($al and $nu and $allele_string =~ /$al/);
          #if both allele_1 and allele_2 with allele part matched by allele_string, then next
          next LINE if $allele_pass and $allele_pass == 2;
	}
	#print "variation_id is $variation_id and allele_1 is $allele_1 allele_2 is $allele_2 and allele_string is $allele_string\n";
        if ($allele_string !~ /$allele_1/ or $allele_string !~ /$allele_2/) {
	
	  if ($allele_1 =~ /\(|\)/g) {
            $allele_1 =~ /\((\S+)\)(\d+)/;
            if ($1 and $2) {
	      $allele_1 = $1;
	      reverse_comp(\$allele_1);
	      $new_allele_1 = "($allele_1)$2";
	    }
          }
	  else {
	    reverse_comp(\$allele_1);
	    $new_allele_1 = $allele_1;
	  }
	  if ($allele_2 =~ /\(|\)/g) {
            $allele_2 =~ /\((\S+)\)(\d+)/;
	    if ($1 and $2) {
              $allele_2 = $1;
	      reverse_comp(\$allele_2);
	      $new_allele_2 = "($allele_2)$2";
	    }
          }
	  else {
	    reverse_comp(\$allele_2);
	    $new_allele_2 = $allele_2;
	  }
          if ($allele_string =~ /$allele_1/ and $allele_string =~ /$allele_2/) {
	  #print "allele_1 is $allele_1 and allele_string is $allele_string\n";
	    if ($table =~ /pop/i) {
	      print OUT "update $table set allele_1 = \"$new_allele_1\", allele_2 = \"$new_allele_2\" where population_genotype_id = $population_genotype_id;\n";
	    }
	    elsif ($table =~ /multi/) {#for individual_genotype_multiple_base table
	      print OUT "update $table set allele_1 = \"$new_allele_1\", allele_2 = \"$new_allele_2\" where variation_id=$variation_id and subsnp_id = $ss_id and sample_id=$sample_id and allele_1 = \"$old_allele_1\" and allele_2 = \"$old_allele_2\";\n";
	    }
	    elsif ($table =~ /single/) {
              if ($species =~ /hum|homo/i) {
                my $chr_name = $rec_seq_region_ids{$seq_region_id};
                if ($chr_name) {
                  $table_name = "$table\_SubInd_ch$chr_name";
                }
              }
              else {
                $table_name = $table;
              }
              #delete from indiviaul_genotype table rather than sub table to get rid of multi mappings to non ref assembly
              print_buffered($buffer, "$TMP_DIR/$table_name\_del",
			     "DELETE FROM $table WHERE variation_id=$variation_id and sample_id=$sample_id;\n");
              print_buffered( $buffer, "$TMP_DIR/$table_name\.import",
              join ("\t",$variation_id,$ss_id,$new_allele_1,$new_allele_2,$sample_id)."\n");
              #print_buffered( $buffer, "$TMP_DIR/$table_name\_change","update $table_name set allele_1 = \"$new_allele_1\", allele_2 = \"$new_allele_2\" where variation_id=$variation_id and sample_id=$sample_id and allele_1 = \"$old_allele_1\" and allele_2 = \"$old_allele_2\";\n");
            }
	  }
          else {#if after reverse comp allele_1 and allele_2, allele_string still not match them, then :
            if ($table =~ /pop/) {
              print OUT1 "DELETE FROM $table WHERE population_genotype_id=$population_genotype_id;\n";
	    }
	    else {
	      if ($old_allele_1 !~ /y|r/i and $old_allele_2 !~ /y|r/i) {
		print_buffered($buffer,"$TMP_DIR/$table\_remove",
			       "DELETE FROM $table WHERE variation_id=$variation_id and allele_1 = \"$old_allele_1\" and allele_2=\"$old_allele_2\" and sample_id=$sample_id;\n");
	      }
	    }
	  }
	}
      }
    }#close while loop
    if ($table =~ /single/) {
      &print_buffered($buffer);
    }
  }
}

sub print_buffered {
    my $buffer = shift;
    my $filename = shift;
    my $text = shift;

    local *FH;

    if( ! $filename ) {
      # flush the buffer
      foreach my $file (keys %{$buffer}){
        open( FH, ">>$file" ) or die "Could not print to file $! \n";
        print FH $buffer->{ $file };
        close FH;
      }
      %{$buffer} = ();
    } 
    else {
      $buffer->{ $filename } .= $text;
      if( length( $buffer->{ $filename } ) > 10_000 ) {
        open( FH, ">>$filename" ) or die;
        print FH $buffer->{ $filename };
        close FH;
        $buffer->{ $filename } = '';
     }
   }
}

sub change_variation_feature_table {

  my $dbVar = shift;
  open OUT, ">$TMP_DIR/variation_feature_change\_$vdbname";
  debug("Change variation_feature table...");
  my $variation_feature_table = "variation_feature";
  my $sth = $dbVar->prepare(qq{select variation_feature_id,allele_string from $variation_feature_table where allele_string like "%)%("}, {mysql_use_result=>1});

  my ($allele_string,$variation_feature_id,);
  $sth->execute();
  $sth->bind_columns(\$variation_feature_id,\$allele_string);
  while($sth->fetch()) {
  ###3)GCTGG(
    my @alleles = split /\//, $allele_string;
    foreach my $allele (@alleles) {
      if ($allele =~ /^\).*\($/ ) {
        reverse_comp($allele);
      }
    }
    my $new_allele_string = join "/", @alleles;
    print OUT "update variation_feature set allele_string = \"$new_allele_string\" where variation_feature_id=$variation_feature_id;\n";
  }
}
