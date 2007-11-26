#! /usr/local/ensembl/bin/perl

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

our ($TMP_DIR, $TMP_FILE, $SEQ_REGION_ID, $SEQ_REGION_NAME, $species, $VAR_DBNAME);



GetOptions('species=s'   => \$species,
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

#check_allele_table($dbVar);
#check_genotype_table($dbVar);
#check_variation_feature_table($dbVar);
#change_allele_table($dbVar);
#change_allele_table1($dbVar);
#change_allele_table2($dbVar);
#change_allele_table3($dbVar);
#change_allele_gtype_tables($dbVar);
#change_variation_feature_table($dbVar);
#rename_tmp_gtype_tables($dbVar);
drop_tables($dbVar);

sub check_allele_table {

  my $dbVar = shift;

  debug("Checking for allele table...One variation should have at least two alleles");
  open OUT, ">$TMP_DIR/$vdbname\_allele_check";
  open OUT1,">$TMP_DIR/$vdbname\_allele_check_one_allele";

  #do not use this command, it's too big
  #$dbVar->do(qq{create table varid_one_allele select variation_id, count(*) as count from (select variation_id,allele from allele group by variation_id,allele) as tmp_table group by variation_id having count =1});
  #my $varid_ref = $dbVar->db_handle->selectall_arrayref(qq{select count(*) from varid_one_allele});
  #my $count = $varid_ref->[0][0];
  #print "There is $count record in table varid_one_allele\n";

  debug("Checking for allele table...All alleles in allele_string of variation_feature table should match alleles in allele table for the variation"); 

  #only find out variation with one mapping that is vf.map_weight=1 and seq_region_id not in (226063,226064,226065,226066)
  my $sth = $dbVar->prepare(qq{select a.*,vf.allele_string from allele a, variation_feature vf where a.variation_id = vf.variation_id and vf.map_weight=1 and seq_region_id not in (226063,226064,226065,226066)}, {mysql_use_result=>1});
  
  my (%total_var,%matched_var,$variation_id,$allele,$allele_string,$allele_id,$frequency,$sample_id);
  $sth->execute();
  $sth->bind_columns(\$allele_id,\$variation_id,\$allele,\$frequency,\$sample_id,\$allele_string);
  while($sth->fetch()) {
     #print "$variation_id,$allele,$allele_string\n";
     next if ($allele =~ /N|^\).*\($|^\(.*\)|\+/ig);
     #print "$variation_id,$allele,$allele_string again\n";
     $total_var{$variation_id}=1;
     if ($allele_string !~ /$allele/i) {
        $frequency = '\N' if ! defined $frequency;
        $sample_id = '\N' if ! defined $sample_id;
        print OUT "$allele_id\t$variation_id\t$allele\t$frequency\t$sample_id\t$allele_string\n";
     }
     $matched_var{$variation_id}++;
  }
  foreach my $variation_id (keys %total_var) {
    if ($matched_var{$variation_id} and $matched_var{$variation_id} <2) {
      print OUT1 "$variation_id\n";
    }
  }
  #create table allele_change to hold them
  system("mv $TMP_DIR/$vdbname\_allele_check $TMP_DIR/$TMP_FILE");
  create_and_load($dbVar,"$vdbname\_allele_check","allele_id i*","variation_id i*","allele","frequency","sample_id i","allele_string");
  #system("mv $TMP_DIR/$vdbname\_allele_check_one_allele $TMP_DIR/$TMP_FILE");
  #create_and_load($dbVar,"$vdbname\_allele_check_one_allele","variation_id i*");
  #most variations do have two alleles, such as - and (ALU) which has been excluded before counting
  #$dbVar->do(qq{create table varid_one_allele_new select variation_id, count(*) as count from (select a.variation_id,a.allele from allele a, allele_check_one_allele ac where a.variation_id=ac.variation_id group by a.variation_id,a.allele) as tmp_table group by variation_id having count =1});
  #I changed those 8 variations by hand this time:
  #[ens-genomics1]yuan_hum_var_46a>select * from allele_check_one_allele_new ;
  #[ens-genomics1]yuan_hum_var_46a>select * from variation_feature where variation_id=11811030;
  #[ens-genomics1]yuan_hum_var_46a>select * from tmp_var_allele where refsnp_id=41563015;
  #[cbi3]dbSNP_127_human_9606>select * from AlleleFreqBySsPop where subsnp_id=68074678;
  #[ens-genomics1]yuan_hum_var_46a>select * from allele where variation_id=11811030;
  #[ens-genomics1]yuan_hum_var_46a>insert into allele (variation_id,allele,frequency,sample_id)values(11811030,"CGGGGCCTGCGCCTGCGCGCTCAGCGGCCGG",0.9130435,714);
}

sub check_genotype_table {

  my $dbVar = shift;
  #system("mv $TMP_DIR/$vdbname\_tmp_individual_genotype_single_bp_check $TMP_DIR/$TMP_FILE");
  #create_and_load($dbVar,"tmp_individual_genotype_single_bp_check","variation_id i*","allele_1","allele_2","sample_id i","allele_string","seq_region_id");
  #exit;
  debug("Checking for genotype tables...All genotypes in genotype table should match alleles in allele_string of variation_feature table for the variation, also check if we have four genotypes for the same variation");

  foreach my $table ("tmp_individual_genotype_single_bp_more_insert") {  
  #foreach my $table ("tmp_individual_genotype_single_bp","individual_genotype_multiple_bp","population_genotype") {
    debug("checking table $table...");

    open OUT, ">$TMP_DIR/$vdbname\_$table\_check";
    #only find out variation with one mapping that is vf.map_weight=1 and seq_region_id not in(226063,226064,226065,226066
    my $sth = $dbVar->prepare(qq{select gt.*, vf.allele_string, vf.seq_region_id from $table gt, $vdbname.variation_feature vf where gt.variation_id = vf.variation_id and vf.map_weight=1 and vf.seq_region_id not in (226063,226064,226065,226066)}, {mysql_use_result=>1});

    my (%rec_var,$population_genotype_id,$variation_id,$allele_1,$allele_2,$allele_string,$seq_region_id,$frequency,$sample_id);
    $sth->execute();
    if ($table =~ /pop/i) {
      $sth->bind_columns(\$population_genotype_id,\$variation_id,\$allele_1,\$allele_2,\$frequency,\$sample_id,\$allele_string,\$seq_region_id);
    }
    else {
      $sth->bind_columns(\$variation_id,\$allele_1,\$allele_2,\$sample_id,\$allele_string,\$seq_region_id);
    }
    while($sth->fetch()) {
      next if $allele_1 =~ /N/i;
      next if $allele_2 =~ /N/i;
      next if $allele_1 =~ /\+/;
      next if $allele_2 =~ /\+/;
      next if $allele_1 =~ /^\(.*\)/;
      next if $allele_2 =~ /^\(.*\)/;
      next if (!$allele_1 or $allele_1 =~ /\+|N|^\(.*\)$/ig);
      next if (!$allele_2 or $allele_2 =~ /\+|N|^\(.*\)$/ig);
      if ($allele_string and ($allele_string !~ /$allele_1/i or $allele_string !~ /$allele_2/i)) {
        if ($table =~ /pop/i) {
          print OUT "$population_genotype_id\t$variation_id\t$allele_1\t$allele_2\t$frequency\t$sample_id\t$allele_string\t$seq_region_id\n";
        }
        else {
          print OUT "$variation_id\t$allele_1\t$allele_2\t$sample_id\t$allele_string\t$seq_region_id\n";
        }
      }
      #$rec_var{$variation_id}{$allele_1}++;
      #$rec_var{$variation_id}{$allele_2}++;
    }
    
    #foreach my $variation_id (keys %rec_var) {
    #  my @alleles = keys %{$rec_var{$variation_id}};
    #  if (scalar @alleles ==4) {
    #    print "$variation_id has 4 genotypes in table $table\n";
    #  }
    #}

    system("mv $TMP_DIR/$vdbname\_$table\_check $TMP_DIR/$TMP_FILE");
    
    if ($table =~ /pop/i) {
      create_and_load($dbVar,"$vdbname\_$table\_check","population_genotype_id i","variation_id i*","allele_1","allele_2","frequency","sample_id i","allele_string","seq_region_id");
    }
    else {
      create_and_load($dbVar,"$vdbname\_$table\_check2","variation_id i*","allele_1","allele_2","sample_id i","allele_string","seq_region_id");
    }
  }
}

sub check_variation_feature_table {

  my $dbVar = shift;
  open OUT, ">$TMP_DIR/variation_feature_check\_$vdbname";
  debug("Checking for variation_feature table...All alleles in allele_string of variation_feature table should less than four for the variation");  
    
  my $sth = $dbVar->prepare(qq{select variation_id, allele_string from variation_feature where allele_string like "_/_/_/_"}, {mysql_use_result=>1});  
    
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
      $allele = "($2)".reverse($1);
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

sub change_allele_gtype_tables {
  my $dbVar = shift;  
  my $buffer = {}; #hash containing all the files to be parallelized
  my ($table_name,$population_genotype_id,$variation_id,$allele_id,$allele,$allele_1,$allele_2,$frequency,$sample_id,$allele_string,$seq_region_id);

  my (%rec_seq_region_ids,%rec_table_name);
  if ($species =~ /hum|homo/i) {
  my $sth1 = $dbCore->dbc->prepare(qq{SELECT sr.seq_region_id, sr.name
                                       FROM seq_region sr, coord_system cs 
                                       WHERE sr.coord_system_id=cs.coord_system_id 
                                       AND cs.name='chromosome'
                                       AND cs.version='NCBI36'});
  $sth1->execute();
  
  while (my ($seq_reg_id,$name) = $sth1->fetchrow_array()) {
    $rec_seq_region_ids{$seq_reg_id}=$name;
  }
  }
  
  #foreach my $table ("tmp_individual_genotype_single_bp","individual_genotype_multiple_bp","population_genotype","allele") {#we only change things that have map_weight=1
  foreach my $table ("population_genotype") { 
  #foreach my $table ("individual_genotype_multiple_bp","population_genotype","allele") {
    my $allele_column = "allele_1";
    debug("processling table $table");    
    #if ($table !~ /single/) {
      open OUT, ">$TMP_DIR/$table\_change\_$vdbname" or die "can't open file $table\_out : $!";
    #}
    if ($table =~ /allele/) {
      $allele_column = "allele";
    }
      
    $dbVar->do(qq{DELETE FROM $table\_check WHERE $allele_column like "(%)"});
    
    my $sth = $dbVar->prepare(qq{select v.* from $table\_check v}, {mysql_use_result=>1});

    $sth->execute();
  
    if ($table =~ /pop/i) {
      $sth->bind_columns(\$population_genotype_id,\$variation_id,\$allele_1,\$allele_2,\$frequency,\$sample_id,\$allele_string,\$seq_region_id);
    }
    elsif($table =~ /allele/i) {
      $sth->bind_columns(\$allele_id,\$variation_id,\$allele,\$frequency,\$sample_id,\$allele_string);
    }
    else {
      $sth->bind_columns(\$variation_id,\$allele_1,\$allele_2,\$sample_id,\$allele_string,\$seq_region_id);
    }

    while($sth->fetch()) {
      my ($old_allele,$old_allele_1,$old_allele_2);
      next if ($allele_1 and $allele_1 =~ /\+|N|^\(.*\)$/ig);
      next if ($allele_2 and $allele_2 =~ /\+|N|^\(.*\)$/ig);
      next if ($allele and $allele =~ /\+|N|^\(.*\)$/ig);
      my @alleles = split /\//, $allele_string;
      foreach my $a (@alleles) {
	if ($a =~ /\((\S+)\)(\d+)/) {
	  $a = $1;
	  reverse_comp(\$a);
	  $a = "($a)$2";
	}
	else {
	  reverse_comp(\$a);
	}
      }
      my $new_allele_string = join "/", @alleles;
      
      if ($table =~ /allele/) {
	$old_allele = $allele;
	if ($allele_string !~ /$allele/ or ($new_allele_string =~ /$allele/ and $allele_string =~ /$allele/)) {
	  if ($allele =~ /\((\S+)\)(\d+)/) {
	    $allele = $1;
	    reverse_comp(\$allele);
	    $allele = "($allele)$2";
	  }
	  else {
	    reverse_comp(\$allele);
	  }
          if ($allele_string =~ /$allele/) {
	    print OUT "update $table set allele = \"$allele\" where allele_id=$allele_id;\n";
          }
	}
      }
      else {
	next if ($allele_1 =~ /\-/ and $allele_2 =~ /\-/);
	next if ($allele_1 =~/\+|^\(.*\)$/g or $allele_2 =~/\+|^\(.*\)$/g);
	next if ($allele_1 =~/homozygous|indeterminate|SEQUENCE/ig or $allele_2 =~/homozygous|indeterminate|SEQUENCE/ig);
	$old_allele_1 = $allele_1;
	$old_allele_2 = $allele_2;
	if (($allele_string !~ /$allele_1|$allele_2/) or ($new_allele_string =~ /$allele_1|$allele_2/ and $allele_string =~ /$allele_1|$allele_2/)) {
	  if ($allele_1 =~ /\((\S+)\)(\d+)/) {
	    $allele_1 = $1;
	    reverse_comp(\$allele_1);
	    $allele_1 = "($allele_1)$2";
	  }
	  else {
	    reverse_comp(\$allele_1);
	  }
	  if ($allele_2 =~ /\((\S+)\)(\d+)/) {
	    $allele_2 = $1;
	    reverse_comp(\$allele_2);
	    $allele_2 = "($allele_1)$2";
	  }   
	  else {
	    reverse_comp(\$allele_2);
	  }  
          if ($allele_string =~ /$allele_1/ and $allele_string =~ /$allele_2/) {
	    if ($table =~ /pop/i) {
	      print OUT "update $table set allele_1 = \"$allele_1\", allele_2 = \"$allele_2\" where population_genotype_id = $population_genotype_id;\n";
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
              print_buffered( $buffer, "$TMP_DIR/$table_name\_change","update $table_name set allele_1 = \"$allele_1\", allele_2 = \"$allele_2\" where variation_id=$variation_id and sample_id=$sample_id and allele_1 = \"$old_allele_1\" and allele_2 = \"$old_allele_2\";\n");
            }
	    else {#for individual_genotype_multiple_base table
	      print OUT "update $table set allele_1 = \"$allele_1\", allele_2 = \"$allele_2\" where variation_id=$variation_id and sample_id=$sample_id and allele_1 = \"$old_allele_1\" and allele_2 = \"$old_allele_2\";\n";
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

sub rename_tmp_gtype_tables {

  my $dbVar = shift;
  foreach my $table ("SubInd_ch") { 
    my $single_table_ref = $dbVar->db_handle->selectall_arrayref(qq{SHOW tables like "$table%"});
    my @tables = map {$_->[0] } @$single_table_ref;
    foreach my $table (@tables) {
      debug("rename table $table...");
      $dbVar->do(qq{RENAME TABLE $table TO yuan_hum_gtype_48.$table});
    }
  }
}

sub drop_tables {

  my $dbVar = shift;
  foreach my $table ("varid","old","tmp_ids","47") { 
    my $single_table_ref = $dbVar->db_handle->selectall_arrayref(qq{SHOW tables like "%$table%"});
    my @tables = map {$_->[0] } @$single_table_ref;

    my $table_names = join ',',@tables;
    print "table_names is $table_names\n";
    $dbVar->do(qq{DROP TABLES $table_names});
  }
}
