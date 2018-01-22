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
use ImportUtils qw(debug);
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use ImportUtils qw(dumpSQL debug create_and_load load);
use FindBin qw($Bin);

my ($TMP_DIR, $TMP_FILE, $LIMIT);
our ($alias_file, $genotype_file, $species);

GetOptions('tmpdir=s'  => \$ImportUtils::TMP_DIR,
           'tmpfile=s' => \$ImportUtils::TMP_FILE,
		   'alias_file=s' => \$alias_file,
		   'genotype_file=s' => \$genotype_file,
		   'species=s' => \$species,
	  );

warn("Make sure you have a updated ensembl.registry file!\n");

my $registry_file ||= $Bin . "/ensembl.registry";

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

$species ||= "rat";

Bio::EnsEMBL::Registry->load_all( $registry_file );

my $cdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $vdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');


my $dbVar = $vdba->dbc->db_handle;
my $dbCore = $cdba;

my (%rec,%rec_mdc, %rec_mdc_pop, %rec_mdc_ids, $population_ref_sample_id,$population_star_sample_id);


import_strain_name ();
create_db(\%rec_mdc,\%rec_mdc_pop);

sub import_strain_name {

  my $sym;
  
  die "Alias file not defined (use -alias_file [alias file])\n" unless defined $alias_file;
  
  debug("Reading alias file...");

  open IN, $alias_file or die "Could not read from alias file $alias_file\n";

  while (<IN>) {
    chomp;
    if (!/^\#/) {
      my ($mdc_id,$cng_id,$symbol,$country,$origin) = split /\t/, $_;
      #my ($symbol,$null1,$null2,$mdc_id) = split /\t/, $_;

      #print "symbol is $symbol and mdc_id is $mdc_id\n" if ($symbol and $mdc_id);
      next if !$symbol;

      $rec_mdc_ids{$mdc_id}=$symbol;
      my $sym_ids = "$symbol\|$mdc_id";
      push @{$rec{$symbol}}, "$sym_ids\t";
    }
  }
#   foreach my $sym (sort {$a cmp $b} keys %rec) {
#     my @ind = @{$rec{$sym}};
#     my $num_ind = scalar @ind;
#     print "\nsub_population/strain_name $sym, including $num_ind individuls: @ind\n" if $num_ind>1;
#   }

#  exit;
  my $ref_strain = "Ref:BN/Crl";
  my $star_pop_name = "STAR:mdc/cng";
  my $individual_type_id = 2;

  $population_ref_sample_id = get_id("sample",$ref_strain);
  if (!$population_ref_sample_id ) {
    debug("Insert into population and sample tables...");

    $dbVar->do(qq{INSERT INTO sample (name,description) values ("$ref_strain","reference population sample")});
    $population_ref_sample_id = $dbVar->{'mysql_insertid'};
    $dbVar->do(qq{INSERT INTO population (sample_id) values ($population_ref_sample_id)});
  }
  $population_star_sample_id = get_id("sample",$star_pop_name);
  if (!$population_star_sample_id ) {
    debug("Insert into population and sample tables...");

    $dbVar->do(qq{INSERT INTO sample (name,description) values ("$star_pop_name","STAR population samples")});
    $population_star_sample_id = $dbVar->{'mysql_insertid'};
    $dbVar->do(qq{INSERT INTO population (sample_id) values ($population_star_sample_id)});
  }

  foreach my $sym (sort {$a cmp $b} keys %rec) {
    my $population_strain_sample_id;
    my @ind = @{$rec{$sym}};
#     my $num_ind = scalar @ind;
#     my $pop_name = "MDC/CNG:$sym";
#     my $pop_desc = "a population used in STAR project, including $num_ind individuals";
#     #print "\nsub_population/strain_name $sym, including $num_ind individuls: @ind\n" if $sym;
#     $population_strain_sample_id = get_id("sample",$pop_name);
#     if (!$population_strain_sample_id) {
#       debug("Don't have population_strain_sample_id, insert it into sample table");
#       $dbVar->do(qq{INSERT INTO sample (name,size,description) values ("$pop_name",$num_ind,"$pop_desc")});
#       $population_strain_sample_id = $dbVar->{'mysql_insertid'};
#       $dbVar->do(qq{INSERT INTO population (sample_id) values ($population_strain_sample_id)});
#     }
    foreach my $ind (@ind) {
      my $strain_name = $sym;
      #print "strain_name is $strain_name\n";
      my ($symbol,$mdc_id) = split /\|/, $ind; $mdc_id =~ s/\s+//;
      my $ind_desc = "individual from STAR project $ind";
      my $individual_sample_id;
      $individual_sample_id = get_id("sample",$strain_name,$mdc_id);
      if (!$individual_sample_id) {
	debug("Don't have individual_sample_id, insert it into sample table");
	$dbVar->do(qq{INSERT INTO sample (name,description) values ("$strain_name","$ind_desc")});
	$individual_sample_id = $dbVar->{'mysql_insertid'};
	$dbVar->do(qq{INSERT INTO individual (sample_id,gender,individual_type_id) values ($individual_sample_id,"Unknown", $individual_type_id)});
	$dbVar->do(qq{INSERT INTO individual_population (individual_sample_id,population_sample_id) values ($individual_sample_id,$population_star_sample_id)}) if ($strain_name ne $ref_strain) ;
      }
      $rec_mdc{$mdc_id} = $individual_sample_id;
      $rec_mdc_pop{$mdc_id} = $population_star_sample_id;
      #print "$mdc_id and ind is $individual_sample_id and pop is $population_strain_sample_id\n" if ($mdc_id and $mdc_id =~ /MDC-03-86|MDC-03-16/);
    }
  }
}

sub get_id {

  my ($table,$name,$mdc_id) = @_;
  my ($id_ar,$id);
  print "table is $table and name is $name\n" if ($mdc_id and $mdc_id =~ /MDC-03-86|MDC-03-16/);
  if (!$mdc_id) {
    $id_ar = $dbVar->selectall_arrayref(qq{SELECT $table\_id from $table WHERE name = "$name"});
    $id = $id_ar->[0][0];
  }
  else {
    $id_ar = $dbVar->selectall_arrayref(qq{SELECT $table\_id from $table WHERE name = "$name" and description like "%$mdc_id%"});
    $id = $id_ar->[0][0];
  }
  #print "id is $id mdc_id is $mdc_id\n" if ($mdc_id and $mdc_id =~ /MDC-03-86|MDC-03-16/);
  return $id;
}

sub create_db {

  debug("Creating database now...");
  my (%rec_source,%rec_column,@individuals,%gtype_code);
  my @source_names = ("ENSEMBL:star_gtype");
  my $pre_variation_name = "ENSRNOSNP";

  foreach my $source (@source_names) {
    my $source_id;
    $source_id = get_id("source",$source);
    if (!$source_id) {
      $dbVar->do(qq{INSERT INTO source (name) values ("$source")});
      $source_id = $dbVar->{'mysql_insertid'};
    }
    $rec_source{$source} = $source_id;
  }

  my $gtype_file = $genotype_file;
  open IN, "$gtype_file" or die "gtype_file is not exist\n";
  open OUT, ">$gtype_file\.err";
  open VAR, ">$TMP_DIR/variation.star";
  open FEA, ">$TMP_DIR/variation_feature.star";
  open FLA, ">$TMP_DIR/flanking_sequence.star";
  open ALL, ">$TMP_DIR/allele.star";
  open GTY, ">$TMP_DIR/tmp_individual_genotype_single_bp.star";

  my $column = 0;
  my $column1 = 0;
  my $count = 0;
  my $count1 = 0;
  my $id_count=0;
  my (%rec_region_id,%done);

  my $sth = $dbCore->dbc->prepare(qq{SELECT sr.seq_region_id,sr.name
                      FROM   seq_region_attrib sra, attrib_type at, seq_region sr
                      WHERE sra.attrib_type_id=at.attrib_type_id 
                      AND at.code="toplevel"
                      AND sr.seq_region_id = sra.seq_region_id 
  });

  $sth->execute();

  while (my ($seq_region_id,$seq_region_name) = $sth->fetchrow_array()) {
    $rec_region_id{$seq_region_name} = $seq_region_id;
  }

  if ($gtype_file !~ /sub/) {
    $dbVar->do(qq{CREATE TABLE IF NOT EXISTS tmp_individual_genotype_single_bp (
                           variation_id int not null,allele_1 varchar(255),allele_2 varchar(255),sample_id int,
                           key variation_idx(variation_id),
                           key sample_idx(sample_id)
                           ) MAX_ROWS = 100000000
    });

    $dbVar->do(qq{CREATE UNIQUE INDEX uniq_idx ON allele(variation_id,allele,sample_id)});
    $dbVar->do(qq{CREATE UNIQUE INDEX uniq_idx ON tmp_individual_genotype_single_bp(variation_id,allele_1,allele_2,sample_id)});

    $dbVar->do(qq{ALTER TABLE variation add column internal_name varchar(50)});
  }

  my $last_var_id_ref = $dbVar->selectall_arrayref(qq{SELECT MAX(variation_id) from variation});
  my $last_var_id = $last_var_id_ref->[0][0];
  $count = $last_var_id;
  $count1 = $count;


  LINE : while (<IN>) {
    if (/^a1\_External\_ID/) {
      my @all = split /\t/, $_;
      @individuals = @all[7..$#all];
      foreach my $ind_id (@individuals) {
	$rec_column{$column} = $ind_id;
	$ind_id =~ s/\s+//;
	$column++;
        $id_count++;
	print " $id_count $ind_id not exist in aliase file\n" if !$rec_mdc_ids{$ind_id};
      }
      debug("Total column/individual is $column");
    }
    else {
      $count = $count1;
      $count++;
      my ($ext_id,$chr_name,$chr_pos,$flank_seq,$ref_allele,$strain_allele,$source,@gtypes) = split /\t/, $_;
      my $num_gtypes = scalar @gtypes;
      if ($num_gtypes != $column) {
	warn "The genotype_column is not equal total column $column for $chr_name/$chr_pos\n";
	next;
      }
      my ($allele1,$allele2) = split /\//, $strain_allele;
      $gtype_code{0} = "$ref_allele$ref_allele";
      $gtype_code{1} = "$ref_allele$allele2";
      $gtype_code{2} = "$allele2$allele2";
      $gtype_code{5} = "";
      $gtype_code{6} = "";
      #my $source_name = "STAR_GTYPE:$source";
      my $source_name = "ENSEMBL:star_gtype";
      my $source_id = $rec_source{$source_name};
      my $variation_name = "$pre_variation_name$count";
      my $internal_name = "$chr_name\_$chr_pos";
      my $flank5 = substr($flank_seq,0,40);
      my $flank3 = substr($flank_seq,41);

      my $variation_id;
      if (%rec_region_id) {
	my $variation_id_ar = $dbVar->selectall_arrayref(qq{SELECT variation_id FROM variation_feature
                                                        WHERE seq_region_id=$rec_region_id{$chr_name}
                                                        AND seq_region_start = $chr_pos
        });

	$variation_id = $variation_id_ar->[0][0];
	#print "variation_id is $variation_id\n" if $variation_id;
      }
      if (!$variation_id) {
	#print variation_feature table...;
	print FEA join ("\t", $rec_region_id{$chr_name},$chr_pos,$chr_pos,1,$count,"$strain_allele","$variation_name",0,'\N',$source_id,'\N',"INTERGENIC"),"\n";
	#print variation table...
	print VAR "$count\t$source_id\t$variation_name\t$internal_name\n";

	#debug("Insert into flanking_sequence table...");
	print FLA join ("\t", $count,$flank5,$flank3,$chr_pos-40,$chr_pos-1,$chr_pos+1,$chr_pos+40,$rec_region_id{$chr_name},1),"\n";
	$count1 = $count;
      }
      else {
	$count1 = $count-1;
	$count = $variation_id;
      }
      #debug("Insert into allele/gtype table...");
      foreach my $gtype_code (@gtypes) {
	if (!$gtype_code{$gtype_code}) {
	  $column1++;
	  if ($column1 == $column) {
	    $column1 = 0;
	    next LINE;
	  }
	  next;
	}
	else {
	  my $pop_sample_id = $rec_mdc_pop{$rec_column{$column1}} if $rec_mdc_pop{$rec_column{$column1}};
	  my $ind_sample_id = $rec_mdc{$rec_column{$column1}} if $rec_mdc{$rec_column{$column1}};
	  my ($allele_1,$allele_2) = split "", $gtype_code{$gtype_code};
	  if (! $pop_sample_id or ! $ind_sample_id) {
	    print OUT "column $column1 with mdc_id $rec_column{$column1} don't have pop/ind_ids\n";
	  }
	  else {
	    if (!$done{$count}) {#get rid of duplicates
	      print ALL join ("\t", $count,$ref_allele,$population_ref_sample_id),"\n";
	      $done{$count}=1;
	    }
	    print ALL join ("\t", $count,$allele_2,$pop_sample_id),"\n";
	    print GTY join ("\t", $count,$allele_1,$allele_2,$ind_sample_id),"\n";
	  }
	  $column1++;
	  if ($column1 == $column) {
	    $column1 = 0;
	    next LINE;
	  }
	}
      }
    }
  }

  #importing tables...
  foreach my $table ("variation","variation_feature","flanking_sequence","allele","tmp_individual_genotype_single_bp") {
    system("cp $TMP_DIR/$table\.star $TMP_DIR/$TMP_FILE");
    debug("Importing table $table...");

    load($dbVar,qw(variation variation_id source_id name internal_name)) if ($table eq "variation" and !(-z "$TMP_DIR/$TMP_FILE"));
    load($dbVar,qw(variation_feature seq_region_id seq_region_start seq_region_end seq_region_strand variation_id allele_string variation_name map_weight flags source_id validation_status consequence_type)) if ($table eq "variation_feature" and !(-z "$TMP_DIR/$TMP_FILE"));
    load($dbVar,qw(flanking_sequence variation_id up_seq down_seq up_seq_region_start up_seq_region_end down_seq_region_start down_seq_region_end seq_region_id seq_region_strand)) if ($table eq "flanking_sequence" and !(-z "$TMP_DIR/$TMP_FILE"));
    load($dbVar,qw(allele variation_id allele sample_id)) if ($table eq "allele" and !(-z "$TMP_DIR/$TMP_FILE"));
    load($dbVar,qw(tmp_individual_genotype_single_bp variation_id allele_1 allele_2 sample_id)) if ($table eq "tmp_individual_genotype_single_bp" and !(-z "$TMP_DIR/$TMP_FILE"));
 }

}
