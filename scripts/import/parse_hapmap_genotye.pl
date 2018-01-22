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

# this script read hapmap genotyping data into individual, 
#individual_genotype, population, population_genotype tables

use strict;
use warnings;
use DBI;
use DBH;
use Getopt::Long;
use Benchmark;
use Bio::EnsEMBL::DBSQL::DBConnection;
use ImportUtils qw(dumpSQL debug create_and_load load );

my ($TAX_ID, $LIMIT_SQL, $CONTIG_SQL, $TMP_DIR, $TMP_FILE);

my $dbSNP;
my $dbVar;
my $dbCore;
my %rec;

{
  my($dshost, $dsuser, $dspass, $dsport, $dsdbname, # dbSNP db
     $chost, $cuser, $cpass, $cport, $cdbname,      # ensembl core db
     $vhost, $vuser, $vpass, $vport, $vdbname,      # ensembl variation db
     $genotype_file, $limit);
  
  GetOptions('vuser=s'   => \$vuser,
             'vhost=s'   => \$vhost,
             'vpass=s'   => \$vpass,
             'vport=i'   => \$vport,
             'vdbname=s' => \$vdbname,
	     'chost=s'   => \$chost,
             'cuser=s'   => \$cuser,
             'cpass=s'   => \$cpass,
             'cport=i'   => \$cport,
             'cdbname=s' => \$cdbname,
	     'genotypefile=s' => \$genotype_file,
             'tmpdir=s'  => \$ImportUtils::TMP_DIR,
             'tmpfile=s' => \$ImportUtils::TMP_FILE,
             'limit=i'   => \$limit);
  
  $vhost    ||= 'ia64g';
  $vport    ||= 3306;
  $vuser    ||= 'ensadmin';
  $vdbname  ||= 'dr2_homo_sapiens_variation_27';

  $chost    ||= 'ecs2';
  $cuser    ||= 'ensro';
  $cport    ||= 3364;
  $cdbname  ||= 'homo_sapiens_core_27_35a';

  usage('-cdbname argument is required.') if(!$cdbname);
  usage('-vdbname argument is required.') if(!$vdbname);
  usage('-genotypefile argument is required.') if(!$genotype_file);
  usage('-vpass argument is required.') if(!$vpass);

  $TMP_DIR  = $ImportUtils::TMP_DIR;
  $TMP_FILE = $ImportUtils::TMP_FILE;


  $LIMIT_SQL = ($limit) ? " LIMIT $limit " : '';

  $dbVar = DBH->connect
    ("DBI:mysql:host=$vhost;dbname=$vdbname;port=$vport",$vuser, $vpass,
    {'RaiseError' => 1});
  die("Could not connect to variation database: $!") if(!$dbVar);

  $dbCore = Bio::EnsEMBL::DBSQL::DBConnection->new
    (-host   => $chost,
     -user   => $cuser,
     -pass   => $cpass,
     -port   => $cport,
     -dbname => $cdbname);


  import($genotype_file);
}

sub import {

  my $genotype_file = shift;

  our ($population_name, %pop_desc, %ind_name, %ind_snp, %var_ids);

  if ($genotype_file =~ /\_JPT\.txt/) {
    $population_name = "Japanese";
    $pop_desc{$population_name}{'desc'} = "Hapmap_JPT: Japanese in Tokyo, Japan";
    $pop_desc{$population_name}{'super_pop_id'} = 6;
  }
  elsif ($genotype_file =~ /\_YRI\.txt/) {
    $population_name = "Yoruba-30-trios";
    $pop_desc{$population_name}{'desc'} = "Hapmap_YRI: Yoruba in Ibadan, Nigeria";
    $pop_desc{$population_name}{'super_pop_id'} = 4;
  }
  elsif ($genotype_file =~ /\_CEU\.txt/) {
    $population_name = "CEPH-30-trios";
    $pop_desc{$population_name}{'desc'} = "Hapmap_CEU: CEPH (Utah residents with ancestry from northern and western Europe)";
    $pop_desc{$population_name}{'super_pop_id'} = 5;
  }
  elsif ($genotype_file =~ /\_HCB\.txt/) {
    $population_name = "Han_Chinese";
    $pop_desc{$population_name}{'desc'} = "Hapmap_HCB: Han Chinese in Beijing, China";
    $pop_desc{$population_name}{'super_pop_id'} = 6;
  }

  if ($genotype_file =~ /\.gz$/) {
    open GENO, "gunzip -c $genotype_file |" or debug"can't open file $genotype_file: $!\n";
  }
  else {
    open GENO, "$genotype_file |" or debug"can't open file $genotype_file: $!\n";
  }

  
  while (<GENO>) {
    chomp;
    if (/^rs\#/) {
      my @array = split;
      my $count = 1; 
      my @a= splice (@array,11);
      ###there are 11 columns before the individual names
      %ind_name = map {$count++, $_} @a; 
    }
    elsif (/^rs\d+\s+/) {
      my %ind_gty;
      my @array = split;
      my $snp_id= $array[0];
      my $count = 1;
      my @a = splice (@array,11);
      %ind_gty = map {$count++, $_} @a;
      $ind_snp{$snp_id} = \%ind_gty;
    }

  }

  ##to make a hash hold snp_id and variation_id
  my @snp_ids = keys %ind_snp;
  my $snp_id_str = "IN (" .join(',', map{"'$_'"} @snp_ids).");";
  
  my $sth = $dbVar->prepare (qq{SELECT variation_id, name
                                FROM variation where name $snp_id_str});
  $sth->execute();

  while(my ($variation_id, $name) = $sth->fetchrow_array()) {
    $var_ids{$name} = $variation_id;
  }

  my @names = values %ind_name;
  my $ind_name_str = "IN (".join(',', map{"'$_'"} @names). ");";

  my $sub_pop_id = population(\%pop_desc);
  population_genotype(\%ind_snp, \%var_ids, $sub_pop_id);
  individual(\%ind_name, $sub_pop_id, $ind_name_str);
  individual_genotype(\%ind_name, \%ind_snp, \%var_ids, $sub_pop_id, $ind_name_str);

}


sub population {

  my $pop_desc = shift;

  my %pop_desc = %$pop_desc;
  my ($name) = keys %pop_desc;
  my $desc = $pop_desc{$name}->{'desc'};
  my $super_pop_id = $pop_desc{$name}->{'super_pop_id'};

  ##checking the name is already in population table or not
  my $sth = $dbVar->prepare (qq{SELECT population_id from population where name = ?;});
  $sth->execute("$name");
  
  my ($sub_pop_id) = $sth->fetchrow_array;

  if (!$sub_pop_id) {
    debug("loading data into population and population_structure tables...");

    $dbVar->do(qq{INSERT INTO population (name,description) values ("$name","$desc")
		 });
    $sub_pop_id = $dbVar->dbh()->{'mysql_insertid'};
    $dbVar->do(qq{INSERT INTO population_structure (super_population_id, sub_population_id)
		values ($super_pop_id,$sub_pop_id)
	       });
  }

  return $sub_pop_id;
}

sub population_genotype {

  my ($ind_snp, $var_ids, $sub_pop_id) = @_;

  my %ind_snp = %$ind_snp;
  my %var_ids = %$var_ids;
  my $TMP_FILE = "population_genotype.$$";

  open FH, ">$TMP_DIR/$TMP_FILE" or die "can't open file $TMP_DIR/$TMP_FILE\n";
  
  foreach my $snp_id (keys %ind_snp) {
    if ($snp_id and $var_ids{$snp_id}) {
      my %count;
      my %ind_gty = %{$ind_snp{$snp_id}};
      
      
      my $tot_count = keys %ind_gty;
      foreach my $num (keys %ind_gty) {
	$count{$ind_gty{$num}}++;
      }
      foreach my $gtys (keys %count) {
	my $tot_gty = keys %count;
	my ($allele1, $allele2) = split "", $gtys;
	my $freq = $count{$gtys}/$tot_count if ($tot_count != 0);
	#print "$var_ids{$snp_id}\t$allele1\t$allele2\t$freq\t$sub_pop_id\n";
	print FH "$var_ids{$snp_id}\t$allele1\t$allele2\t$freq\t$sub_pop_id\n";
      }
    }
  }
  
  close FH; ###need to close handle before loading data

  #debug("loading data into population_genotype table...\n");

  #load( $dbVar, "population_genotype", "variation_id", "allele_1", "allele_2", "frequency", "population_id" );
}
    
sub individual {

  my ($ind_name, $sub_pop_id, $ind_name_str) = @_;
  
  my %ind_name = %$ind_name;
  my (%ind_ids, $new_ind_found);
  #my $TMP_FILE = "individual.$$";

  open ( FH, ">$TMP_DIR/$TMP_FILE" ) or die "can't open file $TMP_DIR/$TMP_FILE\n";
  
  my $sth = $dbVar->prepare (qq{SELECT individual_id, name
                                FROM individual where population_id = $sub_pop_id
				AND name $ind_name_str});
  $sth->execute();

  while(my ($ind_id, $name) = $sth->fetchrow_array()) {
    $ind_ids{$name} = $ind_id if $ind_id;
  }
  
  foreach my $num (keys %ind_name) {
    if (!$ind_ids{$ind_name{$num}}) {
      $new_ind_found=1;
      print FH "$ind_name{$num}\t$ind_name{$num}\t$sub_pop_id\tUnknown\n";
    }
  }

  close FH; ###need to close handle before loading data

  if ($new_ind_found) {
    debug("loading data into individual table...\n");

    load($dbVar,"individual","name","description","population_id",
	 "gender","father_individual_id","mother_individual_id");
  }
}

sub individual_genotype {

  my ($ind_name, $ind_snp, $var_ids, $sub_pop_id, $ind_name_str) = @_;

  my %ind_name = %$ind_name;
  my %ind_snp = %$ind_snp;
  my %var_ids = %$var_ids;
  my %ind_ids;
  my $TMP_FILE = "individual_genotype.$$";

  open ( FH, ">$TMP_DIR/$TMP_FILE" ) or die "can't open file $TMP_DIR/$TMP_FILE\n";
  
  my $sth = $dbVar->prepare (qq{SELECT individual_id, name
                                FROM individual where population_id = $sub_pop_id
				AND name $ind_name_str});
  $sth->execute();

  while(my ($ind_id, $name) = $sth->fetchrow_array()) {
    $ind_ids{$name} = $ind_id if $ind_id; 
  }

  foreach my $snp_id (keys %ind_snp) {
    if ($snp_id and $var_ids{$snp_id}) {
      my %ind_gty = %{$ind_snp{$snp_id}};
      foreach my $num (keys %ind_gty) {
	my ($allele_1, $allele_2) = split "", $ind_gty{$num};
	print FH "$var_ids{$snp_id}\t$allele_1\t$allele_2\t$ind_ids{$ind_name{$num}}\n";
      }
    }
  }
  
  close FH;
  
  #debug("loading data into individual_genotype table...\n");
  
  #load($dbVar,"individual_genotype","variation_id","allele_1","allele_2","individual_id");
}

sub usage {
  my $msg = shift;

  print STDERR <<EOF;

usage: parse_hapmap_genotype.pl <options>

options:
    -chost <hostname>    hostname of core Ensembl MySQL database (default = ecs2)
    -cuser <user>        username of core Ensembl MySQL database (default = ensro)
    -cpass <pass>        password of core Ensembl MySQL database
    -cport <port>        TCP port of core Ensembl MySQL database (default = 3364)
    -cdbname <dbname>    dbname of core Ensembl MySQL database
    -vhost <hostname>    hostname of variation MySQL database to write to
    -vuser <user>        username of variation MySQL database to write to (default = ensadmin)
    -vpass <pass>        password of variation MySQL database to write to
    -vport <port>        TCP port of variation MySQL database to write to (default = 3306)
    -vdbname <dbname>    dbname of variation MySQL database to write to
    -limit <number>      limit the number of rows for testing
    -tmpdir <dir>        temp directory to use (with lots of space!)
    -tmpfile <filename>   name of temp file to use
    -num_processes <number> number of processes that are running (default = 1)
    -status_file <filename> name of a temp file where all the processes write when they finish
EOF

  die("\n$msg\n\n");
}

