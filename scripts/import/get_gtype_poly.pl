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


#directly use this script if generate from small genome like zebrafish

use strict;
#use DBH;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Data::Dumper;
use Bio::SeqIO;
use FindBin qw( $Bin );
use Getopt::Long;
use ImportUtils qw(dumpSQL debug create_and_load load);

our ($species,$TMP_DIR,$TMP_FILE,$registry_file);

GetOptions(
  'species=s'  => \$species,
  'tmpdir=s'   => \$ImportUtils::TMP_DIR,
  'tmpfile=s'  => \$ImportUtils::TMP_FILE,
  'registry=s' => \$registry_file,
);

#$registry_file ||= $Bin . "/ensembl.registry";
$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;
Bio::EnsEMBL::Registry->load_all( $registry_file );

my $vdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');

my $dbVar = $vdb->dbc;

get_gtype_poly();

sub get_gtype_poly {
  my @sample_ids;
  my $count=0;
  my %sample_id;
  my $sth;

  my $sth = $dbVar->prepare(qq{select distinct(s.sample_id), s.name from tmp_sample_genotype_single_bp t, sample s where s.sample_id = t.sample_id;});

  $sth->execute();
  my ($sample_id,$name,%sample_name);
  $sth->bind_columns(\$sample_id,\$name);

  while ($sth->fetch()) {
    $name =~ s/\/| |\+|\.|\-|\,|\(|\)|\<|\>/\_/g;
    $sample_name{$sample_id} = $name;
    print STDERR "$sample_id $name\n";
    push @sample_ids, $sample_id;
  }

  %sample_id = map {$count++,$_} @sample_ids;
  dumpSQL($dbVar, qq{
    SELECT tg.variation_id,tg.allele_1,tg.sample_id,substring(vf.allele_string,1,1) as ref_allele
    FROM tmp_sample_genotype_single_bp tg, variation_feature vf
    WHERE tg.allele_1 = tg.allele_2
    AND tg.variation_id = vf.variation_id
    });
  
  `sort -k 1 -g -o $TMP_DIR/$TMP_FILE\_s $TMP_DIR/$TMP_FILE`;
  system("mv $TMP_DIR/$TMP_FILE\_s $TMP_DIR/$TMP_FILE");

  open IN, "$TMP_DIR/$TMP_FILE";
  open OUT, ">$TMP_DIR/$TMP_FILE\_GTYPE";

  my ($previous_varid,$previous_ref_allele,@alleles,%rec,%done);

  while (<IN>) {
    my ($variation_id,$allele,$sample_id,$ref_allele) = split;
    next if (!$sample_name{$sample_id});

    if ($previous_varid != $variation_id and $previous_varid != 0) {
      foreach my $num (sort {$a<=>$b} keys %sample_id) {
        my $sampleid = $sample_id{$num};
        my $allele = ($rec{$sampleid}) ? $rec{$sampleid} : '';
        push @alleles, $allele;
      }
      push @alleles, $previous_ref_allele;
      print OUT join ("\t", $previous_varid,@alleles)."\n";
      undef %rec;
      undef %done;
      undef @alleles;
    }

    $previous_varid = $variation_id;
    $previous_ref_allele = $ref_allele;
    if(!$done{$sample_id}) {
      $rec{$sample_id}=$allele;
      $done{$sample_id}=1;
    }
    else {#for more than one allele per sample_id, give it as "?"
      $rec{$sample_id} = "?" if $allele ne $rec{$sample_id};
    }
  }

  #for last variation_id
  foreach my $num (sort {$a<=>$b} keys %sample_id) {
    my $sampleid = $sample_id{$num};
    my $allele = ($rec{$sampleid}) ? $rec{$sampleid} : '';
    push @alleles, $allele;
  }
  push @alleles, $previous_ref_allele;
  print OUT join ("\t", $previous_varid,@alleles)."\n";

  my (@cols,@columns);
  foreach my $num (sort {$a<=>$b} keys %sample_id ) {
    my $strain_name = $sample_name{$sample_id{$num}};
    push @columns, "`$strain_name` varchar(100)"; #column names
    push @cols, "`$strain_name`"; #column names
  }

  my $ref_strain;
  if ($dbVar->dbname =~ /rat|ratt/) {
    $ref_strain = "reference_BN_SsNHsd_Mcwi";
  }
  elsif ($dbVar->dbname =~ /mouse|mus/) {
    $ref_strain = "reference_C57BL_6J";
  }

  push @columns, "`$ref_strain` varchar(100)";
  push @cols, "`$ref_strain`";

  my $column_name = join ",",@columns;
  my $cols_name = join ",",@cols;

  $dbVar->do(qq{DROP TABLE IF EXISTS strain_gtype_poly});
  $dbVar->do(qq{CREATE TABLE IF NOT EXISTS strain_gtype_poly (
                variation_id int(10) unsigned not null,
                $column_name,
                primary key( variation_id ))
             });
  system("mv $TMP_DIR/$TMP_FILE\_GTYPE $TMP_DIR/$TMP_FILE");
  load($dbVar,"strain_gtype_poly","variation_id",$cols_name);
}

