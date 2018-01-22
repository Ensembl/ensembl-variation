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

our ($species,$TMP_DIR,$TMP_FILE);


GetOptions('species=s'          => \$species,
           'tmpdir=s'           => \$ImportUtils::TMP_DIR,
           'tmpfile=s'          => \$ImportUtils::TMP_FILE,
);

my $registry_file ||= $Bin . "/ensembl.registry";
$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;
Bio::EnsEMBL::Registry->load_all( $registry_file );

my $cdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $vdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');

my $dbCore = $cdb->dbc;
my $dbVar = $vdb->dbc;

my (%seq_region_ids);

=head
$dbVar->do(qq{CREATE TABLE tmp_gtype (variation_id int,allele varchar(5),unique index uniq_inx(variation_id,allele))});
$dbVar->do(qq{INSERT IGNORE INTO tmp_gtype select variation_id,allele_1 as allele from tmp_individual_genotype_single_bp_venter_uniq});
$dbVar->do(qq{INSERT IGNORE INTO tmp_gtype select variation_id,allele_2 as allele from tmp_individual_genotype_single_bp_venter_uniq});
$dbVar->do(qq{create table tmp_gtype_uniq like tmp_individual_genotype_single_bp});
$dbVar->do(qq{create table varid_one_allele_tmp_gtype select variation_id, sample_id,count(*) as count from tmp_gtype group by variation_id,sample_id having count=1});
$dbVar->do(qq{alter table varid_one_allele_tmp_gtype add index variation_idx(variation_id,sample_id)});
$dbVar->do(qq{insert into tmp_gtype_uniq select t.variation_id,t.allele as allele_1,t.allele as allele_2,t.sample_id from tmp_gtype t, varid_one_allele_tmp_gtype v where t.variation_id=v.variation_id and t.sample_id=v.sample_id});
$dbVar->do(qq{insert into tmp_gtype_uniq select t.variation_id,substring(group_concat(t.allele),1,1) as allele_1,substring(group_concat(t.allele),3,1) as allele_2,t.sample_id from tmp_gtype t group by t.variation_id,sample_id having count(*) =2});
$dbVar->do(qq{create table tmp_two_alleles select variation_id,count(*) as count from tmp_gtype_uniq group by variation_id having count=2});
$dbVar->do(qq{create table tmp_two_allele_with_refsample select t.variation_id from tmp_two_alleles t, tmp_gtype_uniq t1 where t1.sample_id=20 and t.variation_id=t1.variation_id});
$dbVar->do(qq{alter table tmp_two_allele_with_refsample add index variation_id(variation_id)});
=cut
my $sth = $dbVar->prepare(qq{select t.*, vf.allele_string 
                             from tmp_gtype_uniq t, tmp_two_allele_with_refsample t1, variation_feature vf
                             where t.variation_id=t1.variation_id
                             and vf.variation_id = t.variation_id 
                             order by t1.variation_id},{mysql_use_result => 1});
$sth->execute();
my ($previous_varid, $previous_refallele,%rec);

open OUT, ">$TMP_DIR/VARID_OK";
open OUT1, ">$TMP_DIR/tmp_gtype_uniq.delete";
open OUT2, ">$TMP_DIR/tmp_gtype_uniq.import";

while (my ($variation_id, $allele_1, $allele_2,$sample_id,$allele_string) = $sth->fetchrow_array()) {
  my ($ref_allele) = split /\//, $allele_string;
  #print "$variation_id, $allele_1, $allele_2,$sample_id,$allele_string\n";

  if ($variation_id != $previous_varid and $previous_varid) {
    foreach my $varid (keys %rec) {
      my ($sample_id1,$sample_id2) = keys %{$rec{$varid}};
      my ($al1a,$al2a) = split / /,$rec{$varid}{$sample_id1};
      my ($al1b,$al2b) = split / /,$rec{$varid}{$sample_id2};
      #print "$varid,$previous_refallele,$al1a,$al2a,$sample_id1,$al1b,$al2b,$sample_id2\n";
      if (($al1a =~/$previous_refallele/ or $al2a =~/$previous_refallele/) and $sample_id1 ==20) {
	print OUT "$varid ok\n";
      }
      elsif(($al1b =~/$previous_refallele/ or $al2b =~/$previous_refallele/) and $sample_id2 ==20) {
	print OUT "$varid\n";
      }
      elsif(($al1a =~/$previous_refallele/ or $al2a =~/$previous_refallele/) and $sample_id2 ==20) {
	print OUT1 "delete from tmp_gtype_uniq where variation_id = $varid;\n";
	print OUT2 "$varid\t$al1a\t$al2a\t$sample_id2\n";
	print OUT2 "$varid\t$al1b\t$al2b\t$sample_id1\n";
      }
      elsif (($al1b =~/$previous_refallele/ or $al2b =~/$previous_refallele/) and $sample_id1 ==20) {
	print OUT1 "delete from tmp_gtype_uniq where variation_id = $varid;\n";
	print OUT2 "$varid\t$al1a\t$al2a\t$sample_id2\n";
	print OUT2 "$varid\t$al1b\t$al2b\t$sample_id1\n";
      }
      else {
	print OUT "$varid can't be changed\n";
      }
      undef %rec;
      undef $previous_refallele;
    }
  }
  $previous_varid = $variation_id;
  $previous_refallele = $ref_allele;
  $rec{$previous_varid}{$sample_id} = "$allele_1 $allele_2";
}
#for last variation_id
foreach my $varid (keys %rec) {
  my ($sample_id1,$sample_id2) = keys %{$rec{$varid}};
  my ($al1a,$al2a) = split / /,$rec{$varid}{$sample_id1};
  my ($al1b,$al2b) = split / /,$rec{$varid}{$sample_id2};
  #print "$varid,$previous_refallele,$al1a,$al2a,$sample_id1,$al1b,$al2b,$sample_id2\n";
  if (($al1a =~/$previous_refallele/ or $al2a =~/$previous_refallele/) and $sample_id1 ==20) {
    print OUT "$varid\n";
  }
  elsif(($al1b =~/$previous_refallele/ or $al2b =~/$previous_refallele/) and $sample_id2 ==20) {
    print OUT "$varid\n";
  }
  elsif(($al1a =~/$previous_refallele/ or $al2a =~/$previous_refallele/) and $sample_id2 ==20) {
    print OUT1 "delete from tmp_gtype_uniq where variation_id = $varid;\n";
    print OUT2 "$varid\t$al1a\t$al2a\t$sample_id2\n";
    print OUT2 "$varid\t$al1b\t$al2b\t$sample_id1\n";
  }
  elsif (($al1b =~/$previous_refallele/ or $al2b =~/$previous_refallele/) and $sample_id1 ==20) {
    print OUT1 "delete from tmp_gtype_uniq where variation_id = $varid;\n";
    print OUT2 "$varid\t$al1a\t$al2a\t$sample_id2\n";
    print OUT2 "$varid\t$al1b\t$al2b\t$sample_id1\n";
  }
}



