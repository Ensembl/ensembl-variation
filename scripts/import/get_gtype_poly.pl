#!/usr/bin/env perl
# Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
  developers list at <dev@ensembl.org>.

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


GetOptions('species=s'          => \$species,
           'tmpdir=s'           => \$ImportUtils::TMP_DIR,
           'tmpfile=s'          => \$ImportUtils::TMP_FILE,
		   'registry=s'         => \$registry_file,
);

$registry_file ||= $Bin . "/ensembl.registry";
$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;
Bio::EnsEMBL::Registry->load_all( $registry_file );

my $cdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $vdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
my $sadb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,"sara_$species");

my $dbCore = $cdb->dbc;
my $dbVar = $vdb->dbc;
my $dbSara = $sadb->dbc;

my (%seq_region_ids);
 

my $sthc = $dbCore->prepare(qq{
                select sra.seq_region_id, sr.name 
                from seq_region_attrib sra, attrib_type at, seq_region sr
                where sra.attrib_type_id = at.attrib_type_id 
                and at.code = "toplevel" 
                and sr.seq_region_id = sra.seq_region_id });
$sthc->execute();

while (my ($seq_region_id, $seq_region_name) = $sthc->fetchrow_array()) {
  $seq_region_ids{$seq_region_id} = $seq_region_name;
}

get_gtype_poly();
#star_frequency(\%seq_region_ids);
#base_coverage(\%seq_region_ids);

sub get_gtype_poly {

  my @sample_ids;
  my $count=0;
  my %sample_id;
  my $sth;

  #get rid of duplicated sample names first
  my $name_ref = $dbVar->db_handle->selectall_arrayref(qq{select name from (select name, count(*) as count from individual group by name having count >1) as name_count});
  my @names = map{$_->[0]} @$name_ref;
  my $name_string = (@names) ? "WHERE name NOT IN (".join(',',map{"'$_'"} @names).");" : '';
  $sth = $dbVar->prepare(qq{SELECT individual_id from individual $name_string});

  $sth->execute();
  my ($sample_id,$name,%sample_name);
  $sth->bind_columns(\$sample_id,\$name);

  while ($sth->fetch()) {
    $name =~ s/\/| |\+|\.|\-|\,|\(|\)|\<|\>/\_/g;
    $sample_name{$sample_id} = $name;
    push @sample_ids, $sample_id;
  }

  %sample_id = map {$count++,$_} @sample_ids;

  #my @nums = keys %sample_id;
  #print "There are @nums ",scalar @nums," sample_ids\n";

  #foreach my $sample_id (keys %sample_name) {
  #  print "num is $num and sample_id is $sample_id{$num}\n";
  #  print "name is $sample_name{$sample_id}\n";
  #}

  dumpSQL($dbVar,qq{SELECT tg.variation_id,tg.allele_1,tg.individual_id,substring(vf.allele_string,1,1) as ref_allele
                    FROM tmp_individual_genotype_single_bp tg, variation_feature vf
                    WHERE tg.allele_1 = tg.allele_2
                    AND tg.variation_id = vf.variation_id
                    #LIMIT 10
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
    push @columns, "$sample_name{$sample_id{$num}} varchar(100)"; #column names
    push @cols, "$sample_name{$sample_id{$num}}";#column values
  }

  my $ref_strain;
  if ($dbVar->dbname =~ /rat|ratt/) {
    $ref_strain = "reference_BN_SsNHsd_Mcwi";
  }
  elsif ($dbVar->dbname =~ /mouse|mus/) {
     $ref_strain = "reference_C57BL_6J";
   }

  push @columns, "$ref_strain varchar(100)";
  push @cols, "$ref_strain";

  my $column_name = join ",",@columns;
  my $cols_name = join ",",@cols;
  #print "column_name is $column_name\n";

  $dbVar->do(qq{DROP TABLE IF EXISTS strain_gtype_poly});
  $dbVar->do(qq{CREATE TABLE IF NOT EXISTS strain_gtype_poly (
              variation_id int(10) unsigned not null,
              $column_name,
              primary key( variation_id ))
              });
  system("mv $TMP_DIR/$TMP_FILE\_GTYPE $TMP_DIR/$TMP_FILE");
  #print "$cols_name\n";
  load($dbVar,"strain_gtype_poly","variation_id",$cols_name);
}

sub star_frequency {

  my $seq_region_ids = shift;
  my %seq_region_ids = %$seq_region_ids;

  my $tra  = $cdb->get_TranscriptAdaptor();
  my $sa = $cdb->get_SliceAdaptor();

=head
  my (%rec_cons,%rec_sample,%count);
  my $ids= $dbVar->db_handle->selectall_arrayref(qq{SELECT variation_id from variation where source_id =3});
  my %variation_ids = map {$_->[0],1 } @$ids;
  my $ids_vs= $dbVar->db_handle->selectall_arrayref(qq{SELECT variation_id from variation_synonym where source_id =3});
  %variation_ids = map {$_->[0],1 } @$ids_vs;
  my ($variation_id,$consequence_type);
  my $sth = $dbVar->prepare(qq{SELECT vf.variation_id,vf.consequence_type 
                               FROM variation_feature vf
                               });
  $sth->execute();
  $sth->bind_columns(\$variation_id,\$consequence_type);
  while ($sth->fetch()) {
    if ($variation_ids{$variation_id}) {
      if ($consequence_type =~ /CODING/) {
	$consequence_type = "cds";
      }
      elsif ($consequence_type =~ /5PRIME_UTR/) {
	$consequence_type = "5PRIME_UTR";
      }
      elsif ($consequence_type =~ /3PRIME_UTR/) {
	$consequence_type = "3PRIME_UTR";
      }
      elsif ($consequence_type =~ /INTRONIC/) {
	$consequence_type = "INTRONIC";
      }
      
      $rec_cons{$variation_id}{$consequence_type}++;
    }
  }

  my $sth1 = $dbVar->prepare(qq{SELECT tg.variation_id,tg.sample_id
                               FROM tmp_individual_genotype_single_bp tg, sample s
                               WHERE s.sample_id = tg.sample_id
                               AND s.sample_id in (68,69,70,71)
                           });
  $sth1->execute();
  my ($var_id,$sample_id);
  $sth1->bind_columns(\$var_id,\$sample_id);
  while($sth1->fetch()) {
    $rec_sample{$sample_id}{$var_id}++;
  }

  foreach my $sample_id (keys %rec_sample) {
    foreach my $variation_id (keys %{$rec_sample{$sample_id}}) {
      foreach my $cons (keys %{$rec_cons{$variation_id}}) {
	$count{$sample_id}{$cons}++;
      }
    }
  }

  foreach my $sample_id (keys %count) {
    foreach my $cons (keys %{$count{$sample_id}}) {
      print "There are $count{$sample_id}{$cons} variations in $cons region in sample with $sample_id\n";
    }
  }

  undef %rec_cons;
  undef %rec_sample;
  undef %count;
=cut
  open OUT, ">$TMP_DIR/$TMP_FILE";
  foreach my $seq_region_id (keys %seq_region_ids) {
    #print "seq_region_id is $seq_region_id\n";
    my $slice = $sa->fetch_by_seq_region_id($seq_region_id);
    #print "slice is $slice\n";
    my @genes = sort {$a->start<=>$b->start} @{$slice->get_all_Genes_by_type('protein_coding','ensembl',1)};
    my ($n,$gene_start);
    my $gene_end=0;
    for ($n=0;$n<@genes;$n++) {
      print OUT join "\t",$seq_region_id,"intergenic",$gene_end+1,$genes[$n]->start-1,"\n" if ($gene_end+1 < $genes[$n]->start-1);
      $gene_end = $genes[$n]->end;
    }
    print OUT join "\t",$seq_region_id,"intergenic",$gene_end+1,$slice->end,"\n";

    foreach my $gene (@genes) {
      #print "gene_stable_id is ",$gene->stable_id,"\n";
      my @transcripts_start = sort {$a->start<=>$b->start} @{$gene->get_all_Transcripts};
      my @transcripts_end = sort {$b->end<=>$a->end} @transcripts_start;
      #print "transcript_stable_id is ",$transcript->stable_id,"\n";
      my @all_exons = sort {$a->seq_region_start<=>$b->seq_region_start} @{$gene->get_all_Exons()};
      my @all_introns;
#      for(my $i=0; $i < scalar(@all_exons)-1; $i++){
#        my $intron = new Bio::EnsEMBL::Intron($all_exons[$i],$all_exons[$i+1]);
#	 push(@all_introns, $intron);
#      }
      my $five_utr_start = $transcripts_start[0]->start();
      my $five_utr_end = $transcripts_start[0]->coding_region_start() - 1;
      my $three_utr_start = $transcripts_end[0]->coding_region_end() + 1;
      my $three_utr_end = $transcripts_end[0]->end();

      if (@all_introns) {
	foreach my $intron(@all_introns) {
	  my $seq_region_start = $intron->seq_region_start;
	  my $seq_region_end = $intron->seq_region_end;
	  ($seq_region_start,$seq_region_end) = ($seq_region_end,$seq_region_start) if $seq_region_start > $seq_region_end;
	  print OUT join "\t",$seq_region_id,"intron",$seq_region_start,$seq_region_end,"\n" if $seq_region_start < $seq_region_end;
	}
      }

      if (@all_exons) {
	foreach my $exon (@all_exons) {
	  my $seq_region_start = $exon->seq_region_start;
	  my $seq_region_end = $exon->seq_region_end;
	  ($seq_region_start,$seq_region_end) = ($seq_region_end,$seq_region_start) if $seq_region_start > $seq_region_end;
	  print OUT join "\t",$seq_region_id,"cds",$seq_region_start,$seq_region_end,"\n" if $seq_region_start < $seq_region_end;
	}
      }
      
      for (my $n=0;$n<@all_exons-1;$n++) {
        print OUT join "\t",$seq_region_id,"intron",$all_exons[$n]->seq_region_end+1,$all_exons[$n+1]->seq_region_start-1,"\n" if ($all_exons[$n]->seq_region_end+1 < $all_exons[$n+1]->seq_region_start-1);
      }

      print OUT join "\t",$seq_region_id,"three_utr",$three_utr_start,$three_utr_end,"\n" if ($three_utr_start < $three_utr_end);
      print OUT join "\t",$seq_region_id,"five_utr",$five_utr_start,$five_utr_end,"\n" if ($five_utr_start < $five_utr_end);
    }
    #exit;

  }

  #load($dbVar,"star_consequence_type_region","seq_region_id","consequence_type","start","end");
  create_and_load($dbVar,"star_consequence_type_region2","seq_region_id i*","consequence_type","start i*","end");

}

sub base_coverage {

  my $seq_region_ids = shift;
  my %seq_region_ids = %$seq_region_ids;
  my ($individual_name, $seq_region_id,$seq_region_start, $seq_region_end,%rec_ssaha);

  my $sth1 = $dbSara->prepare(qq{ select individual_name,target_seq_region_id,target_start, target_end 
                                 from ssahaSNP_feature 
                                 order by target_seq_region_id,target_start});
  $sth1->execute();
  $sth1->bind_columns(\$individual_name,\$seq_region_id,\$seq_region_start,\$seq_region_end);
  while ($sth1->fetch()) {
    push @{$rec_ssaha{$individual_name}{$seq_region_id}},"$seq_region_start-$seq_region_end";
  }

  my ($start,$end,$total_length,%rec_ind,%rec_type);
  foreach my $type ("intron","cds","three_utr","five_utr","intergenic") {
  #foreach my $type ("intron") {
  print "processing $type...\n";
    foreach my $individual (keys %rec_ssaha) {
    #foreach my $individual ("SS/Jr") {
      print "processing $individual...\n";
      foreach my $seq_region_id (keys %{$rec_ssaha{$individual}}) {
      #foreach my $seq_region_id (127175) {  
        #print "processing $seq_region_id...\n";
        dumpSQL($dbVar,qq{ select start,end from star_consequence_type_region2
        where seq_region_id = $seq_region_id 
        and consequence_type = "$type"
        #order by start
        } );
        `sort -k 1 -g -o $TMP_DIR/$TMP_FILE\_s $TMP_DIR/$TMP_FILE`;
        system("mv $TMP_DIR/$TMP_FILE\_s $TMP_DIR/$TMP_FILE");
        open IN, "$TMP_DIR/$TMP_FILE"; 
        while (<IN>) {
	  chomp;
	  my ($start,$end) = split;
	  my @regions = @{$rec_ssaha{$individual}{$seq_region_id}};
	  foreach my $region (@regions) {
	    my ($seq_region_start,$seq_region_end) = split /\-/, $region;
            last if $seq_region_start > $end ;
            next if $seq_region_end < $start ;
            #print "start is $start,end is $end,seq_region_start is $seq_region_start,seq_region_end is $seq_region_end\n";
	    if ($seq_region_start <= $start and $seq_region_end >= $end) {
	      my $length = $end-$start+1;
	      $total_length += $length;
	    }
	    elsif ($seq_region_start <= $start and $seq_region_end <= $end) {
	      my $length = $seq_region_end-$start+1;
	      $total_length += $length;
	    }
	    elsif ($seq_region_start >= $start and $seq_region_end <= $end) {
	      my $length = $seq_region_end-$seq_region_start+1;
	      $total_length += $length;
	    }
	    elsif ($seq_region_start >= $start and $seq_region_end >= $end) {
	      my $length = $end-$seq_region_start+1;
	      $total_length += $length;
	    }
	  }
        }
      }
      $rec_ind{$type}{$individual} = $total_length;
      undef $total_length;
    }
  }

  foreach my $type (keys %rec_ind) {
    foreach my $individual (keys %{$rec_ind{$type}}) {
      print "strain $individual has total $rec_ind{$type}{$individual} bases for $type\n";
    }
  }
}
