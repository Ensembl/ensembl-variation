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

#array files are download from affy websit, files are in directory:
#/ecs4/scratch2/yuan/array/Homo_sapiens,should be moved to work dir
#now moved to:/ecs4/work2/yuan/microarray/mapping_array
#mapping db is in:ecs2:3362:nath_homo_sapiens_core_40_36b
#var db is in :ia64g:yuan_affy_var_40
#affy_probeset_name<->rsids are in :ia64g:yuan_hum_aff_40

#&make_array_input;
&make_variation_db_file;

sub make_array_input {

  ##read input file:Mapping50K_Hind_probe_tab
  #SNP_A-1742900 1004 51 -4 TGCAGCTATGACCCAAAATTTGATG r PM A 
  #SNP_A-1742900 176  19 -4 CATCAAATGTTGGGTCATAGCTGCA f PM C
  my (%rec_id,%rec_allele);

  while (<Mapping*probe_id>) {
    print "file is $_\n";
    my $file_name = $_;
    my $out_file_name = "$_\.fa";
    #open OUT, ">$out_file_name" or die "can't open $out_file_name : $!";
    #open OUT1, ">>all_file_dir";
    open OUT2, ">>all_affy_alleles_new1";

    open IN, "$file_name" or die "can't open $file_name : $!";
       
    while (<IN>) {
      if (/^SNP/) {
	chomp;
	my $probe_set_id = $_;
	#print "probe_set_id is #$probe_set_id#\n";
	$rec_id{$probe_set_id}=1;
      }
    }
    my $seq_file = $file_name;
    my $array_name = $file_name;
    $seq_file =~s/id$/tab/;
    $array_name =~ s/\_probe\_id$//;
    print "seq_file name is $seq_file\n";
    open SEQ, "$seq_file" or die "can't open $seq_file: $!";
    while (<SEQ>) {
      my ($probe_set,$x,$y,$num,$seq,$dir,$match,$base) = split;
      #print "probe_set is $probe_set\n";
      if ($rec_id{$probe_set}) {
	#my $allele;
	#if ($dir eq 'f') {
	#  $allele = substr($seq,12+$num,1);
	#}
	#elsif ($dir eq 'r') {
	#  $allele = $base;
	#}
	my $full_name = "$array_name\_$probe_set";
	#$rec_allele{$full_name}{$allele}++;
	$rec_allele{$full_name}{$base}++;
	#print "full_name is $full_name and allele is $allele\n";
	#my $name = "$array_name\:$probe_set\:$x\:$y\-$num\-$base\-$dir\;";
	#print OUT ">$name\n$seq\n";
	#print OUT1 "$array_name\t$probe_set\t$x\:$y\-$num\-$base\;\t$num\t$base\t$dir\n";
      }
    }
  }
  #check alleles only have 2 different base
  foreach my $full_name (keys %rec_allele) {
    my @alleles = keys %{$rec_allele{$full_name}};
    if (scalar @alleles ==2) {
      my $allele_string = join "/",@alleles;
      print OUT2 "$full_name\t$allele_string\n";
    }
  }
}

sub make_variation_db_file {
  
  #aff_snp combine info from olig_feature,oligo_probe,oligo_array and the one I generated which contain direction/base info which was missed when I gave data to Nathan.
  open IN, "aff_snp";
  open IN1, "all_allele";
  open VAR, ">variation.aff_new";
  open VF, ">variation_feature.aff_new";
  #open ALLELE, ">allele.aff";
  open FL, ">flanking_sequence.aff_new";

  my (%rec,%rec_base,$variation_id,$variation_feature_id,$variation_synonym_id,$new_seq_region_start,$new_base);
  my $source_id=4;

  while (<IN1>) {
    my ($array_name,$probeset,$allele_string) = split;
    my $var_name = "$array_name\_$probeset";
    $rec_base{$var_name}=$allele_string;
  }

  while (<IN>) {
    my ($array_name,$probeset,$null,$num,$base,$dir,$null1,$seq_region_id,$seq_region_start,$seq_region_end,$seq_region_strand) = split;
    my $var_name = "$array_name\_$probeset";
    if ($dir eq 'f') {
      $new_seq_region_start = $seq_region_start +12 + $num;
    }
    elsif ($dir eq 'r') {
      $new_seq_region_start = $seq_region_start +12 - $num;
    }
    #if ($probeset eq "SNP_A-2038704") {print "$_ AND $old_base $new_base\n";}
    my $pos = "$seq_region_id\_$new_seq_region_start";
    $rec{$var_name}{$pos}++;
  }

  foreach my $var_name (keys %rec) {
    my @start_keys = keys %{$rec{$var_name}};
    if (@start_keys ==1 ) {
      my $new_id_start = $start_keys[0];
      my ($seq_region_id,$new_start) = split /\_/, $new_id_start;
      my $seq_region_strand=1; #make seq_region_strand alwayls=1 because $base is always in forward direction
	my $allele_string = $rec_base{$var_name};
	print "var_name is $var_name and allele_string is $allele_string\n";
	$variation_id++;
	$variation_feature_id++;
	$variation_synonym_id++;
	my $map_weight = 1;
	my $up_seq_region_start = $new_start -100;
	my $up_seq_region_end = $new_start-1;
	my $down_seq_region_start = $new_start + 1;
	my $down_seq_region_end = $new_end + 100;
	print VAR "$variation_id\t$source_id\t$var_name\t\\N\t\\N\n";
	print VF "$variation_feature_id\t$seq_region_id\t$new_start\t$new_start\t$seq_region_strand\t$variation_id\t$allele_string\t$var_name\t$map_weight\t\\N\t$source_id\t\\N\t\\N\n";
	print FL "$variation_id\t\\N\t\\N\t$up_seq_region_start\t$up_seq_region_end\t$down_seq_region_start\t$down_seq_region_end\t$seq_region_id\t$seq_region_strand\n";
    }
  }
}

