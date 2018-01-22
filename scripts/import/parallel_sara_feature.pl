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

use strict;


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

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

our ($TMP_DIR, $TMP_FILE, $SEQ_REGION_ID, $SEQ_REGION_NAME, $species, $VAR_DBNAME, $JOB, $alignment_file,%CACHE,$index_filename);

{

  GetOptions('species=s'   => \$species,
	     'tmpdir=s'  => \$ImportUtils::TMP_DIR,
	     'tmpfile=s' => \$ImportUtils::TMP_FILE,
	     'seq_region_id=i'   => \$SEQ_REGION_ID,
	     'seq_region_name=s' => \$SEQ_REGION_NAME,
	     'job=s'     => \$JOB,
	     'alignment_file=s' =>\$alignment_file,
		 'index_file=s' => \$index_filename,
	     );


  my $registry_file ||= $Bin . "/ensembl.registry";

  Bio::EnsEMBL::Registry->load_all( $registry_file );

  my $cdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
  my $vdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
  my $vdb_sara = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'sara_hum');

  my $dbVar = $vdba->dbc;
  #my $dbCore = $cdba;
  my $dbSara = $vdb_sara->dbc;

  #we need to use cache to speed up fetch read sequence and quality data
  my $MAX_READS_CACHE_SIZE = 500;
  tie(%CACHE, 'Bio::EnsEMBL::Utils::Cache', $MAX_READS_CACHE_SIZE);
  
  $VAR_DBNAME = $dbVar->dbname();
  $TMP_DIR  = $ImportUtils::TMP_DIR;
  $ImportUtils::TMP_FILE .= "\_$SEQ_REGION_ID";
  $TMP_FILE = $ImportUtils::TMP_FILE;

  ssahasnp_feature($dbSara) if ($JOB eq "snp_pos");
  flanking_qual($dbSara) if ($JOB eq "flank");
  gtype_allele($dbSara) if ($JOB eq "gtype_allele");
  read_coverage() if ($JOB eq "read_coverage");
}

sub ssahasnp_feature {

  my ($dbSara) = @_;

  debug("Make query to variation_feature table");
  my $sth = $dbSara->prepare(qq{SELECT v.variation_id,vf.seq_region_start,vf.allele_string 
                                FROM  $VAR_DBNAME.variation v, $VAR_DBNAME.variation_feature vf
                                WHERE v.variation_id = vf.variation_id
                                AND vf.seq_region_id = $SEQ_REGION_ID
                                #ORDER BY seq_region_start
                               });
  my ($variation_id,$seq_region_id,$seq_region_start,$allele_string,%rec_var,%rec_seq_start);
  $sth->execute();
  $sth->bind_columns(\$variation_id,\$seq_region_start,\$allele_string);

  while ($sth->fetch()) {
    $rec_var{$variation_id} = $allele_string;
    $rec_seq_start{$seq_region_start} = $variation_id;
  }

  combine_feature($dbSara,\%rec_var,\%rec_seq_start);
  snp_pos($dbSara);

}

sub combine_feature {
  my $dbSara = shift;
  my $rec_var= shift;
  my $rec_seq_start = shift;
  my %rec_var = %$rec_var;
  my %rec_seq_start = %$rec_seq_start;
  #sort seq_start to speed up two big array comparison
  my @sorted_start = sort {$a<=>$b} keys %rec_seq_start;
  my @seq_region_ids;

  debug("Dumping ssahaSNP_feature");
  dumpSQL($dbSara, qq{SELECT * FROM ssahaSNP_feature WHERE target_seq_region_id = $SEQ_REGION_ID ORDER BY target_start});

  open SF,  "$TMP_DIR/$TMP_FILE";
  open OUT, ">$TMP_DIR/snp_feature_$SEQ_REGION_ID";
  LINE : while (<SF>) {
    my ($feature_id,$query_name,$query_start,$query_end,$null,$null1,$query_strand_old,$target_seq_region_id,$target_name,$target_start,$target_end,$target_strand,$score,$cigar_string,$individual_name) = split;

    for (my $i=0;$i<@sorted_start;$i++) {
      my $seq_region_start = $sorted_start[$i];
      my $variation_id = $rec_seq_start{$seq_region_start};
      my ($allele_string) = split /\_/, $rec_var{$variation_id};
      if ($seq_region_start > $target_start and $seq_region_start < $target_end) {
	print OUT "$feature_id\t$query_name\t$query_start\t$query_end\t$query_strand_old\t$target_name\t$target_start\t$target_end\t$target_strand\t$target_seq_region_id\t$variation_id\t$seq_region_start\t$allele_string\t$cigar_string\n";
       }
      elsif ($seq_region_start >$target_end) {
	next LINE;
      }
      elsif ($seq_region_start < $target_start) {
	my $shifted = shift @sorted_start;
	my $first = $sorted_start[0];
	$i=-1; #make sure always checking from the first start in @sorted_start array
      }
    }
  }
  system("mv $TMP_DIR/$TMP_FILE $TMP_DIR/ssahasnp_feature_$SEQ_REGION_ID");
  system ("mv $TMP_DIR/snp_feature_$SEQ_REGION_ID $TMP_DIR/$TMP_FILE");

  debug("Creating conmbined_feature table");
  create_and_load($dbSara,"combine_feature_$SEQ_REGION_ID","feature_id","query_name","query_start","query_end","query_strand_old","target_name","target_start","target_end","target_strand","target_seq_region_id i*","variation_id","vf_seq_region_start","allele_string","cigar_string");

  debug("Done for loading combine_feature_$SEQ_REGION_ID");

}

sub snp_pos {

  my $dbSara = shift;

  my ($feature_id,$query_name,$query_start,$query_end,$query_strand_old,$target_name,$target_start,$target_end,$target_strand,$target_seq_region_id,$variation_id,$seq_region_start,$allele_string,$cigar_string);

  debug("Prepare table snp_pos_$SEQ_REGION_ID");

  my $sth = $dbSara->prepare(qq{SELECT feature_id,query_name,query_start,query_end,query_strand_old,target_name,target_start,target_end,target_strand,target_seq_region_id,variation_id,vf_seq_region_start,allele_string,cigar_string
                             FROM combine_feature_$SEQ_REGION_ID
                             where variation_id=240663 and query_name="17000131958406"
			     #FROM combine_feature 
			     #WHERE target_seq_region_id = $SEQ_REGION_ID
   });
  $sth->execute();

  $sth->bind_columns(\$feature_id,\$query_name,\$query_start,\$query_end,\$query_strand_old,\$target_name,\$target_start,\$target_end,\$target_strand,\$target_seq_region_id,\$variation_id,\$seq_region_start,\$allele_string,\$cigar_string);

  open OUT, ">$TMP_DIR/$TMP_FILE" or die "can't open tmp_file $!";

  #in cigar line, query_start and query_end have been changed to forward strand if they are in reverse, but this will change in next release
  my $query_strand=1;

  while ($sth->fetch()){
    my ($start1,$start2);
    #print "cigar_string is $cigar_string\n";
    next if ! $cigar_string;
    my @blocks = ( $cigar_string =~ /(\d*[MDI])/g );
    if ($query_strand ==1) {
      $start1 = $query_start;
    }
    else {
      $start1 = $query_end;
    }
    if ($target_strand ==1) {
      $start2 = $target_start;
    }
    else {
      $start2 = $target_end;
    }

    my ($q_start,$q_end,$t_start,$t_end);

    foreach my $block ( @blocks ){
      my ($length) = ( $block =~ /^(\d*)/ );
      $length =1 if $length eq "";

      if ( $block =~ /M$/ ){
	if ($query_strand ==1) {
	  $q_start = $start1;
	  $q_end = $start1 + $length - 1;
	  $start1 = $q_end + 1;
	}
	else {
	  $q_end = $start1;
	  $q_start = $start1 - $length +1;
	  $start1 = $q_start -1;
	}

	if ($target_strand ==1) {
	  $t_start = $start2;
	  $t_end = $start2 + $length - 1 ;
	  $start2 = $t_end + 1;
	}
	else {
	  $t_end = $start2;
	  $t_start = $start2 - $length +1;
	  $start2 = $t_start -1;
	}
	#print STDERR "block: $block\tfeat1: $q_start-$q_end\tfeat2: $t_start-$t_end\n";
      }
      elsif ( $block =~ /I$/ ){
	if ($query_strand ==1) {
	  $start1 += $length;
	}
	else {
	  $start1 -= $length;
	}
	#print STDERR "block: $block\tfeat1: $q_start-$q_end\tfeat2: $t_start-$t_end START1 IS $start1\n";
      }
      elsif ( $block =~ /D$/ ){
	if ($query_strand ==1) {
	  $start2 += $length;
	}
	else {
	  $start2 -= $length;
	}
	#print STDERR "block: $block\tfeat1: $q_start-$q_end\tfeat2: $t_start-$t_end START2 IS $start2\n" ;
      }
      #calculate snp position in query sequence
      if ($seq_region_start >$t_start and $seq_region_start <$t_end) {
	my $snp_pos;
	my $distance = $seq_region_start-$t_start+1;
	if ($query_strand==$target_strand) {
	  $snp_pos = $q_start+$distance-1;
	  print OUT "$feature_id\t$query_name\t$query_strand_old\t$snp_pos\t$variation_id\t$target_seq_region_id\t$seq_region_start\t$allele_string\n";
	}
	else {
	  $snp_pos = $q_end-$distance+1;
	  print OUT "$feature_id\t$query_name\t$query_strand_old\t$snp_pos\t$variation_id\t$target_seq_region_id\t$seq_region_start\t$allele_string\n";
	}
	last;
      }
    }
  }
  create_and_load($dbSara,"snp_pos_$SEQ_REGION_ID","feature_id i*","query_name *","query_strand_old i","snp_pos i","variation_id i*","target_seq_region_id i*","seq_region_start i*","allele_string");

  debug("Done for loading snp_pos_$SEQ_REGION_ID");	
}

sub flanking_qual {

  my $dbSara = shift;

  #my $tableref = $dbSara->db_handle->selectall_arrayref(qq{show tables like "flanking_qual_$SEQ_REGION_ID"});
  #if ( $tableref->[0][0]) {
  #  debug("table flanking_qual_$SEQ_REGION_ID already exist");
  #  return;
  #}
  my $index_file = $index_filename;
  system("lsrcp $index_file\.pag /tmp/index_file\_$SEQ_REGION_ID\.pag");
  system("lsrcp $index_file\.dir /tmp/index_file\_$SEQ_REGION_ID\.dir");
  unlink "$TMP_DIR/flanking_qual_error" if (-e "$TMP_DIR/flanking_qual_error");

  open OUT, ">$TMP_DIR/flank\_$TMP_FILE";
  open ERR, ">>$TMP_DIR/flanking_qual_error";

  debug("Dumping snp_pos");
  dumpSQL($dbSara,qq{SELECT feature_id, variation_id, query_name, query_strand_old, snp_pos, allele_string
                                            FROM snp_pos_$SEQ_REGION_ID
                                            ORDER BY query_name 
  					 });

  open IN, "$TMP_DIR/$TMP_FILE";
  my ($feature_id,$variation_id,$query_name,$query_strand_old,$snp_pos,$allele_string);
  LINE : while (<IN>){
    my ($feature_id,$variation_id,$query_name,$query_strand_old,$snp_pos,$allele_string) = split;
    #my $t0 = [gettimeofday];
    my ($flank_5_qual,$flank_3_qual,$snp_base,$snp_qual) = get_flanking_seq($query_name,$query_strand_old,$snp_pos);
    #my $t1 = [gettimeofday];
    #print "Time to get flank from database ", tv_interval($t0,$t1),"\n";
    if (! $snp_qual or $snp_qual < 23 ) {
      print ERR "snp_qual < 23:$feature_id\t$variation_id\t$query_name\t$query_strand_old\t$snp_pos\t$allele_string\n";
      next LINE;
    }
    foreach my $qual (@$flank_5_qual,@$flank_3_qual) {
      #some of the snp near the end, flanking bases are less than 5
      if ($qual < 15) {
	print ERR "flank_qual < 15:$feature_id\t$variation_id\t$query_name\t$query_strand_old\t$snp_pos\t$allele_string\n";
        next LINE;
      }
    }
    print OUT "$feature_id\t$variation_id\t",join (" ",@$flank_5_qual,),"\t", join (" ",@$flank_3_qual),"\t$snp_base\t$snp_qual\t$query_name\t$snp_pos\n";
  }

  system("mv $TMP_DIR/flank\_$TMP_FILE $TMP_DIR/$TMP_FILE");
  create_and_load($dbSara,"flanking_new_qual_$SEQ_REGION_ID","feature_id i*","variation_id i*","flanking_5_qual","flanking_3_qual","snp_base","snp_qual","query_name *","snp_pos i*");
 
  debug("Done for loading flanking_qual_$SEQ_REGION_ID");
  unlink "/tmp/$index_file\_$SEQ_REGION_ID"; 
  return;
}

sub get_pfetch_sequence_and_quality {
    my $query_name = shift;
    my $index_filename = "/tmp/index_file_$SEQ_REGION_ID";
    my $fastq_index = Bio::Index::Fastq->new(-filename => $index_filename);

    if (! defined($fastq_index)) {
	$fastq_index = Bio::Index::Fastq->new(-filename => $index_filename);
    }

    if (exists($CACHE{$query_name})) {
        return $CACHE{$query_name};
    }
    
    my $seq = $fastq_index->fetch($query_name);
    die "could not fetch $query_name" unless ($seq);
    
    $CACHE{$query_name} = $seq;

    return $seq;
  }

sub get_flanking_seq {
  my ($query_name,$query_strand_old,$snp_pos) = @_;
  my ($qual_5,$qual_3,$snp_base,$snp_qual,$start,$end);
  #note here $snp_pos is in forward strand snp position
  #my $t0 = [gettimeofday];

  local($\) = undef;
  my $seqobj = &get_pfetch_sequence_and_quality($query_name);

  my $len = $seqobj->length;

  $start = $snp_pos-5;
  $end = $snp_pos+5;
  $start = 1 if $snp_pos -5 < 0;
  $end = $len if $snp_pos+5 > $len;

  if ($query_strand_old != -1) {
    $snp_base = $seqobj->baseat($snp_pos);
    $snp_qual = $seqobj->qualat($snp_pos);
    $qual_5 = $seqobj->subqual($start,$snp_pos-1);
    $qual_3 = $seqobj->subqual($snp_pos+1,$end);
  }
  else {
    my $seq = $seqobj->seq();
    reverse_comp(\$seq);
    $snp_base = substr($seq,$snp_pos-1,1);
    #$seqstr = $seqobj->baseat($snp_pos);
 
    #$seqstr =~ tr/ACGTacgt/TCGAtcga/;
    my $qual = $seqobj->qual();
    my @quals = reverse @{$seqobj->qual()};
    $snp_qual = $quals[$snp_pos-1];
    #$snp_qual = $seqobj->qualat($snp_pos);
    $qual_5 = [ @quals[$start-1..$snp_pos-2] ];
    $qual_3 = [ @quals[$snp_pos..$end-1] ];
  }


  #my $t1 = [gettimeofday];
  #print "Time to run pfetch: ",tv_interval($t0,$t1),"\n";

  return ($qual_5,$qual_3,$snp_base,$snp_qual);
  
}

sub gtype_allele {

  my $dbSara = shift;
  my $allele_file = "allele_$SEQ_REGION_ID";	 
  open GTY, ">$TMP_DIR/gtype_$TMP_FILE" or die "can't open tmp_file : $!";
  open ALE, ">$TMP_DIR/$allele_file" or die "can't open allele_file : $!";
  open ERR, ">>$TMP_DIR/gtype_allele_error";

  dumpSQL($dbSara,qq{SELECT f.variation_id, f.snp_base, s.allele_string, m.strain_name
                                          FROM flanking_qual_$SEQ_REGION_ID f, snp_pos_$SEQ_REGION_ID s, query_match_length_strain m
                                          WHERE s.variation_id=f.variation_id
                                          AND f.query_name = m.query_name
					  #GROUP BY f.variation_id,f.snp_base,s.allele_string,m.strain_name
                                         });

  system("uniq $TMP_DIR/$TMP_FILE > $TMP_DIR/$TMP_FILE.su");
  system("mv $TMP_DIR/$TMP_FILE.su $TMP_DIR/$TMP_FILE");

  open IN, "$TMP_DIR/$TMP_FILE" or die "Can't open tmp file\n";
  my (%rec_var_allele,%var_allele_string);
  
  while (<IN>) {
    my ($variation_id,$snp_base,$allele_string,$strain_name) = split;
    $snp_base = uc($snp_base);
    next if $rec_var_allele{$strain_name}{$variation_id}{$snp_base};
    $rec_var_allele{$strain_name}{$variation_id}{$snp_base}++;
    my ($allele_1,$allele_2) = split /\//, $allele_string;
    $allele_1 = uc($allele_1);
    $allele_2 = uc($allele_2);
    $var_allele_string{$strain_name}{$variation_id}{$allele_1}='ref_allele';
    $var_allele_string{$strain_name}{$variation_id}{$allele_2}='non_ref_allele';
    $var_allele_string{$strain_name}{$variation_id}{'ref_allele'} = $allele_1;
    $var_allele_string{$strain_name}{$variation_id}{'non_ref_allele'} = $allele_2;
  }

  foreach my $strain_name (keys %rec_var_allele) {
    foreach my $variation_id (keys %{$rec_var_allele{$strain_name}}) {
      my @snp_bases = keys %{$rec_var_allele{$strain_name}{$variation_id}};
      my $ref_allele = $var_allele_string{$strain_name}{$variation_id}{'ref_allele'};
      my $non_ref_allele = $var_allele_string{$strain_name}{$variation_id}{'non_ref_allele'};
      if (scalar @snp_bases >2) {
	print ERR "snp_bases >2 :$strain_name\t$variation_id\t@snp_bases\t$ref_allele/$non_ref_allele\n";
	next;
      }
      elsif (@snp_bases == 1) {
	if ($var_allele_string{$strain_name}{$variation_id}{$snp_bases[0]}) {
	  if ($snp_bases[0] eq "$ref_allele") {
	    print GTY "$variation_id\t$snp_bases[0]\t$snp_bases[0]\t$strain_name\n";
	    print ALE "$variation_id\t$snp_bases[0]\t$strain_name\n";
	  }
	  elsif ($snp_bases[0] eq "$non_ref_allele") {
	    print GTY "$variation_id\t$snp_bases[0]\t$snp_bases[0]\t$strain_name\n";
	    print ALE "$variation_id\t$snp_bases[0]\t$strain_name\n";
	    print ALE "$variation_id\t$ref_allele\trefstrain\n";
	  }
	}
	else {
	  print ERR "strain $strain_name with variation_id $variation_id don't have base $snp_bases[0] in $ref_allele/$non_ref_allele\n";
	}
      }
      elsif (@snp_bases == 2) {
	if ( $var_allele_string{$strain_name}{$variation_id}{$snp_bases[0]} and $var_allele_string{$strain_name}{$variation_id}{$snp_bases[1]} ) {
	  print GTY "$variation_id\t$snp_bases[0]\t$snp_bases[1]\t$strain_name\n";
	  print ALE "$variation_id\t$snp_bases[0]\t$strain_name\n";
	  print ALE "$variation_id\t$snp_bases[1]\t$strain_name\n";
	}
	else {
	  print ERR "strain $strain_name with variation_id $variation_id don't have both bases $snp_bases[0]/$snp_bases[1] in $ref_allele/$non_ref_allele\n";
	}
      }
    }
  }

  system("mv $TMP_DIR/gtype_$TMP_FILE $TMP_DIR/$TMP_FILE");
  create_and_load($dbSara,"tmp_individual_genotype_single_bp_$SEQ_REGION_ID","variation_id i*","allele_1","allele_2","sample_name");

  system("mv $TMP_DIR/$allele_file $TMP_DIR/$TMP_FILE");

  create_and_load($dbSara,"gtype_allele_$SEQ_REGION_ID","variation_id i*","allele","sample_name");

  debug("Done for gnotype and allele table");
}

sub read_coverage {

  open IN, "$alignment_file" or die "can't open alignment_file : $!";
  while (<IN>) {
    #ALIGNMENT 764 gnl|ti|900540082 10-1-110718848 5 794 39881618 39882407 F 790 99.75 795 SD
    my @all = split;
    if (@all == 13) {
      my ($null1,$score,$read_name,$target_name,$read_start,$read_end,$target_start,$target_end,$dir,$match_length,$identity,$read_length,$strain_name) = split;
      my ($null2,$null3,$chr) = split /\:/, $target_name;
      next if ($chr ne $SEQ_REGION_NAME);
      ($target_end, $target_start) = ($target_start, $target_end) if ($target_end < $target_start);
       my $my_score = $score/$read_length if $read_length != 0;
      if ($read_length == 0) {
	print "read_length=0 $_\n";
      }
      if ($my_score > 0.5) {
	open MAP, ">>$TMP_DIR/$chr\.mapped" or die "can't open mapping file:$!";
	print MAP "$strain_name\t$target_start\t$target_end\n";
      }
    }
  }
}

1;
