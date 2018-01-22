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
#use DBH;
use Time::HiRes qw(tv_interval gettimeofday);
use Bio::EnsEMBL::Utils::Cache;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);
use Bio::EnsEMBL::Utils::Sequence qw(expand reverse_comp);
use ImportUtils qw(debug load create_and_load dumpSQL);

our ($TMP_DIR, $TMP_FILE, $READS_FILE, $SEQ_REGION_ID, $SEQ_REGION_NAME, $LIMIT, $species, $source_id, $VAR_DBNAME, $JOB, $alignment_file,%CACHE,%REC_SEQ,%REC_QUAL);

{

  GetOptions('species=s'   => \$species,
	     'source_id=i'  => \$source_id,
	     'tmpdir=s'  => \$ImportUtils::TMP_DIR,
	     'tmpfile=s' => \$ImportUtils::TMP_FILE,
	     'reads_file=s' => \$READS_FILE,
             'seq_region_id=i'   => \$SEQ_REGION_ID,
	     'seq_region_name=s' => \$SEQ_REGION_NAME,
             'limit=s'  => \$LIMIT,
	     'job=s'     => \$JOB,
	     'alignment_file=s' =>\$alignment_file,
	     );

  my $registry_file ||= $Bin . "/ensembl.registry";

  Bio::EnsEMBL::Registry->load_all( $registry_file );

  my $cdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
  my $vdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
  my $vdb_sara = Bio::EnsEMBL::Registry->get_DBAdaptor($species,"sara\_$species");

  my $dbVar = $vdba->dbc;
  my $dbCore = $cdba;
  my $dbSara = $vdb_sara->dbc;

  #we need to use cache to speed up fetch read sequence and quality data
  my $MAX_READS_CACHE_SIZE = 500;
  tie(%CACHE, 'Bio::EnsEMBL::Utils::Cache', $MAX_READS_CACHE_SIZE);
  
  $VAR_DBNAME = $dbVar->dbname();
  $TMP_DIR  = $ImportUtils::TMP_DIR;
  $SEQ_REGION_ID = ($SEQ_REGION_ID) ? "_$SEQ_REGION_ID" : '';
  $ImportUtils::TMP_FILE .= "$SEQ_REGION_ID";
  $TMP_FILE = $ImportUtils::TMP_FILE;

  ssahasnp_feature($dbSara) if ($JOB eq "snp_pos");
  flanking_qual($dbSara) if ($JOB eq "flank");
  gtype_allele($dbSara) if ($JOB eq "gtype_allele");
  read_coverage() if ($JOB eq "read_coverage");
}

sub ssahasnp_feature {

  my ($dbSara) = @_;
  my ($SEQ_REGION_ID_NEW) = $SEQ_REGION_ID =~ /^\_(\d+)$/;
  debug("Make query to variation_feature table");
  my $sth = $dbSara->prepare(qq{SELECT v.variation_id,vf.seq_region_start,vf.allele_string 
                                FROM  $VAR_DBNAME.variation v,$VAR_DBNAME.variation_feature vf
                                WHERE v.source_id = $source_id
                                AND v.variation_id = vf.variation_id
                                AND vf.seq_region_id = $SEQ_REGION_ID_NEW
                                AND vf.map_weight=1
                                #ORDER BY seq_region_start #sorted later
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
  my ($SEQ_REGION_ID_NEW) = $SEQ_REGION_ID =~ /^\_(\d+)$/;

  #sort seq_start to speed up two big array comparison
  my @sorted_start = sort {$a<=>$b} keys %rec_seq_start;
  my @seq_region_ids;

  debug("Dumping ssahaSNP_feature");
  dumpSQL($dbSara, qq{SELECT * FROM ssahaSNP_feature WHERE target_seq_region_id = $SEQ_REGION_ID_NEW});
  `sort -k 8 -g -o $TMP_DIR/$TMP_FILE\_s $TMP_DIR/$TMP_FILE`;#order the file by the target_start
  system("mv $TMP_DIR/$TMP_FILE\_s $TMP_DIR/$TMP_FILE");
  open SF,  "$TMP_DIR/$TMP_FILE";
  open OUT, ">$TMP_DIR/snp_feature$SEQ_REGION_ID";
  LINE : while (<SF>) {
    my ($feature_id,$query_name,$query_start,$query_end,$query_strand,$target_seq_region_id,$target_name,$target_start,$target_end,$target_strand,$score,$cigar_string,$individual_name) = split;

    for (my $i=0;$i<@sorted_start;$i++) {
      my $seq_region_start = $sorted_start[$i];
      my $variation_id = $rec_seq_start{$seq_region_start};
      my $allele_string = $rec_var{$variation_id};
      if ($seq_region_start > $target_start and $seq_region_start < $target_end) {
	print OUT "$feature_id\t$query_name\t$query_start\t$query_end\t$query_strand\t$target_name\t$target_start\t$target_end\t$target_strand\t$target_seq_region_id\t$variation_id\t$seq_region_start\t$allele_string\t$cigar_string\n";
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
  system("mv $TMP_DIR/$TMP_FILE $TMP_DIR/ssahasnp_feature$SEQ_REGION_ID");
  system ("mv $TMP_DIR/snp_feature$SEQ_REGION_ID $TMP_DIR/$TMP_FILE");

  debug("Creating conmbined_feature table");
  create_and_load($dbSara,"combine_feature$SEQ_REGION_ID","feature_id","query_name","query_start","query_end","query_strand","target_name","target_start","target_end","target_strand","target_seq_region_id i*","variation_id i*","vf_seq_region_start","allele_string","cigar_string");

  debug("Done for loading combine_feature$SEQ_REGION_ID");
  unlink ("$TMP_DIR/ssahasnp_feature$SEQ_REGION_ID");
}

sub snp_pos {

  my $dbSara = shift;

  my ($feature_id,$query_name,$query_start,$query_end,$query_strand,$target_name,$target_start,$target_end,$target_strand,$target_seq_region_id,$variation_id,$seq_region_start,$allele_string,$cigar_string);

  debug("Prepare table snp_pos$SEQ_REGION_ID");

  my $tableref = $dbSara->db_handle->selectall_arrayref(qq{show tables like "snp_pos$SEQ_REGION_ID"});  
  if ( $tableref->[0][0]) {    
    debug("table snp_pos$SEQ_REGION_ID already exist");    
    $dbSara->do(qq{drop table snp_pos$SEQ_REGION_ID});  
  }

  dumpSQL($dbSara, qq{SELECT feature_id,query_name,query_start,query_end,query_strand,target_name,target_start,target_end,target_strand,target_seq_region_id,variation_id,vf_seq_region_start,allele_string,cigar_string
                             FROM combine_feature$SEQ_REGION_ID
                             #where variation_id=240663 and query_name="17000131958406"
   });

  open SF,  "$TMP_DIR/$TMP_FILE";
  open OUT, ">$TMP_DIR/snp_pos$SEQ_REGION_ID";

  while (<SF>) {
    my ($feature_id,$query_name,$query_start,$query_end,$query_strand,$target_name,$target_start,$target_end,$target_strand,$target_seq_region_id,$variation_id,$seq_region_start,$allele_string,$cigar_string) = split /\t/, $_;

    #in cigar line, query_start and query_end have been changed to forward strand if they are in reverse, but this will change in next release
    if ($query_strand==-1) {
      ($query_start,$query_end) = ($query_end,$query_start);
    }
    my ($start1,$start2);
    #print "cigar_string is $cigar_string\n";
    next if ! $cigar_string;
    my @blocks = ( $cigar_string =~ /(\d*[MDI])/g );
    if ($query_strand ==1) {
      $start1 = $query_start;
    } else {
      $start1 = $query_end;
    }
    if ($target_strand ==1) {
      $start2 = $target_start;
    } else {
      $start2 = $target_end;
    }

    my ($q_start,$q_end,$t_start,$t_end);

    foreach my $block ( @blocks ) {
      my ($length) = ( $block =~ /^(\d*)/ );
      $length =1 if $length eq "";

      if ( $block =~ /M$/ ) {
	if ($query_strand ==1) {
	  $q_start = $start1;
	  $q_end = $start1 + $length - 1;
	  $start1 = $q_end + 1;
	} else {
	  $q_end = $start1;
	  $q_start = $start1 - $length +1;
	  $start1 = $q_start -1;
	}

	if ($target_strand ==1) {
	  $t_start = $start2;
	  $t_end = $start2 + $length - 1 ;
	  $start2 = $t_end + 1;
	} else {
	  $t_end = $start2;
	  $t_start = $start2 - $length +1;
	  $start2 = $t_start -1;
	}
	#print STDERR "block: $block\tfeat1: $q_start-$q_end\tfeat2: $t_start-$t_end\n";
      } elsif ( $block =~ /I$/ ) {
	if ($query_strand ==1) {
	  $start1 += $length;
	} else {
	  $start1 -= $length;
	}
	#print STDERR "block: $block\tfeat1: $q_start-$q_end\tfeat2: $t_start-$t_end START1 IS $start1\n";
      } elsif ( $block =~ /D$/ ) {
	if ($target_strand ==1) {
	  $start2 += $length;
	} else {
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
	} else {
	  $snp_pos = $q_end-$distance+1;
	}
        print OUT "$feature_id\t$query_name\t$query_strand\t$snp_pos\t$variation_id\t$target_seq_region_id\t$seq_region_start\t$allele_string\n";
	last;
      }
    }
  }

  system("mv $TMP_DIR/snp_pos$SEQ_REGION_ID $TMP_DIR/$TMP_FILE");
  create_and_load($dbSara,"snp_pos$SEQ_REGION_ID","feature_id i*","query_name *","query_strand i","snp_pos i","variation_id i*","target_seq_region_id i*","seq_region_start i*","allele_string");

  debug("Done for loading snp_pos$SEQ_REGION_ID");
}

sub flanking_qual {

  my $dbSara = shift;

  my $tableref = $dbSara->db_handle->selectall_arrayref(qq{show tables like "flanking_qual$SEQ_REGION_ID"});
  if ( $tableref->[0][0]) {
    debug("table flanking_qual$SEQ_REGION_ID already exist");
    return;
  }

  open OUT, ">/tmp/flank\_$TMP_FILE";
  open ERR, ">$TMP_DIR/flanking_qual_error\_$TMP_FILE";
  my $LIMIT_SQL = ($LIMIT) ? " LIMIT $LIMIT " : '';

  debug("make_seq_qual_hash...");
  make_seq_qual_hash("$READS_FILE");
  debug("Dumping snp_pos");

  if (! -e "$TMP_DIR/$TMP_FILE") { 
    dumpSQL($dbSara,qq{SELECT feature_id, variation_id, query_name, query_strand, snp_pos, allele_string
                       FROM snp_pos$SEQ_REGION_ID
                       $LIMIT_SQL
                       #ORDER BY query_name 
     		    });
  }
  open IN, "$TMP_DIR/$TMP_FILE" or die "$TMP_DIR/$TMP_FILE not exist";
  LINE : while (<IN>){
    my ($feature_id,$variation_id,$query_name,$query_strand,$snp_pos,$allele_string) = split;
    #my $t0 = [gettimeofday];
    my $flank = get_flanking_seq($query_name,$query_strand,$snp_pos);
    if ($flank) {
      my ($flank_5_qual,$flank_3_qual,$snp_base,$snp_qual) = @$flank;
    
      #my $t1 = [gettimeofday];
      #print "Time to fetch sequence and get flanking qual  ", tv_interval($t0,$t1),"\n";
      if (! $snp_qual or $snp_qual < 23 ) {
        print ERR "snp_qual < 23\t$feature_id\t$variation_id\t$query_name\t$query_strand\t$snp_pos\t$allele_string\n";
        next LINE;
      }
      foreach my $qual (@$flank_5_qual,@$flank_3_qual) {
        #some of the snp near the end, flanking bases are less than 5
        if ($qual < 15) {
	  print ERR "flank_qual < 15\t$feature_id\t$variation_id\t$query_name\t$query_strand\t$snp_pos\t$allele_string\n";
          next LINE;
        }
      }  
      print OUT "$feature_id\t$variation_id\t",join (" ",@$flank_5_qual,),"\t", join (" ",@$flank_3_qual),"\t$snp_base\t$snp_qual\t$query_name\t$query_strand\t$snp_pos\n";
    }
    else {
      print "No return for seq_qual call for $query_name\n";
    }
  }
  my $table_name;
  if (!$LIMIT_SQL) {
    $table_name = "flanking_qual$SEQ_REGION_ID";
  }
  else {
    $table_name = "flanking_qual_$TMP_FILE";
  }
  system("mv /tmp/flank\_$TMP_FILE $TMP_DIR/$TMP_FILE");
  create_and_load($dbSara,"$table_name","feature_id i*","variation_id i*","flanking_5_qual","flanking_3_qual","snp_base","snp_qual","query_name *","query_strand","snp_pos i*");
  system("mv $TMP_DIR/flanking_qual_error\_$TMP_FILE $TMP_DIR/$TMP_FILE");
  create_and_load($dbSara,"failed\_$table_name","description","feature_id i*","variation_id i*","query_name","query_strand i","snp_pos i","allele_string");
  
  debug("Done for loading $table_name");
  #unlink "$index_name" if $TMP_DIR !~ /turing/; 
  return;
}

sub make_seq_qual_hash {
  my $reads_file = shift;
  my $name;
  #while (<$reads_file>) {
  #while (<$reads_file\*fastq /turing/mouse129_extra/yuan/watson/output_dir/missed1_dir/missed_name1\*fastq>) {
    #print "file is $_\n";
    my $file = $reads_file;
    open FASTQ, $file or die "Can't open $file";
    my $found;
    while (<FASTQ>) {
      chomp;
      if (/^\@(.*)$/) {
	$name = $1;
        ($name) = split /\s+/,$name;
      }
#      elsif (/^\+(.*)/) {
#	$name = $1;
#        ($name) = split /\s+/,$name;
#      }
      elsif ($name and /(^\!.*)$/) {
	$REC_QUAL{$name} = $1;
	undef $name;
      }
      elsif ($name and /[^?%*=<!+]/) {
	$REC_SEQ{$name} = $_;
	#undef $name;
      }
    }
    close FASTQ;
  #}
}

sub get_pfetch_sequence_and_quality {
    my $index_filename = shift;
    my $query_name = shift;
    
    #my $index_filename = "/tmp/index_file_sara";
    #my $index_filename = "/turing/mouse129_extra/yuan/human_celera/fastq/index_file";
    #my $index_filename = "/gvar/hum-snp5/yuan/fastq/index_file";
    #my $index_filename = "/ecs2/scratch6/yuan/rat/CELERA/reads/index";
    my $fastq_index = Bio::Index::Fastq->new(-filename => $index_filename);

    if (exists($CACHE{$query_name})) {
        return $CACHE{$query_name};
    }
    
    my $seq = $fastq_index->fetch($query_name);
    die "could not fetch $query_name" unless ($seq);
    
    $CACHE{$query_name} = $seq;

    return $seq;
  }

sub get_flanking_seq {
  my ($query_name,$query_strand,$snp_pos) = @_;
  my ($qual_5,$qual_3,$snp_base,$snp_qual,$start,$end);
  #my $t0 = [gettimeofday];

  die "hash %REC_SEQ is empty check you reads_file" unless (keys %REC_SEQ >1);

  local($\) = undef;
  #my $seqobj = &get_pfetch_sequence_and_quality($index_name,$query_name);
  my $seq = $REC_SEQ{"$query_name"};
  my $qual = $REC_QUAL{"$query_name"};
  if (!$seq or !$qual) {
    print "query $query_name has a snp_pos $snp_pos, don't have seq or qual\n";
    return;
  }
  #my $t1 = [gettimeofday];
  #print "Time to fetch seq only : ",tv_interval($t0,$t1),"\n";
  my $len = length($seq);

  $start = $snp_pos-5;
  $end = $snp_pos+5;
  $start = 1 if $snp_pos -5 <= 0;
  $end = $len if $snp_pos+5 > $len;
  #print "seq is $seq and snp_pos is $snp_pos\n";
  $snp_base = substr($seq,$snp_pos-1,1);
  $snp_qual = substr($qual,$snp_pos-1,1);
  $qual_5 = substr($qual,$start-1,5);
  $qual_3 = substr($qual,$snp_pos,5);
  if ($query_strand ==-1) {
    $snp_base =~ tr/ACGTacgt/TGCAtgca/;
    ($qual_5,$qual_3) = ($qual_3,$qual_5);
  }
              
  my @qual_5 = split '',$qual_5;
  my @qual_3 = split '',$qual_3;
  
  if ($query_strand ==-1) {
    @qual_5 = reverse @qual_5;
    @qual_3 = reverse @qual_3;
  }
  
  @qual_5 = map{unpack("C",$_)-33} @qual_5;
  @qual_3 = map{unpack("C",$_)-33} @qual_3;
  $snp_qual = unpack("C",$snp_qual)-33;

  $qual_5 = [ @qual_5 ];
  $qual_3 = [ @qual_3 ];

  #my $t2 = [gettimeofday];
  #print "Time to fetch and calculate flanking qual: ",tv_interval($t0,$t2),"\n";

  return ([$qual_5,$qual_3,$snp_base,$snp_qual]);
  #exit;
}

sub gtype_allele {

  my $dbSara = shift;
  my $LIMIT_SQL = ($LIMIT) ? " LIMIT $LIMIT " : '';
  
  my $allele_file = "allele_$TMP_FILE";	 
  open GTY, ">$TMP_DIR/gtype_$TMP_FILE" or die "can't open tmp_file : $!";
  open ALE, ">$TMP_DIR/$allele_file" or die "can't open allele_file : $!";
  open ERR, ">$TMP_DIR/gtype_allele_error_$TMP_FILE";

  my $aref = $dbSara->db_handle->selectall_arrayref(qq{SHOW TABLES LIKE "flanking_qual$SEQ_REGION_ID"});
  if (!$aref->[0][0]) {
    return ;
  }
  my $tableref = $dbSara->db_handle->selectall_arrayref(qq{show tables like "tmp_individual_genotype_single_bp$SEQ_REGION_ID"});
  if ( $tableref->[0][0]) {
    debug("table tmp_individual_genotype_single_bp$SEQ_REGION_ID already exist");
    return;
  }
                
  dumpSQL($dbSara,qq{SELECT f.variation_id, f.snp_base, s.allele_string, m.individual_name
                                          FROM flanking_qual$SEQ_REGION_ID f, snp_pos$SEQ_REGION_ID s, ssahaSNP_feature m
                                          WHERE s.variation_id=f.variation_id
                                          AND f.feature_id = s.feature_id
                                          AND f.query_name = s.query_name
                                          AND f.query_name = m.query_name
					  $LIMIT_SQL
                                          #AND f.variation_id=922
                                          #GROUP BY f.variation_id,f.snp_base,s.allele_string,m.strain_name
                                         });

  system("uniq $TMP_DIR/$TMP_FILE > $TMP_DIR/$TMP_FILE.su");
  system("mv $TMP_DIR/$TMP_FILE.su $TMP_DIR/$TMP_FILE");

  open IN, "$TMP_DIR/$TMP_FILE" or die "Can't open tmp file\n";
  my (%rec_var_allele,%var_allele_string,%ref_allele_done);
  
  while (<IN>) {
    my ($variation_id,$snp_base,$allele_string,$strain_name) = split;
    $snp_base = uc($snp_base);
    #$rec_var_allele{$strain_name}{$variation_id}++;
    $rec_var_allele{$strain_name}{$variation_id}{$snp_base}++;
    next if $rec_var_allele{$strain_name}{$variation_id}{$snp_base}>1;
    my ($allele_1,$allele_2,$allele_3,$allele_4) = split /\//, $allele_string;
    $allele_1 = uc($allele_1);
    $allele_2 = uc($allele_2);
    $allele_3 = uc($allele_3);
    $allele_4 = uc($allele_4);
    $var_allele_string{$strain_name}{$variation_id}{$allele_1}='ref_allele';
    $var_allele_string{$strain_name}{$variation_id}{$allele_2}='non_ref_allele1';
    $var_allele_string{$strain_name}{$variation_id}{$allele_3}='non_ref_allele2';
    $var_allele_string{$strain_name}{$variation_id}{$allele_4}='non_ref_allele3';
    $var_allele_string{$strain_name}{$variation_id}{'ref_allele'} = $allele_1;
    $var_allele_string{$strain_name}{$variation_id}{'non_ref_allele1'} = $allele_2;
    $var_allele_string{$strain_name}{$variation_id}{'non_ref_allele2'} = $allele_3;
    $var_allele_string{$strain_name}{$variation_id}{'non_ref_allele3'} = $allele_4;
  }

  foreach my $strain_name (keys %rec_var_allele) {
    foreach my $variation_id (keys %{$rec_var_allele{$strain_name}}) {
      my @snp_bases = keys %{$rec_var_allele{$strain_name}{$variation_id}};
      my ($ref_allele,$non_ref_allele1,$non_ref_allele2,$non_ref_allele3);
      $ref_allele = $var_allele_string{$strain_name}{$variation_id}{'ref_allele'};
      $non_ref_allele1 = $var_allele_string{$strain_name}{$variation_id}{'non_ref_allele1'};
      $non_ref_allele2 = $var_allele_string{$strain_name}{$variation_id}{'non_ref_allele2'};
      $non_ref_allele3 = $var_allele_string{$strain_name}{$variation_id}{'non_ref_allele3'};
      if (scalar @snp_bases >2) {
	my @new_snp_bases ;
	#if more than 2 snp_bases, check any snp_base only have one read support, get rid of it
	foreach my $base (@snp_bases) {
	  push @new_snp_bases, $base if $rec_var_allele{$strain_name}{$variation_id}{$base}>1;
	}
	if (scalar @new_snp_bases >2 ) {
	  print ERR "snp_bases >2\t$strain_name\t$variation_id\t@new_snp_bases\t$ref_allele/$non_ref_allele1,$non_ref_allele2,$non_ref_allele3\n";
	  next;
	}
	elsif (@new_snp_bases == 1) {
	  if ($var_allele_string{$strain_name}{$variation_id}{$new_snp_bases[0]}) {
	    print GTY "$variation_id\t$new_snp_bases[0]\t$new_snp_bases[0]\t$strain_name\n";
	    print ALE "$variation_id\t$new_snp_bases[0]\t$strain_name\n";
	    print ALE "$variation_id\t$ref_allele\trefstrain\n" if !$ref_allele_done{$variation_id};
	    $ref_allele_done{$variation_id}=1;
	  }
	  else {
	    print ERR "snp_base not in allele_string\t$strain_name\t$variation_id\t$new_snp_bases[0]\t$ref_allele/$non_ref_allele1,$non_ref_allele2,$non_ref_allele3\n";
	  }
	}
	elsif (@new_snp_bases == 2 and $species !~ /mus|mouse/i) {
	  if ( $var_allele_string{$strain_name}{$variation_id}{$new_snp_bases[0]} and $var_allele_string{$strain_name}{$variation_id}{$new_snp_bases[1]} ) {
	    print GTY "$variation_id\t$new_snp_bases[0]\t$new_snp_bases[1]\t$strain_name\n";
	    print ALE "$variation_id\t$new_snp_bases[0]\t$strain_name\n";
	    print ALE "$variation_id\t$new_snp_bases[1]\t$strain_name\n";
	  }
	  else {
	    print ERR "snp_base not in allele_string\t$strain_name\t$variation_id\t@new_snp_bases\t$ref_allele/$non_ref_allele1,$non_ref_allele2,$non_ref_allele3\n";
	  }
	}
      }

      elsif (@snp_bases == 1) {#problem of snp allele failed, so not a snp anymore????see variation_id=943938 for tetraodon
        #if (scalar keys %rec_var_allele ==1 and $snp_bases[0] ne $ref_allele or scalar keys %rec_var_allele > 1) {#if only one snp_base survived, it has to be different from ref_allele, we give homozygoes genotype, but this only for run with single individual, for multi individual, allow single base same as reference base
        #if ($snp_bases[0] eq $ref_allele) {#for this case (run against dbSNP SNPs), we expecting all bases are same as reference base also for multi individuals, this individual may have alleles all same as reference one
	if ($var_allele_string{$strain_name}{$variation_id}{$snp_bases[0]}) {##for multi individual
	  print GTY "$variation_id\t$snp_bases[0]\t$snp_bases[0]\t$strain_name\n";
          print ALE "$variation_id\t$snp_bases[0]\t$strain_name\n";
          print ALE "$variation_id\t$ref_allele\trefstrain\n" if !$ref_allele_done{$variation_id};
	  $ref_allele_done{$variation_id}=1;
        }
        else {
	  print ERR "snp_base not in allele_string\t$strain_name\t$variation_id\t$snp_bases[0]\t$ref_allele/$non_ref_allele1,$non_ref_allele2,$non_ref_allele3\n";
	}
      }

      elsif (@snp_bases == 2 and $species !~ /mus|mouse/i) {
	if ( $var_allele_string{$strain_name}{$variation_id}{$snp_bases[0]} and $var_allele_string{$strain_name}{$variation_id}{$snp_bases[1]} ) {
	  print GTY "$variation_id\t$snp_bases[0]\t$snp_bases[1]\t$strain_name\n";
	  print ALE "$variation_id\t$snp_bases[0]\t$strain_name\n";
	  print ALE "$variation_id\t$snp_bases[1]\t$strain_name\n";
	}
	else {
	  print ERR "snp_base not in allele_string\t$strain_name\t$variation_id\t@snp_bases\t$ref_allele/$non_ref_allele1,$non_ref_allele2,$non_ref_allele3\n";
	}
      }
      else {
	print ERR "snp_bases = 2\t$strain_name\t$variation_id\t@snp_bases\t$ref_allele/$non_ref_allele1,$non_ref_allele2,$non_ref_allele3\n";
      }
    }
  }

  my ($gtype_table_name, $allele_table_name,$failed_table_name);
  if (! $LIMIT_SQL) {
    $gtype_table_name = "tmp_individual_genotype_single_bp$SEQ_REGION_ID";
    $allele_table_name = "gtype_allele$SEQ_REGION_ID";
    $failed_table_name = "failed_gtype$SEQ_REGION_ID";
  }
  else {
    $gtype_table_name = "tmp_individual_genotype_single_bp_$TMP_FILE";
    $allele_table_name = "gtype_allele_$TMP_FILE";
    $failed_table_name = "failed_gtype_$TMP_FILE";
  }

  system("mv $TMP_DIR/gtype_$TMP_FILE $TMP_DIR/$TMP_FILE");
  create_and_load($dbSara,"$gtype_table_name","variation_id i*","allele_1","allele_2","sample_name");

  system("mv $TMP_DIR/$allele_file $TMP_DIR/$TMP_FILE");

  create_and_load($dbSara,"$allele_table_name","variation_id i*","allele","sample_name");

  system("mv $TMP_DIR/gtype_allele_error_$TMP_FILE $TMP_DIR/$TMP_FILE");

  create_and_load($dbSara,"$failed_table_name","description","strain_name","variation_id i*","gtype","allele_string");

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
