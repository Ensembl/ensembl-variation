#! /usr/local/bin/perl

use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::FeaturePair;
use ImportUtils qw(dumpSQL debug create_and_load load );
use FindBin qw( $Bin );
use Bio::EnsEMBL::Variation::Utils::Sequence qw(unambiguity_code);
use Bio::EnsEMBL::Utils::Sequence qw(expand reverse_comp);
use Data::Dumper;

my ($species,$input_file_name, $chr, $output_dir,$target_dir, $input_dir);

GetOptions('species=s'       => \$species,
	   'input_file_name=s'    => \$input_file_name,
	   'chr=s'           => \$chr,
	   'input_dir=s'      => \$input_dir,
           'output_dir=s'      => \$output_dir,
	   'target_dir=s'      => \$target_dir);

my $registry_file;
$registry_file ||= $Bin . "/ensembl.registry";
Bio::EnsEMBL::Registry->load_all( $registry_file );

$species ||='human';

my $cdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
#my $vdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
my $dbCore = $cdb;
#my $dbVar = $vdb->dbc;

my $sa = $dbCore->get_SliceAdaptor();
#$input_dir ||= "/turing/mouse129_extra/yuan/LRG/input_dir";
#$target_dir ||= "/turing/mouse129_extra/yuan/human";
#$output_dir ||= "/turing/mouse129_extra/yuan/LRG/output_dir";

$input_dir ||= "/lustre/work1/ensembl/yuan/SARA/LRG/input_dir";
$output_dir ||= "/lustre/work1/ensembl/yuan/SARA/LRG/output_dir";
$target_dir ||= "/lustre/work1/ensembl/yuan/SARA/human/ref_seq_hash";

my $queue = "-q normal -R'select[mem>10000] rusage[mem=10000]'";

my $output_file_name = "ssaha2_output\_$input_file_name";
my $input_file = "$input_dir/$input_file_name";
my $output_file ="$output_dir/$output_file_name";


my ($subject, %rec_seq, %rec_find, %input_length, %done);

if ($chr) {
  $subject = "$target_dir/$chr\.fa";
}
else {
  $subject = "$target_dir/ref";
}

my $name;
open INPUT, "$input_file" or die "can't open $input_file ";
while (<INPUT>) {#to get input name/sequence hash;
  if (/^\>/) {
    s/^\>|\n//g;
    $name = $_;
  }
  else {
    s/\n//;
    $rec_seq{$name} .= $_;
  }
}

foreach my $name (keys %rec_seq) {#get query seq object
  #print "name is $name and seq is $rec_seq{$name}\n";
  my $seqobj = Bio::PrimarySeq->new(-id => $name, -seq => $rec_seq{$name});
  $rec_seq{$name} = $seqobj;
}

#bsub_ssaha_job($queue,$input_file,$subject);

my $call = "bsub -P ensembl-variation -K -w 'done($input_file\_ssaha_job)' -J waiting_process sleep 1"; #waits until all ssaha jobs have finished to continue
#system($call);

parse_ssaha2_out ();

sub bsub_ssaha_job {
  
  my ($queue, $input_file, $subject) = @_ ;

  #my $ssaha_command = "/nfs/acari/yuan/ensembl/src/ensembl-variation/scripts/ssahaSNP/ssaha2/ssaha2_v1.0.9_ia64/ssaha2";
  my $ssaha_command = "/nfs/acari/yuan/ensembl/src/ensembl-variation/scripts/ssahaSNP/ssaha2/ssaha2_v1.0.9_x86_64/ssaha2";
  #$ssaha_command .= " -align 0 -kmer 2 -seeds 2 -cut 1000 -output vulgar -depth 10 -best 1 $input_dir/test.fa $input_file";
  $ssaha_command .= " -align 1 -kmer 12 -seeds 4 -cut 1000 -output vulgar -depth 10 -best 1 -save $subject $input_file";
  my $call = "bsub -J $input_file\_ssaha_job -P ensembl-variation $queue -e $output_dir/error_ssaha -o $output_file $ssaha_command";
  system ($call);
}


sub parse_ssaha2_out {

  open OUT, ">$output_dir/mapping_file_$input_file_name" or die "can't open output_file : $!\n";

  open SSAHA, "$output_file" or die "can't open $output_file : $!\n";

  while (<SSAHA>) {
    #print ;
    next if ($_ !~ /^vulgar/);
    if (/^vulgar\:\s+(\S+)\s+(\d+)\s+(\d+)\s+(\+|\-)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\+|\-)\s+(\d+)\s+(.*)$/) {
      #next unless $1 eq "rs16822290"
      #print;
      my $q_strand = $4;
      my $t_strand = $8;

      my $h={};
      $h->{'q_id'} = $1;
      $h->{'q_start'} = $2;
      $h->{'q_end'} = $3;
      $h->{'t_id'} = $5;
      $h->{'t_start'} = $6;
      $h->{'t_end'} = $7;
      $h->{'t_strand'} = $t_strand;
      $h->{'score'} = $9;
      $h->{'match'} = $10;
      $q_strand = ($q_strand =~ /\-/) ? -1 : 1;
      $t_strand = ($t_strand =~ /\-/) ? -1 : 1;
      $h->{'q_strand'} = $q_strand;
      $h->{'t_strand'} = $t_strand;
      push @{$rec_find{$h->{'q_id'}}}, $h if $h->{'q_id'};
    }
  }

  my ($pairs,$feature_pairs);
  foreach my $q_id (keys %rec_find) {
    my @h = sort {$b->{'score'}<=>$a->{'score'}} @{$rec_find{$q_id}};
    if ($h[0]->{'score'} > $h[1]->{'score'} and @h>=2) {
      ($pairs,$feature_pairs) = find_results($h[0]);
    }
    elsif ($h[1]->{'score'} > $h[2]->{'score'} and @h>=3) {
      ($pairs,$feature_pairs) = find_results($h[0],$h[1]);
    }
    elsif ($h[2]->{'score'} > $h[3]->{'score'} and @h>=4) {
      ($pairs,$feature_pairs) = find_results($h[0],$h[1],$h[2]);
    }
    else {
      print "\tmore than 3 hits having same score for $q_id\n";
    }
  }

  get_annotations($feature_pairs);

  my ($total_seq,$no);

  $total_seq = keys %rec_seq;

  foreach my $q_id (keys %rec_seq) {
    if (!$done{$q_id}) {
      $no++;
    }
  }

  print "$no out of $total_seq are not mapped\n";
}

sub find_results {

  my @pairs;
  my ($h1,$h2,$h3) = @_;
  my @fps;

  foreach my $h ($h1,$h2,$h3) {
    next if ! $h->{'q_id'};
    my $q_id = $h->{'q_id'};
    my $q_start = $h->{'q_start'};
    my $q_end = $h->{'q_end'};
    my $q_strand = $h->{'q_strand'};
    my $t_id = $h->{'t_id'};
    my $t_start = $h->{'t_start'};
    my $t_end = $h->{'t_end'};
    my $t_strand = $h->{'t_strand'};
    my $score = $h->{'score'};
    my $match = $h->{'match'};

    my ($f_q_start,$f_q_end) ;#feature start-end used to make feature_pairs
    if ($q_strand ==1) {
      ($f_q_start,$f_q_end) = ($q_start,$q_end);
    }
    else {
      ($f_q_start,$f_q_end) = ($q_end,$q_start);
    }

    my ($seq_region_name) = split /\-/, $t_id;
    print "seq_region_name is $seq_region_name and t_start is $t_start, t_end is $t_end\n";
    my $slice = $sa->fetch_by_region('chromosome',$seq_region_name,$t_start,$t_end, 1);
    #print "slice is ",ref($slice),"\n";
    print "slice seq_region_name is ",$slice->seq_region_name,'-',length($slice->seq),'-',$slice->length,'-',$slice->start,'-',$slice->end,"\n";
    #print "slice seq is ",$slice->seq,"\n";
    my $q_seqobj = $rec_seq{$q_id};

    #warning that query_seq length != matched query_length
    if (length($q_seqobj->seq) != abs($q_end - $q_start) + 1){
      print "q_start is $q_start and q_end is $q_end and length is ",length($q_seqobj->seq),"\n";
      die("query sequence length not equal matched query sequence length");
    }

    my @match_components = split /\s+/, $match;

    my ($no_gap,$mis_match,$full_match,$new_q_start, $new_q_end, $new_t_start, $new_t_end);

    print "$q_id,$q_start,$q_end,$q_strand,$t_id,$t_start,$t_end,$t_strand,$score,@match_components\n";
    if (scalar @match_components ==3 and $match_components[0] eq 'M') {
      $no_gap=1;
    }

    while (@match_components) {

      my $type = shift @match_components;
      my $query_match_length = shift @match_components;
      my $target_match_length = shift @match_components;
      my ($tmp_q_start,$tmp_q_end);
      
      if ($type eq 'M') {#go through each base to check SNP,target strand always 1
	my ($q_seq);
	if ($q_strand == 1) {
	  $new_q_start = $q_start + $query_match_length - 1 ;
	  $q_seq = substr($q_seqobj->seq,$q_start-1,$query_match_length);
	}
	else {
	  $new_q_start = $q_start - $query_match_length + 1;
	  $q_seq = substr($q_seqobj->seq,$new_q_start-1,$query_match_length);
	  reverse_comp(\$q_seq);
	}
	$new_t_start = $t_start + $target_match_length -1;
	print "$q_start,$new_q_start,$t_start,$new_t_start\n";
	my $t_seq = $slice->seq;
	print "$q_start,$new_q_start,$t_start,$new_t_start and length of q_seq ",length($q_seq),'-',length($t_seq),"\n";
	#print "q_seq_is $q_seq\n";
	#print "t_seq is $t_seq\n";
	my $q_count = 1;
	my $t_count = 1;

	my %q_seqs = map {$q_count++,$_} split '', $q_seq;
	my %t_seqs = map {$t_count++,$_} split '', $t_seq;
	my ($sub_q_end);
	my $sub_q_start = $q_start;#initial sub_q_start
	my $sub_t_start = $t_start;
	foreach my $count (sort {$a<=>$b} keys %q_seqs) {
	  if ($q_seqs{$count} !~ /$t_seqs{$count}/i) {
	    #next;
	    print "count is $count\n";
	    $mis_match=1;
	    my $sub_t_end = $t_start + ($count -1) - 1;
	    if ($q_strand ==1) {
	      $sub_q_end = $q_start + ($count -1) - 1;
	    }
	    else {
	      $sub_q_end = $q_start - ($count -1) + 1;
	    }
	    if ($q_strand==1) {
	      $tmp_q_start = $sub_q_start;
	      $tmp_q_end = $sub_q_end;
	      $sub_q_start = $sub_q_end+2;#+1 cover SNP base, +1 cover next base
	    }
	    else {
	      $tmp_q_start = $sub_q_end;
	      $tmp_q_end = $sub_q_start;
	      $sub_q_start = $sub_q_end-2;
	    }
	
	    ($tmp_q_start,$tmp_q_end) = ($tmp_q_end,$tmp_q_start) if($tmp_q_end<$tmp_q_start);
	    push @pairs, [$type,$tmp_q_start,$tmp_q_end,$sub_t_start,$sub_t_end,$q_strand];

	    $sub_t_start = $sub_t_end+2;
	  }
	}

	$q_start = $sub_q_start;
	$t_start = $sub_t_start;
	
	print "q_start is $q_start and t_start is $t_start\n";
        print "in M, q_start is $q_start and q_end is $new_q_start and t_start is $t_start and t_end is $new_t_start\n";

	if  ($q_strand==1) {
	  $tmp_q_start = $q_start;
	  $tmp_q_end = $new_q_start;
	}
	else {
	  $tmp_q_start = $new_q_start;
	  $tmp_q_end = $q_start;
	}
	($tmp_q_start,$tmp_q_end) = ($tmp_q_end,$tmp_q_start) if($tmp_q_end<$tmp_q_start);
	push @pairs, [$type,$tmp_q_start,$tmp_q_end,$t_start,$new_t_start,$q_strand];
      }
      elsif ($type eq 'G') {
	if ($q_strand ==1) {
	  $new_q_start = $q_start + $query_match_length + 1;
	}
	else {
	  $new_q_start = $q_start - $query_match_length - 1;
	}
	$new_t_start = $t_start + $target_match_length +1;
	if ($q_strand ==-1) {
	  $tmp_q_start = $new_q_start+1;
	  $tmp_q_end = $q_start-1;
	}
	else {
	  $tmp_q_start = $q_start+1;
	  $tmp_q_end = $new_q_start-1;
	}
	($tmp_q_start,$tmp_q_end) = ($tmp_q_end,$tmp_q_start) if($tmp_q_end<$tmp_q_start);
	push @pairs, [$type,$tmp_q_start,$tmp_q_end,$t_start+1,$new_t_start-1,$q_strand];

	print "in G, q_start is $q_start and q_end is $new_q_start and t_start is $t_start and t_end is $new_t_start\n";
      }

      $q_start = $new_q_start;
      $t_start = $new_t_start;
    }

    $done{$q_id} =1 if $pairs[0]->[0];


    #make featurePair here, 
    if ($no_gap and ! $mis_match) {
      $full_match=1;
    }

    my $fp = Bio::EnsEMBL::FeaturePair->new(-start    => $f_q_start,
					    -end      => $f_q_end,
					    -strand   => $q_strand,
					    -display_id => $h->{'q_id'},
					    -hstart   => $h->{'t_start'},
					    -hend     => $h->{'t_end'},
					    -hstrand  => $h->{'t_strand'},
					    -hseqname => $h->{'t_id'},
					    -slice    => $slice,
					    -type     => $full_match,
					   );
    $fp->seqname($h->{'q_id'});
    push @fps, $fp;

    foreach my $pair (sort {$a->[3]<=>$b->[3]} @pairs) {
      print "pairs are ",$pair->[0],'-',$pair->[1],'-',$pair->[2],'-',$pair->[3],'-',$pair->[4],"\n";
    }
  }
  return \@pairs,\@fps;
}

sub seqname {

  my $self = shift;
  $self->{'seqname'} = shift if(@_);
  return $self->{seqname};

}
sub get_annotations {
  my $feature_pairs = shift;
  my $lrg_name = 'LRG2';

  foreach my $fp (@$feature_pairs) {
    my $full_match = $fp->type;print "full_match is $full_match\n";
    my ($q_start,$q_end,$t_start,$t_end,$q_strand);
    my $slice = $fp->slice; print "slice is ",ref($slice),"\n";

    $q_start = $fp->start;
    $q_end = $fp->end;
    $t_start = $fp->hstart;
    $t_end = $fp->hend;
    $q_strand = $fp->strand;

    if (!$full_match) {#it's not full match, insert data in several tables
      my $csa = $dbCore->get_CoordSystemAdaptor();
      my $cs = $csa->fetch_all_by_name('LRG');
      my $cs_id = $cs->[0]->dbID();
      my $q_seq_region_id;
      my $q_seqobj = $rec_seq{$fp->seqname};
      my $q_seq = $q_seqobj->seq;
      my $q_seq_length = length($q_seq);
      my $t_seq_region_id = $sa->get_seq_region_id($slice);
      my $lrg_name_ref =$dbCore->dbc->db_handle->selectall_arrayref(qq{SELECT name from seq_region WHERE name="$lrg_name"});
      if (! $lrg_name_ref->[0][0]) {
	$dbCore->dbc->do(qq{INSERT INTO seq_region(name,coord_system_id,length)values("$lrg_name",$cs_id,$q_seq_length)});
	$q_seq_region_id = $dbCore->dbc->db_handle->{'mysql_insertid'};
	$dbCore->dbc->do(qq{INSERT INTO dna(seq_region_id,sequence)values($q_seq_region_id,"$q_seq")});
      }
      else {
	my $q_seqid_ref = $dbCore->dbc->db_handle->selectall_arrayref(qq{SELECT seq_region_id from seq_region WHERE name="$lrg_name"});
	$q_seq_region_id = $q_seqid_ref->[0][0];
      }

      if ($q_end-$q_start == $t_end-$t_start) {
	#$dbCore->dbc->do(qq{INSERT IGNORE INTO assembly(asm_seq_region_id,cmp_seq_region_id,asm_start,asm_end,cmp_start,cmp_end,ori)values($t_seq_region_id,$q_seq_region_id,$t_start,$t_end,$q_start,$q_end,$q_strand)});
      }
      else {
	throw("distance between query and target is not the same, there is a indel");
      }
    }

      my $hslice = $sa->fetch_by_region('LRG',"$lrg_name");
      print "q_seq from hslice is ",$hslice->start,'-'.$hslice->end,"\n";
      #my $exp_slice      = $hslice->expand( 10000, 10000);
      #print "the expanded slice length is ",length($exp_slice->seq),"\n";
      my $contig_slice = $sa->fetch_by_region('contig',"AC015909.14.1.217746");
      my $chr_slice = $hslice->project('chromosome');
      #print Dumper($chr_slice);
      my @genes = @{$$chr_slice[0]->to_Slice->get_all_Genes()};
      my @transcripts =  @{$$chr_slice[0]->to_Slice->get_all_Transcripts()};
      print "gene name is ", $genes[0]->stable_id,'-',$genes[0]->start,'-',$genes[0]->end,"\n" if defined $genes[0];
      print "trans name is ", $transcripts[0]->stable_id,'-',$transcripts[0]->start,'-',$transcripts[0]->end,"\n" if defined $transcripts[0];
      my @gene_dbentries = @{$genes[0]->get_all_DBEntries()};
      foreach my $dbe (@gene_dbentries) {
	my @gene_syns = @{$dbe->get_all_synonyms()};
	print "all_synonyms are @gene_syns\n";
	print $dbe->db_display_name(),'-',$dbe->description(),'-',$genes[0]->external_db(),"\n";
      }
      my @exons = $transcripts[0]->get_all_Exons();#having problem
      print "exons name again is ", $exons[0]->stable_id,'-',$exons[0]->start,'-',$exons[0]->end,"\n" if defined $exons[0];
      my @transcript = $genes[0]->get_all_Transcripts();#having problem
      print "trans name again is ", $transcript[0]->stable_id,'-',$transcript[0]->start,'-',$transcript[0]->end,"\n" if defined $transcript[0];
      #print Dumper($chr_slice);
      #print Dumper($contig_slice);
  }

}



