#! /usr/local/bin/perl

use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::EnsEMBL::MappedSliceContainer;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::FeaturePair;
use ImportUtils qw(dumpSQL debug create_and_load load );
use FindBin qw( $Bin );
use Bio::EnsEMBL::Variation::Utils::Sequence qw(unambiguity_code);
use Bio::EnsEMBL::Utils::Sequence qw(expand reverse_comp);
use Data::Dumper;

my ($species,$input_file_name, $chr, $output_dir,$target_dir, $input_dir);

our $SSAHA_BIN = 'ssaha2';

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
my $asma = Bio::EnsEMBL::Registry->get_adaptor($species,"core","assemblymapper");
my $csa = Bio::EnsEMBL::Registry->get_adaptor($species,"core","coordsystem");





my $queue = "-q normal -R'select[mem>5000] rusage[mem=5000]'";

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

# bsub_ssaha_job($queue,$input_file,$subject);

my $call = "bsub -P ensembl-variation -K -w 'done($input_file\_ssaha_job)' -J waiting_process sleep 1"; #waits until all ssaha jobs have finished to continue
# system($call);

parse_ssaha2_out ();

sub bsub_ssaha_job {
  
  my ($queue, $input_file, $subject) = @_ ;

  my $ssaha_command = "$SSAHA_BIN";
  #$ssaha_command .= " -align 0 -kmer 2 -seeds 2 -cut 1000 -output vulgar -depth 10 -best 1 $input_dir/test.fa $input_file";
  $ssaha_command .= " -align 1 -kmer 12 -seeds 5 -cut 1000 -output vulgar -depth 10 -best 1 -save $subject $input_file";
  my $call = "bsub -J $input_file\_ssaha_job $queue -e $output_dir/error_ssaha -o $output_file $ssaha_command";
  system ($call);
}


sub parse_ssaha2_out {

#   open OUT, ">$output_dir/mapping_file_$input_file_name" or die "can't open output_file : $!\n";

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
    if (@h ==1) {
      ($feature_pairs) = find_results($h[0]);
    }
    elsif ($h[0]->{'score'} > $h[1]->{'score'} and @h>=2) {
      ($feature_pairs) = find_results($h[0]);
    }
    elsif ($h[1]->{'score'} > $h[2]->{'score'} and @h>=3) {
      ($feature_pairs) = find_results($h[0],$h[1]);
    }
    elsif ($h[2]->{'score'} > $h[3]->{'score'} and @h>=4) {
      ($feature_pairs) = find_results($h[0],$h[1],$h[2]);
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

  my (@pairs);
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
    
    my ($tmp_q_start, $tmp_q_end);

    my ($f_q_start,$f_q_end) ;#feature start-end used to make feature_pairs
    if ($q_strand ==1) {
      ($f_q_start,$f_q_end) = ($q_start,$q_end);
    }
    else {
      ($f_q_start,$f_q_end) = ($q_end,$q_start);
    }

    my ($seq_region_name) = split /\-/, $t_id;
<<<<<<< run_lrg_ssaha2.pl
<<<<<<< run_lrg_ssaha2.pl
    print "cheesey seq_region_name is $seq_region_name and t_start is $t_start, t_end is $t_end\n";
    my $slice = $sa->fetch_by_region('chromosome',$seq_region_name,$t_start,$t_end, 1);
=======
    print "seq_region_name is $seq_region_name and t_start is $t_start, t_end is $t_end\n";
    #my $slice = $sa->fetch_by_region('chromosome',$seq_region_name,$t_start,$t_end, 1);
=======

>>>>>>> 1.6
    my $slice = $sa->fetch_by_region('chromosome',$seq_region_name);
<<<<<<< run_lrg_ssaha2.pl
>>>>>>> 1.4
    #print "slice is ",ref($slice),"\n";
    print "slice seq_region_name is ",$slice->seq_region_name,'-',length($slice->seq),'-',$slice->length,'-',$slice->start,'-',$slice->end,"\n";
    #print "slice seq is ",$slice->seq,"\n";
=======

>>>>>>> 1.6
    my $q_seqobj = $rec_seq{$q_id};

    #warning that query_seq length != matched query_length
    if ($q_seqobj and length($q_seqobj->seq) != abs($q_end - $q_start) + 1){
      print "q_start is $q_start and q_end is $q_end and length is ",length($q_seqobj->seq),"\n";
      die("query sequence length not equal matched query sequence length");
    }

    my @match_components = split /\s+/, $match;

    my ($full_match,$new_q_start, $new_q_end, $new_t_start, $new_t_end);

    $full_match=1; #initial full_match

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
	my $tmp_qs = ($q_strand ==1) ? $q_start : $new_q_start;
	my $tmp_qe = ($q_strand ==1) ? $new_q_start : $q_start;
	push @pairs, ['DNA',$tmp_qs,$tmp_qe,$t_start,$new_t_start,$q_strand];

	#print "$q_start,$new_q_start,$t_start,$new_t_start\n";
	my $t_seq = substr($slice->seq,$t_start-1,$new_t_start-$t_start+1);
	#print "$q_start,$new_q_start,$t_start,$new_t_start and length of q_seq ",length($q_seq),'-',length($t_seq),"\n";
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
	    #if there is a mismatch, we need to record the base based in query sequence
            if ($q_strand==-1) {
              reverse_comp(\$q_seqs{$count});
              reverse_comp(\$t_seqs{$count});
            }
	    print "count is $count\n";
	    $full_match =0;
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
	    push @pairs, [$type,$tmp_q_start,$tmp_q_end,$sub_t_start,$sub_t_end,$q_strand,uc($q_seqs{$count}),uc($t_seqs{$count})];
            print "in mis-match,q_seq is $q_seqs{$count} and t_seq is $t_seqs{$count}\n";
	    $sub_t_start = $sub_t_end+2;
	  }
	}

	$q_start = $sub_q_start;
	$t_start = $sub_t_start;
	
	#print "q_start is $q_start and t_start is $t_start\n";
        #print "in M, q_start is $q_start and q_end is $new_q_start and t_start is $t_start and t_end is $new_t_start\n";

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
        my $q_seq;
	if ($q_strand ==1) {
	  $new_q_start = $q_start + $query_match_length + 1;
          $q_seq = substr($q_seqobj->seq,$q_start,$query_match_length);
          $q_seq = ($q_seq) ? $q_seq : '-';
	}
	else {
	  $new_q_start = $q_start - $query_match_length - 1;
          $q_seq = substr($q_seqobj->seq,$new_q_start,$query_match_length);
          #reverse_comp(\$q_seq) if $q_seq;#query seq keep same, change target seq
          $q_seq = ($q_seq) ? $q_seq : '-';

	}
	$new_t_start = $t_start + $target_match_length +1;
        my $t_seq = substr($slice->seq, $t_start, $target_match_length);
        if ($q_strand ==-1) {
          reverse_comp(\$t_seq) if $q_seq;
        }
        $t_seq = ($t_seq) ? $t_seq : '-';

	if ($q_strand ==-1) {
	  $tmp_q_start = $new_q_start+1;
	  $tmp_q_end = $q_start-1;
	}
	else {
	  $tmp_q_start = $q_start+1;
	  $tmp_q_end = $new_q_start-1;
	}
	($tmp_q_start,$tmp_q_end) = ($tmp_q_end,$tmp_q_start) if($tmp_q_end<$tmp_q_start);
	push @pairs, [$type,$tmp_q_start,$tmp_q_end,$t_start+1,$new_t_start-1,$q_strand,uc($q_seq),uc($t_seq)];

	print "in G, q_start is $q_start and q_end is $new_q_start and t_start is $t_start and t_end is $new_t_start and q_seq is $q_seq and t_seq is $t_seq\n";
      }

      $q_start = $new_q_start;
      $t_start = $new_t_start;
    }

    $done{$q_id} =1 if $pairs[0]->[0];


    my $fp = Bio::EnsEMBL::FeaturePair->new(-start    => $f_q_start,
					    -end      => $f_q_end,
					    -strand   => $q_strand,
					    -display_id => $h->{'q_id'},
					    -hstart   => $h->{'t_start'},
					    -hend     => $h->{'t_end'},
					    -hstrand  => $h->{'t_strand'},
					    -hseqname => $h->{'t_id'},
					    -slice    => $slice,
					   );
    $fp->seqname($h->{'q_id'});
    $fp->type(\@pairs);
    $fp->identical_matches($full_match);

    push @fps, $fp;

  }

  return \@fps;
}

sub seqname {

  my $self = shift;
  $self->{'seqname'} = shift if(@_);
  return $self->{seqname};

}

sub type {

  my $self = shift;
  $self->{'type'} = shift if(@_);
  return $self->{type};

}

sub identical_matches {

  my $self = shift;
  $self->{'identical_matches'} = shift if(@_);
  return $self->{'identical_matches'};

}

sub get_annotations {
  my $feature_pairs = shift;
  my $lrg_name = 'LRG3';

  foreach my $fp (@$feature_pairs) {
    my ($q_start,$q_end,$t_start,$t_end,$q_strand);
    my $slice = $fp->slice;
    my $q_seqobj = $rec_seq{$fp->seqname};
    my $q_seq = $q_seqobj->seq;
    $q_start = $fp->start;
    $q_end = $fp->end;
    $t_start = $fp->hstart;
    $t_end = $fp->hend;
    $q_strand = $fp->strand;
    my $pairs = $fp->type;

    my $sub_slice = $slice->sub_Slice($t_start,$t_end);
    my $full_match = $fp->identical_matches;

    if ($full_match) {
      $sub_slice->seq_region_name($lrg_name);
      my @genes = @{$sub_slice->get_all_Genes()};
      my @transcripts = @{$genes[0]->get_all_Transcripts()};
      print "gene name is ", $genes[0]->stable_id,'-',$genes[0]->start,'-',$genes[0]->end,"\n" if defined $genes[0];
      print "trans name is ", $transcripts[0]->stable_id,'-',$transcripts[0]->start,'-',$transcripts[0]->end,"\n" if defined $transcripts[0];
      my $hnum_exons;

      foreach my $exon (@{$transcripts[0]->get_all_Exons }) {
	print "  ", $exon->stable_id,"\t",$exon->start,," ", $exon->end, "\n";
	$hnum_exons++;
      }
      print "There are ", $hnum_exons," exons\n";
    }
    else {
      my $csa = $dbCore->get_CoordSystemAdaptor();
      my $cs = $csa->fetch_all_by_name('LRG');
      my $cs_id = $cs->[0]->dbID();
      my $q_seq_region_id;
      my $q_seq_length = length($q_seq);
      my $t_seq_region_id = $sa->get_seq_region_id($slice);
      my $lrg_name_ref =$dbCore->dbc->db_handle->selectall_arrayref(qq{SELECT name from seq_region WHERE name="$lrg_name"});
      if (! $lrg_name_ref->[0][0]) {
	$dbCore->dbc->do(qq{INSERT INTO seq_region(name,coord_system_id,length)values("$lrg_name",$cs_id,$q_seq_length)});
	$q_seq_region_id = $dbCore->dbc->db_handle->{'mysql_insertid'};
      }
      else {
	my $q_seqid_ref = $dbCore->dbc->db_handle->selectall_arrayref(qq{SELECT seq_region_id from seq_region WHERE name="$lrg_name"});
	$q_seq_region_id = $q_seqid_ref->[0][0];
      }
      #my $dna_ref = $dbCore->dbc->db_handle->selectall_arrayref(qq{SELECT seq_region_id from dna WHERE seq_region_id = $q_seq_region_id});
      #if (! $dna_ref->[0][0]) {
	#$dbCore->dbc->do(qq{INSERT INTO dna(seq_region_id,sequence)values($q_seq_region_id,"$q_seq")});
      #}

      foreach  my $pair (sort {$a->[3]<=>$b->[3]} @$pairs) {
	print "pairs are ",$pair->[0],'-',$pair->[1],'-',$pair->[2],'-',$pair->[3],'-',$pair->[4],'-',$pair->[5],'-',$pair->[6],'-',$pair->[7],"\n";

	if ($pair->[0] eq 'DNA') {
	  if ($pair->[2]-$pair->[1] == $pair->[4]-$pair->[3]) {
	    $dbCore->dbc->do(qq{INSERT IGNORE INTO assembly(asm_seq_region_id,cmp_seq_region_id,asm_start,asm_end,cmp_start,cmp_end,ori)values($t_seq_region_id,$q_seq_region_id,$pair->[3],$pair->[4],$pair->[1],$pair->[2],$q_strand)});
	  }
	  else {

	    die("distance between query and target is not the same, there is a indel");
	  }
	}
      }

      my $cs1 = $csa->fetch_by_name("Chromosome","NCBI36");
      my $cs2 = $csa->fetch_by_name("LRG");


      my $asm = $asma->fetch_by_CoordSystems($cs1,$cs2);

      #foreach my $name ($lrg_name){
	my $asm = $asma->fetch_by_CoordSystems($cs1,$cs2);
	$asm->flush;
	my $lrg_slice = $sa->fetch_by_region("LRG",$lrg_name);
  
	print "$name start = ".$lrg_slice->start."\tend= ".$lrg_slice->end."\n";
  
	#foreach my $gene (@{$lrg_slice->get_all_Genes}){
	#  print "\tLRG\t".$gene->stable_id."\n";
	#}

      ###needs projection to make transfer to work!!!

	my $min = 99999999;
	my $max = -9999999;
	my $chrom;
	my $strand;

	foreach my $segment (@{$lrg_slice->project('chromosome')}) {
	  my $from_start = $segment->from_start();
	  my $from_end    = $segment->from_end();
	  my $to_name    = $segment->to_Slice->seq_region_name();
	  $chrom = $to_name;

	  my $to_start    = $segment->to_Slice->start();
	  my $to_end    = $segment->to_Slice->end();
	  if($to_start > $max){
	    $max = $to_start;
	  }
	  if($to_start < $min){
	    $min = $to_start;
	  }
	  if($to_end > $max){
	    $max = $to_end;
	  }
	  if($to_end <  $min){
	    $min = $to_end;
	  }
	  my $ori        = $segment->to_Slice->strand();
	  $strand = $ori;   
    
	  print "$from_start-$from_end  => $to_name $to_start-$to_end ($ori) \n";
	}
	
	print "################\n";
	print "STUFF STARTS NOW\n";
	print "################\n";


	##initial run, there should be only one LRG covering the region
#       my $hslice = $sa->fetch_by_region('LRG',"$lrg_name");
#       print "hslice_start ",$hslice->start," and hslice_end ",$hslice->end,"\n";
#       foreach my $gene (@{$hslice->get_all_Genes()}){
# 	print "gene_name is ",$gene->stable_id,"\n";
#       }
#       print "q_seq from hslice is ",$hslice->start,'-'.$hslice->end,"\n";
#       print "length of q_seq is ",length($q_seq), " and length of hseq is ", length($hslice->seq),"\n";

#       if ($hslice->seq eq $q_seq) {
# 	print "hseq is same as q_seq\n";
#       }
#       else {
# 	print "hseq is different from q_seq\n";
# 	print "q_seq is ",$q_seq,"\n";
# 	print "hseq is ",$hslice->seq,"\n";
#       }

=head
      my $msc = Bio::EnsEMBL::MappedSliceContainer->new(
							-SLICE => $hslice
						       );
      my $asa = $dbCore->get_AssemblySliceAdaptor();
      $msc->set_AssemblySliceAdaptor($asa);
      #$msc->attach_AssemblySlice('NCBI36');
      $msc->attach_LRG('LRG');

      foreach my $mapped_slice (@{$msc->get_all_MappedSlices()}){
	print "mapped_name is ",$mapped_slice->name,"\n";
	#$mapped_slice->get_all_Genes;


	foreach my $mapped_slice (@{$msc->get_all_MappedSlices()}){
	  print "mapped_name is ",$mapped_slice->name,"\n";
	  if ($mapped_slice->seq eq $q_seq) {
	    print "mapped_seq is same as q_seq\n";
	  }
	  else {
	    print "mapped_seq is different from q_seq\n";
	    #print "mapped_seq is ",$mapped_slice->seq,"\n";
	  }

	  print "mapped_slice name is ",$mapped_slice->name,"\t",$mapped_slice->start,"\t",$mapped_slice->end,"\n";
	  my $num_exons = 0;
	  my $num_trans = 0;

	  foreach my $exon (@{ $mapped_slice->get_all_Exons }) {#working here???
	    print "  ", $exon->stable_id,"\t",$exon->start,," ", $exon->end, "\n";
	    $num_exons++;
	  }
	  print "There are ", $num_exons," exons\n"; 
	


	  #need mapped_slice->get_all_Exons to work for multi seq_region_name defined as LRG
	  foreach my $gene (@{ $mapped_slice->get_all_Genes }) {#not working here
	    print "gene is ",ref($gene),"\n";
	    my @gene_dbentries = @{$gene->get_all_DBEntries()};
	    foreach my $dbe (@gene_dbentries) {
	      my @gene_syns = @{$dbe->get_all_synonyms()};
	      print "all_synonyms are @gene_syns\n";
	      print $dbe->db_display_name(),'-',$dbe->description(),'-',$gene->external_db(),"\n";
	    }
	    foreach my $transcript (@{$mapped_slice->get_all_Transcripts()}) {
	      print "Tanscript cdna_start_end  ", $transcript->stable_id,"\t",$transcript->cdna_coding_start,," ", $transcript->cdna_coding_end, "\n";
	      $num_trans++;
	      print "There are ", $num_trans," transcripts\n";
	      my $translation = $transcript->translation();
	      print "The protein sequence is ",$translation->seq,"\n";
	      my @transl_dbentries = @{$translation->get_all_DBEntries()};
	      foreach my $dbe (@transl_dbentries) {
		print $dbe->db_display_name(),'-',$dbe->get_all_synonyms(),'-',$dbe->description(),'-',$translation->external_db,"\n";
	      }
	    }
	  }
	  #=cut
	}
      }
=cut

        #my $min = 45616138;
        #my $max = 45638999;
        #my $strand = -1;
        #my $chrom = 17;
        my $slice = $sa->fetch_by_region("chromosome",$chrom, $min, $max, $strand);
	print "Slice ".$slice->seq_region_name."\t".$slice->start."\t".$slice->end."\n";
	foreach my $gene (@{$slice->get_all_Genes}){

      #my @genes = @{$sub_slice->get_all_Genes()};
      #foreach my $gene (@genes) {
	my $num_exons = 0;
	my $num_trans = 0;
	
	print "before transfprm g start-end ",$gene->start,'-',$gene->end,"\n";
	#my $new_gene = $gene->transform('LRG');
	my $new_gene = $gene->transfer($lrg_slice);
	print "after transform g start-end ",$new_gene->start,'-',$new_gene->end,"\n" if $new_gene;
	my @transcripts = @{$new_gene->get_all_Transcripts()};
	print "trs start-end ",$transcripts[0]->start,'-',$transcripts[0]->end,"\n" if $transcripts[0];
	print "new_gene name is ", $new_gene->stable_id,'-',$new_gene->start,'-',$new_gene->end,"\n" if defined $new_gene;
	#my @genes = @{$hslice->get_all_Genes()};
	my @transcripts = @{$new_gene->get_all_Transcripts()} if $new_gene;
	#print "gene name is ", $genes[0]->stable_id,'-',$genes[0]->start,'-',$genes[0]->end,"\n" if defined $genes[0];
	print "trans name is ", $transcripts[0]->stable_id,'-',$transcripts[0]->start,'-',$transcripts[0]->end,"\n" if defined $transcripts[0];
	my $hnum_exons;

	foreach my $exon (@{$transcripts[0]->get_all_Exons }) {
	  print "  ", $exon->stable_id,"\t",$exon->start,," ", $exon->end, "\n";
	  $hnum_exons++;
	}
	print "There are ", $hnum_exons," exons\n";

	foreach my $transcript (@{$new_gene->get_all_Transcripts()}) {
	    print "Tanscript cdna_start_end  ", $transcript->stable_id,"\t",$transcript->cdna_coding_start,," ", $transcript->cdna_coding_end, "\n";
	    $num_trans++;
	    print "There are ", $num_trans," transcripts\n";
	    my $translation = $transcript->translation();
	    print "The protein sequence is ",$translation->seq,"\n";
	    my @transl_dbentries = @{$translation->get_all_DBEntries()};
	    foreach my $dbe (@transl_dbentries) {
	      print $dbe->db_display_name(),'-',$dbe->get_all_synonyms(),'-',$dbe->description(),"\n";
	    }
	  }
      }
    }
  }
}
#}


