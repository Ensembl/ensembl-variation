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
use Getopt::Long;

my ($start, $end, $chr, $output_dir, $input_dir, $target_file, $run, $parse, $split, $rerun, $ssahabuild, $ssaha2);

GetOptions('start=i'         => \$start,
	   'end=i'           => \$end,
           'output_dir=s'      => \$output_dir,   
	   'input_dir=s'       => \$input_dir,
	   'target_file=s'      => \$target_file,
           'run'             => \$run,
           'parse'           => \$parse,
		   'ssahabuild=s' => \$ssahabuild, # location of ssahaBuild binary
		   'ssaha2=s' => \$ssaha2, # location of ssaha2 binary
	   'split'           => \$split,  #option to split reads in groups to reduce run time
	   'rerun=i'           => \$rerun, #to rerun it if some jobs failed i is the split start number
           );


$output_dir ||=".";
my $seed;
###kmer is always set to 12, so ajust seed for different species
$seed = 2 if $input_dir =~ /zfish/;
my $skip = 3;

##if input from array, this is short sequence, run it in turing (big mem) or use exonerate(small memory 1024 MB)
##[exonerate_bin_dir]/exonerate-1.4.0 --query [ssaha_input_dir]/1_query_seq --target [blast_db_dir]/softmasked_dusted.fa --showalignment no --showvulgar yes
##nathan's affy array parameters to ensure one mismatch : --bestn 100 --dnahspthreshold 116 --fsmmemory 256 --dnawordlen 25 --dnawordthreshold 11
##my affy parameters to ensure full match : [exonerate_bin_dir]/exonerate-1.4.0 --query [ssaha_input_dir]/1_query_seq --target [blast_db_dir]/softmasked_dusted.fa --showalignment yes --showvulgar no --ryo "vulgar: %S %V %pi\n" --dnahspthreshold 156
##dnahspthreshold 116 for 25 mer 156 for 33 mer for exact match test it for each run

$seed = "2 -kmer 6 -ckmer 6 -cmatch 10 -tags 1 -score 12 -skip 1 -sense 1" if $input_dir =~ /array/;
$seed = 5 if $input_dir !~ /zfish|array/;
print "seed is $seed\n";

my $output_file = "ssaha2_output";
my $queue_long = "-q long -M10000000 -R'select[mem>10000] rusage[mem=10000]'";## for watson reads
my $queue_normal = "-q normal -M6000000 -R'select[mem>6000] rusage[mem=6000]'"; ##for new farm abi reads


#my $queue = $queue_long;
my $queue = $queue_normal;

my (@chr_names, %snp_pos, %input_length, %done);

run_ssaha2() if $run;
parse_ssaha2() if $parse;

sub run_ssaha2 {

  my ($subj_dir,$tar_file) = $target_file =~ /^(.*)\/(.*)$/;
  my ($subname) = $tar_file =~ /^(.*)\.*.*$/;
  my $subject = "$subj_dir/$subname";
  print "tar_file is $tar_file and subject is $subject and target_file is $target_file\n";

  if (! -f "$subject\.body") {
    print "Submitting ssaha2Build job...\n";
    my $ssaha2build = "bsub $queue_long -J 'ssaha2build' -o $target_file\_out $ssahabuild -kmer 12 -skip $skip -save $subject $target_file";
    system("$ssaha2build");

    my $call = "bsub -q normal -K -w 'done('ssaha2build')' -J waiting_process sleep 1"; #waits until all variation features have finished to continue
    system($call);
		
  }

  if ($start and $end) {
    #input name will change to \\\$LSB_JOBINDEX\_query_seq for flanking sequence mapping
    #my $input = "$input_dir/Watson_\\\$LSB_JOBINDEX\.fastq" if $input_dir;
    my $input = "$input_dir/\\\$LSB_JOBINDEX\.query_seq" if $input_dir;#need a dot in input, not a '_' after LSB_JOBINDEX
    #next if (-z "$input");
    print "input is $input and split $split\n";
    if (! $split) {#mainly for flanking_sequence mapping
      print "Submitting ssaha2 job array...\n";
      bsub_ssaha_job_array($start,$end,$input,$queue,$target_file);
    }
    else {
      my $input;
      my @files = glob("$input_dir/*fastq") ; #mainly for sequencing reads
      $input = $files[0];

      my ($input_file) = $input =~ /^.*\/(.*)$/;
	
      my $num = `grep "^\@" $input |grep -v -c "[?%*]"`;
      $num =~ s/^\s+|\s+$//g;
      print "number of sequence entries in file is $num \n";

      my $n;
      #my $size = $num;#this is run whole file in one go
      my $size = 50000;		#this is to split file in this size
      my $count;
      $input = "$input_dir/Watson_\\\$LSB_JOBINDEX\.fastq";

      if (!$rerun) {
	for ($n = 1;$n<=$num;$n+=$size) {
	  my $e = $n+$size-1;
	  $e = $num if ($e > $num);
	  $count++;
	  if (! -e "$output_dir/ssaha_out_$input_file\:$n") { #check for ssaha2_output is exists or not
	    print "Submitting ssaha2 job $count ...\n";
	    bsub_ssaha_job ($start,$end,$n,$e,$queue,$input,$target_file);
	  }
	}
      }
      else {
	$n = $rerun;
	my $e = $n+$size-1;
	$e = $num if ($e > $num);
	bsub_ssaha_job ($start,$end,$n,$e,$queue,$input,$target_file);
      }
    }
  }
}


my $call = "bsub -K -w 'done(ssaha_out*)' -J waiting_process sleep 1"; #waits until all ssaha jobs have finished to continue
#system($call);

sub parse_ssaha2 {

  foreach my $num ($start..$end) {
    my $input = "$input_dir/$num\.query_seq" if $input_dir;
    my $out_file = "$output_dir/ssaha.$num\.out";
    print "input is $input and out_file is $out_file\n";
    next if (-z "$out_file");
    parse_ssaha2_out ($num, $input,"$out_file");
  }
}


sub bsub_ssaha_job_array {
  my ($start,$end, $input_file, $queue, $target_file) = @_ ;

  my ($subj_dir,$tar_file) = $target_file =~ /^(.*)\/(.*)$/;
  my ($subname) = $tar_file =~ /^(.*)\.*.*$/;
  my $subject = "$subj_dir/$subname";

  #for 454 reads
  #my $ssaha_command = "$ssaha2 -rtype 454 -best 1 -output cigar -name -save $subject $input_file";
  #for watson for comparison
  #my $ssaha_command = "$ssaha2 -seeds 5 -cut 5000 -memory 300 -best 1 -output cigar -name -save $subject $input_file";

  #for normal flanking mapping
  my $ssaha_command = "$ssaha2 -align 0 -kmer 12 -skip $skip -seeds $seed -cut 5000 -output vulgar -depth 5 -best 1 -tags 1 -name -save $subject $input_file";
  my $call = "bsub -J'ssaha_out_[$start-$end]%50' $queue -e $output_dir/ssaha.%I.err -o $output_dir/ssaha.%I.out ";
  $call .= " $ssaha_command";
  print "$call\n";
  system ($call);

}

sub bsub_ssaha_job {

  my ($start, $end, $n, $e, $queue, $input_file, $target_file) = @_ ;
	

  my ($subj_dir,$tar_file) = $target_file =~ /^(.*)\/(.*)$/;
  my ($subname) = $tar_file =~ /^(.*)\.*.*$/;
  my $subject = "$subj_dir/$subname";

  my ($ssaha_command, $count);

  print "target_file is $target_file\n";

  #for normal mapping with long sequence (more than 2 kb)
  #$ssaha_command = "$ssaha2 -align 0 -kmer 12 -seeds $seed -cut 5000 -output vulgar -depth 5 -best 1 -tags 1 -save $subject $input_file";
  #for abi reads
  #$ssaha_command = "$ssaha2 -start $n -end $end -seeds 5 -score 250 -tags 1 -best 1 -output cigar -name -memory 300 -cut 5000 -save $subject $input_file";
  #for 454 reads -454 = "-skip 5 -seeds 2 -score 30 -sense 1 -cmatch 10 -ckmer 6" (use skip 5 instead of skip 3 to save memory)
  $ssaha_command = "$ssaha2 -start $n -end $e -skip 5 -seeds 2 -score 30 -sense 1 -cmatch 10 -ckmer 6 -best 1 -output cigar -name -save $subject $input_file";
  #for watson reads
  #$ssaha_command = "$ssaha2 -start $n -end $e -seeds 5 -cut 5000 -memory 300 -best 1 -output cigar -name -save $subject $input_file";
  #for sam output
  #$ssaha_command = "$ssaha2 -start $n -end $end -seeds 15 -score 250 -tags 1 -best 1 -output sam -name -memory 300 -cut 5000 -save $subject $input_file";
  print "job_num is ", ++$count, " and start is $start and out is ssaha_out_$start\n";
  #my $call = "bsub -J $input_dir\_ssaha_job_$start $queue -e $output_dir/error_ssaha_$start -o $output_dir/ssaha_out_$start ";

  my $call = "bsub -J'ssaha_out_[$start-$end]%50' $queue -e $output_dir/ssaha.%I:$n.err -o $output_dir/ssaha.%I:$n.out ";

  $call .= " $ssaha_command";
  
  system ($call);
  print "$call\n";
}




sub parse_ssaha2_out {

  my ($start, $input_file, $ssaha2_output) = @_;

  print "input_file is $input_file ssaha2_output is $ssaha2_output\n";

  open INPUT, "$input_file" or die "can't open $input_file : $!\n";
  open OUT, ">$output_dir/mapping_file_$start" or die "can't open output_file : $!\n";
  open OUT1, ">$output_dir/out_mapping_$start"  or die "can't open out_mapping file : $!\n";

  my ($name,%rec_seq,%rec_find);

  while (<INPUT>) {
    if (/^\>/) {
      s/^\>|\n//g;
      $name = $_;
    }
    else {
      s/\n//;
      $rec_seq{$name} .= $_;
    }
  }

  foreach my $name (keys %rec_seq) {
    my ($fseq,$lseq) = split /[W\-\_]+/,$rec_seq{$name};
    my $snp_posf = length($fseq)+1;
    my $snp_posl = length($lseq)+1;
    #my $tot_length = $snp_posf + $snp_posl-2;
    my $tot_length = length($rec_seq{$name});
    $snp_pos{$name} = "$snp_posf\_$snp_posl";
    $snp_pos{$name} = "$snp_posf\_$snp_posl";
    $input_length{$name} = $tot_length;
		#print "$name: ".$snp_pos{$name}.' > '.$rec_seq{$name}."\n";
  }


  open SSAHA, "$ssaha2_output" or die "can't open $output_file : $!\n";
  open MAP, ">>$output_dir/THREE_MAPPINGS";

  while (<SSAHA>) {
    #print ;
    if (/^vulgar\:\s+(\S+)\s+(\d+)\s+(\d+)\s+(\+|\-)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\+|\-)\s+(\d+)\s+(.*)$/) {
      #print "q_id is $1 and q_start is $2 and q_end is $3 and q_strand is $4 and t_id is $5 and t_start is $6 and t_end is $7 and t_strand is $8 and score is $9 and match is $10\n";
      #next unless $1 eq "rs31190187";
      my $h={};
      $h->{'q_id'} = $1;
      $h->{'q_start'} = $2;
      $h->{'q_end'} = $3;
      $h->{'q_strand'} = $4;
      $h->{'t_id'} = $5;
      $h->{'t_start'} = $6;
      $h->{'t_end'} = $7;
      $h->{'t_strand'} = $8;
      $h->{'score'} = $9;
      $h->{'match'} = $10;

      push @{$rec_find{$h->{'q_id'}}}, $h if $h->{'q_id'};
    }
  }

  foreach my $q_id (keys %rec_find) {
    my @h = sort {$b->{'score'}<=>$a->{'score'}} @{$rec_find{$q_id}};
    #print "There are ",scalar @h," hits and q_id is $q_id\n";
    
    # Get rid of mappings that score less than the best scoring mapping
    my @mappings = grep {$_->{'score'} == $h[0]->{'score'}} @h;
    
    find_results(@mappings);

    # Comment out the old code checking number of mappings etc. This should be done by the post-processing scripts
=head
    if (scalar @h==1) {
      find_results($h[0]);
    }
    elsif ($h[0]->{'score'} > $h[1]->{'score'} and @h=>2) {
      find_results($h[0]);
    }
    elsif ($h[0]->{'score'} == $h[1]->{'score'} and @h==2) {
      find_results($h[0],$h[1]);
    }
    elsif ($h[1]->{'score'} > $h[2]->{'score'} and @h=>3) {
      find_results($h[0],$h[1]);
    }
    elsif ($h[1]->{'score'} == $h[2]->{'score'} and @h==3) {
      find_results($h[0],$h[1],$h[2]);
    }
    elsif ($h[2]->{'score'} > $h[3]->{'score'} and @h=>4) {
      find_results($h[0],$h[1],$h[2]);
    }
    #for MHC region
    elsif (($h[3]->{'t_id'}=~/6|MHC/ and $h[3]->{'score'} <= $h[2]->{'score'} and @h==4) or ($h[4] and $h[4]->{'t_id'}=~/6|MHC/ and $h[4]->{'score'} <= $h[2]->{'score'} and @h==5) or ($h[5] and $h[5]->{'t_id'}=~/6|MHC/ and $h[5]->{'score'} <= $h[2]->{'score'} and @h==6) or ($h[6] and $h[6]->{'t_id'}=~/6|MHC/ and $h[6]->{'score'} <= $h[2]->{'score'} and @h==7) or ($h[7] and $h[7]->{'t_id'}=~/6|MHC/ and $h[7]->{'score'} <= $h[2]->{'score'} and @h ==8)) {
      find_results(@h);
    }
    else {
      #needs to be written to the failed_variation table
      print MAP "$q_id\n";
    }
=cut
  }

  my ($total_seq,$no);

  $total_seq = keys %rec_seq;

  foreach my $q_id (keys %rec_seq) {
    if (!$done{$q_id}) {
      $no++;
    }
  }

  print OUT1 "$no out of $total_seq are not mapped\n";

  close OUT;
  close OUT1;
  close MAP;
}

sub find_results {

  my (@h) = @_;

  LINE : foreach my $h (@h) {
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

    my @match_components = split /\s+/, $match;
      
    my ($new_q_start, $new_q_end, $new_t_start, $new_t_end);
      
    #print "$q_id,$q_start,$q_end,$q_strand,$t_id,$t_start,$t_end,$t_strand,$score,@match_components\n";
      
    while (@match_components) {
	
      my $type = shift @match_components;
      my $query_match_length = shift @match_components;
      my $target_match_length = shift @match_components;
      
      #print "type is $type and query_match_length is $query_match_length and target_match_length is $target_match_length\n";
      next LINE if (! defined $query_match_length or ! defined $target_match_length);
      my ($snp_q_pos, $snp_t_start, $snp_t_end, $snp_strand, $pre_target_match_length);
      my ($snp_posf,$snp_posl) = split /\_/, $snp_pos{$q_id};
      $snp_q_pos = $snp_posf;
      
      if ($q_strand eq '+') {
        $new_q_start = $q_start + $query_match_length - 1 if ($type eq 'M');
        $new_q_start = $q_start + $query_match_length + 1 if ($type eq 'G');
        $new_t_start = $t_start + $target_match_length -1 if ($type eq 'M');
        $new_t_start = $t_start + $target_match_length +1 if ($type eq 'G');
      }
      elsif ($q_strand eq '-') {
        $new_q_start = $q_start - $query_match_length + 1 if ($type eq 'M');
        $new_q_start = $q_start - $query_match_length - 1 if ($type eq 'G');
        $new_t_start = $t_start + $target_match_length -1 if ($type eq 'M');
        $new_t_start = $t_start + $target_match_length +1 if ($type eq 'G');
      }
                                              
      #print "$q_id: [$snp_q_pos:$q_start,$new_q_start],$new_t_start and $t_start ; ".$snp_pos{$q_id}."\n";
      if (($snp_q_pos >= $q_start and $snp_q_pos <= $new_q_start) or ($snp_q_pos >= $new_q_start and $snp_q_pos <= $q_start)) {
	if ($type eq 'M') {
	  $snp_t_start = $t_start + abs($snp_q_pos-$q_start);
	  $snp_t_end = $snp_t_start;
		#print "$q_id: $snp_q_pos,$q_start => $snp_t_start\n";
          #print "abs is ",abs($snp_q_pos-$q_start),"\n";
          #print "$snp_t_start and $snp_t_end\n";
	}
	elsif ($type eq 'G') {
	  if ($target_match_length ==0) {
	    $snp_t_start = $t_start+1;
	    $snp_t_end = $snp_t_start -1;
	  }
	  elsif ($query_match_length ==0) {
	    $snp_t_start = $t_start+1;
	    $snp_t_end = $snp_t_start + $target_match_length;
	  }
	}
      }
      #print "$snp_t_start && $snp_t_end && $q_strand\n";      
      if ($snp_t_start && $snp_t_end && $q_strand ) {
	#print "$snp_t_start && $snp_t_end && $q_strand\n";
        my $final_score = $score/$input_length{$q_id};
	#print OUT "MORE_HITS\t" if ($h[2]);
	print OUT "$q_id\t$t_id\t$snp_t_start\t$snp_t_end\t$q_strand\t$final_score\n";
	$done{$q_id}=1;
	last;
      }
      
      $q_start = $new_q_start;
      $t_start = $new_t_start;
    }
  }
}

	
	
