#! /usr/local/ensembl/bin/perl

use strict;
use Getopt::Long;

my ($start, $end, $chr, $output_dir, $input_dir, $target_file, $run, $parse);

GetOptions('start=i'         => \$start,
	   'end=i'           => \$end,
           'output_dir=s'      => \$output_dir,   
	   'input_dir=s'       => \$input_dir,
	   'target_file=s'      => \$target_file,
           'run'             => \$run,
           'parse'           => \$parse,
           );


$output_dir ||=".";
my $seed;
###kmer is always set to 12, so ajust seed for different species
$seed = 2 if $input_dir =~ /zfish/;
##if input from array, this is short sequence, run it in turing (big mem) or use exonerate(small memory 1024 MB)
##/software/ensembl/bin/exonerate-1.4.0 --query /lustre/work1/ensembl/yuan/array/GenomeWideSNP_6/mapping/input_dir/1_query_seq --target /lustre/blastdb/Ensembl/Human/NCBI36/genome/softmasked_dusted.fa --showalignment no --showvulgar yes
##nathan's affy array parameters to ensure one mismatch : --bestn 100 --dnahspthreshold 116 --fsmmemory 256 --dnawordlen 25 --dnawordthreshold 11
##my affy parameters to ensure full match : /software/ensembl/bin/exonerate-1.4.0 --query /lustre/work1/ensembl/yuan/array/GenomeWideSNP_6/mapping/input_dir/1_query_seq --target /lustre/blastdb/Ensembl/Human/NCBI36/genome/softmasked_dusted.fa --showalignment yes --showvulgar no --ryo "vulgar: %S %V %pi\n" --dnahspthreshold 156
##dnahspthreshold 116 for 25 mer 156 for 33 mer for exact match test it for each run

$seed = "2 -kmer 6 -ckmer 6 -cmatch 10 -tags 1 -score 12 -skip 1 -sense 1" if $input_dir =~ /array/;
$seed = 5 if $input_dir !~ /zfish|array/;
print "seed is $seed\n";

my $output_file = "ssaha2_output";
my $queue_long = "-q long -M7000000 -R'select[mem>7000] rusage[mem=7000]'";
my $queue_normal = "-q long -M5000000 -R'select[mem>5000] rusage[mem=5000]'"; ##for new farm


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
  #print "Submitting ssaha2Build job...\n";
  #my $ssaha2build = "bsub $queue_long -J 'ssaha2build' -o $target_file\_out /nfs/users/nfs_y/yuan/ensembl/src/ensembl-variation/scripts/ssahaSNP/ssaha2/ssaha2_v1.0.9_x86_64/ssaha2Build -save $subject $target_file";
  #system("$ssaha2build");

  #my $call = "bsub -q normal -K -w 'done('ssaha2build')' -J waiting_process sleep 1"; #waits until all variation features have finished to continue

  #system($call);

  if ($start and $end) {#mainly for flanking_sequence mapping
    #for (my $start=$start; $start<=$end; $start++) {
    #input name will change to \\\$LSB_JOBINDEX\_query_seq for flanking sequence mapping
      my $input = "$input_dir/Watson_\\\$LSB_JOBINDEX\.fastq" if $input_dir;
      next if (-z "$input");
      if (! -e "$output_dir/ssaha_out_$start") {#check for ssaha2_output is exists or not
	print "Submitting ssaha2 job...\n";
	bsub_ssaha_job_array($start,$end,$input,$queue,$target_file);
      }
    #}
  }
  else {
    while (<$input_dir/*fastq>) { #mainly for sequencing reads
      my $input = $_;
      my ($input_file) = $input =~ /^.*\/(.*)$/;

      my $num = `grep "^\@" $input |grep -v -c "[?%*]"`;
      $num =~ s/^\s+|\s+$//g;
      print "num is $num again\n";

      my $n;
      #my $size = $num;#this is run whole file in one go
      my $size = 50000;		#this is to split file in this size
      my $count;

      for ($n = 1;$n<=$num;$n+=$size) {
	my $end = $n+$size-1;
	$end = $num if ($end > $num);
	$count++;
	if (! -e "$output_dir/ssaha_out_$input_file\:$n") { #check for ssaha2_output is exists or not
	  print "Submitting ssaha2 job $count ...\n";
	  bsub_ssaha_job ("$input_file\:$n",$end,$input,$queue,$target_file);
	}
      }
    }
  }
}


my $call = "bsub -K -w 'done(ssaha_out*)' -J waiting_process sleep 1"; #waits until all ssaha jobs have finished to continue
system($call);

sub parse_ssaha2 {

  if ($start and $end) {
    for (my $start=$start; $start<=$end; $start++) {
      my $input = "$input_dir/$start\_query_seq" if $input_dir;
      next if (-z "$output_dir/ssaha_out_$start");
      parse_ssaha2_out ($start, $input,"$output_dir/ssaha_out_$start");
    }
  }
}

sub bsub_ssaha_job_array {
  my ($start,$end, $input_file, $queue, $target_file) = @_ ;

  my ($subj_dir,$tar_file) = $target_file =~ /^(.*)\/(.*)$/;
  my ($subname) = $tar_file =~ /^(.*)\.*.*$/;
  my $subject = "$subj_dir/$subname";

  #for 454 reads
  my $ssaha_command = "/nfs/users/nfs_y/yuan/ensembl/src/ensembl-variation/scripts/ssahaSNP/ssaha2/pileup_v0.5/ssaha2/ssaha2-2.3_x86_64 -rtype 454 -best 1 -output cigar -name -save $subject $input_file";

  my $call = "bsub -J'ssaha_out_[$start-$end]%2' $queue -e $output_dir/ssaha-%J.%I.err -o $output_dir/ssaha-%J.%I.out ";

  $call .= " $ssaha_command";
  
  system ($call);
  print "$call\n";
}

sub bsub_ssaha_job {
  
  my ($start, $end, $input_file, $queue, $target_file) = @_ ;

  my ($input_file_name,$n) = split /\:/, $start if $start =~ /\:/;
  
  my ($subj_dir,$tar_file) = $target_file =~ /^(.*)\/(.*)$/;
  my ($subname) = $tar_file =~ /^(.*)\.*.*$/;
  my $subject = "$subj_dir/$subname";

  my ($ssaha_command, $count);

  print "target_file is $target_file\n";

  #for normal mapping with long sequence (more than 2 kb)
  #$ssaha_command = "/nfs/users/nfs_y/yuan/ensembl/src/ensembl-variation/scripts/ssahaSNP/ssaha2/ssaha2_v1.0.9_x86_64/ssaha2 -align 0 -kmer 12 -seeds $seed -cut 5000 -output vulgar -depth 5 -best 1 -tags 1 -save $subject $input_file";
  #for abi reads
  #$ssaha_command = "/nfs/users/nfs_y/yuan/ensembl/src/ensembl-variation/scripts/ssahaSNP/ssaha2/pileup_v0.5/ssaha2/ssaha2-2.3_x86_64 -start $n -end $end -seeds 15 -score 250 -tags 1 -best 1 -output cigar -name -memory 300 -cut 5000 -save $subject $input_file";
  #for 454 reads
  $ssaha_command = "/nfs/users/nfs_y/yuan/ensembl/src/ensembl-variation/scripts/ssahaSNP/ssaha2/pileup_v0.5/ssaha2/ssaha2-2.3_x86_64 -start $n -end $end -454 -best 1 -output cigar -name -save $subject $input_file";
  #for watson reads
  #$ssaha_command = "/nfs/users/nfs_y/yuan/ensembl/src/ensembl-variation/scripts/ssahaSNP/ssaha2/pileup_v0.5/ssaha2/ssaha2-2.3_x86_64 -start $n -end $end -seeds 5 -cut 5000 -memory 300 -best 1 -output cigar -name -save $subject $input_file";
  #for sam output
  #$ssaha_command = "/nfs/users/nfs_y/yuan/ensembl/src/ensembl-variation/scripts/ssahaSNP/ssaha2/ssaha2_x86_64.bin -start $n -end $end -seeds 15 -score 250 -tags 1 -best 1 -output sam -name -memory 300 -cut 5000 -save $subject $input_file";
  print "job_num is ", ++$count, " and start is $start and out is ssaha_out_$start\n";
  my $call = "bsub -J $input_dir\_ssaha_job_$start $queue -e $output_dir/error_ssaha_$start -o $output_dir/ssaha_out_$start ";

  $call .= " $ssaha_command";
  
  system ($call);
  print "$call\n";
}




sub parse_ssaha2_out {

  my ($start, $input_file, $ssaha2_output) = @_;

  print "input_file is $input_file ssaha2_output is $ssaha2_output\n";

  open INPUT, "$input_file" or die "can't open $input_file : $!\n";
  open OUT, ">$output_dir/mapping_file_$start" or die "can't open output_file : $!\n";

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
    my ($fseq,$lseq) = split /[NWACTG\-\_]+/,$rec_seq{$name};
    my $snp_posf = length($fseq)+1;
    my $snp_posl = length($lseq)+1;
    my $tot_length = $snp_posf + $snp_posl-2;
    my $tot_length = length($rec_seq{$name});
    $snp_pos{$name} = "$snp_posf\_$snp_posl";
    $snp_pos{$name} = "$snp_posf\_$snp_posl";
    $input_length{$name} = $tot_length;
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
  }

  my ($total_seq,$no);

  $total_seq = keys %rec_seq;

  foreach my $q_id (keys %rec_seq) {
    if (!$done{$q_id}) {
      $no++;
    }
  }

  print "$no out of $total_seq are not mapped\n";

  close OUT;
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
                                              
      #print "$snp_q_pos,$q_start,$new_q_start and $t_start\n";
      if ($snp_q_pos >= $q_start and $snp_q_pos <= $new_q_start or $snp_q_pos >= $new_q_start and $snp_q_pos <= $q_start) {
	if ($type eq 'M') {
	  $snp_t_start = $t_start + abs($snp_q_pos-$q_start);
	  $snp_t_end = $snp_t_start;
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

	
	
