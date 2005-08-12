#! /usr/local/ensembl/bin/perl

use strict;
use Getopt::Long;

my ($start, $end, $input_dir, $target_dir);

GetOptions('start=i'         => \$start,
	   'end=i'           => \$end,
	   'input_dir=s'       => \$input_dir,
	   'target_dir=s'      => \$target_dir);


$input_dir ||= "/ecs2/scratch4/yuan/mouse/release32/snp_seq";
$target_dir ||= "/ecs2/scratch4/yuan/mouse/release32/target";

my $seed;
###kmer is always set to 12, so ajust seed for different species
$seed = 3 if $input_dir =~ /zfish/;
$seed = 5 if $input_dir !~ /zfish/;

my $output_file = "ssaha2_output2";
my $queue = "-q normal -R'select[largedata && mem>1000] rusage[mem=1000]'";
my $queue1 = "-q normal -R'select[largedata]'";
my $queue2 = "-q normal"; ##for new farm

$queue = $queue2;

my (@chr_names, %snp_pos, %input_length);

while (<$target_dir/*fasta>) {
  chomp;
  push @chr_names, $_;
}


if ($start and $end) {
  my $input = "$input_dir/$start\_query_seq" if $input_dir;
  if (! -e "$start\_total_ssaha_out" and ! -e "total_ssaha_out") {
    for ($start=$start; $start<=$end; $start++) {
      bsub_ssaha_job ($start,$input,$queue2,\@chr_names);
    }
    
    my $call = "bsub -K -w 'done($input_dir\_ssaha_job*)' -J waiting_process sleep 1"; #waits until all ssaha jobs have finished to continue
    system($call);

    system("cat ssaha_out_* > total_ssaha_out");
    unlink<ssaha_out_*>;
  }

  print "PARSING SSAHA2 OUTPUT...\n";
  
  parse_ssaha2_out ($input,"total_ssaha_out");
  
}


sub bsub_ssaha_job {
  
  my ($start, $input_file, $queue, $chr_names) = @_ ;

  my @chr_names = @$chr_names;

  my ($ssaha_command, $count);

  foreach my $chr_name (@chr_names) {
    my @chrs = split /\//, $chr_name;
    my $chr = $chrs[-1]; 

    my $target_file = "$chr_name"; ###need to check for exact file_name
    print "target_file is $target_file\n";
    $ssaha_command = "/usr/local/ensembl/bin/ssaha2 $input_file $target_file -align 0 -kmer 12 -seeds $seed -cut 100 -depth 5 -output vulgar";
    print "job_num is ", ++$count, " and start is $start and chr_name is $chr_name out is ssaha_out_$start\_$chr\n";
    my $call = "bsub -J $input_dir\_ssaha_job_$start $queue -e error_ssaha -o ssaha_out_$start\_$chr $ssaha_command";
    system ($call);
  }
}




sub parse_ssaha2_out {

  my ($input_file, $ssaha2_output) = @_;

  print "input_file is $input_file ssaha2_output is $ssaha2_output\n";

  open INPUT, "$input_file" or die "can't open $input_file : $!\n";
  open OUT, ">mapping_file" or die "can't open output_file : $!\n";

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
    my ($fseq,$lseq) = split /[WACTG\-\_]+/,$rec_seq{$name};
    my $snp_posf = length($fseq)+1;
    my $snp_posl = length($lseq)+1;
    my $tot_length = $snp_posf + $snp_posl;
    $snp_pos{$name} = "$snp_posf\_$snp_posl";
    $input_length{$name} = $tot_length;
  }


  open SSAHA, "$ssaha2_output" or die "can't open $output_file : $!\n";

  while (<SSAHA>) {
    #print ;
    if (/^(\S+)\s+(\d+)\s+(\d+)\s+(\+|\-)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\+|\-)\s+(\d+)\s+(.*)$/) {
      #next unless $1 eq "rs16822290";
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
    if ($h[0]->{'score'} > $h[1]->{'score'}) {
      find_results($h[0]);
    }
    elsif ($h[1]->{'score'} > $h[2]->{'score'}) {
      find_results($h[0],$h[1]);
    }
    elsif ($h[2]->{'score'} > $h[3]->{'score'}) {
      find_results($h[0],$h[1],$h[2]);
    }
    else {
      print "\tmore than 3 hits having same score for $q_id\n";
    }
  }
}

sub find_results {

  my ($h1,$h2,$h3) = @_;

  foreach my $h ($h1,$h2,$h3) {
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

    $q_start = $q_start+1;
    $t_start = $t_start+1;
      
    my @match_components = split /\s+/, $match;
      
    my ($new_q_start, $new_q_end, $new_t_start, $new_t_end);
      
    #print "$q_id,$q_start,$q_end,$q_strand,$t_id,$t_start,$t_end,$t_strand,$score,@match_components\n";
      
    while (@match_components) {
	
      my $type = shift @match_components;
      my $query_match_length = shift @match_components;
      my $target_match_length = shift @match_components;
      
      $new_t_start = $t_start + $target_match_length -1;
      my ($snp_q_pos, $snp_t_start, $snp_t_end, $snp_strand, $pre_target_match_length);
      my ($snp_posf,$snp_posl) = split /\_/, $snp_pos{$q_id};
      if ($q_strand eq '-') {
	$snp_q_pos = $snp_posl;
	$snp_strand = -1;
      }
      else {
	$snp_q_pos = $snp_posf;
	$snp_strand = 1;
      }
      
      $new_q_start = $q_start + $query_match_length - 1 if ($type eq 'M');
      $new_q_start = $q_start + $query_match_length + 1 if ($type eq 'G');
      $new_t_start = $t_start + $target_match_length -1 if ($type eq 'M');
      $new_t_start = $t_start + $target_match_length +1 if ($type eq 'G');
      
      #print "$snp_q_pos,$q_start,$new_q_start and $t_start\n";
      if ($snp_q_pos >= $q_start and $snp_q_pos <= $new_q_start) {
	if ($type eq 'M') {
	  $snp_t_start = $t_start + ($snp_q_pos-$q_start);
	  $snp_t_end = $snp_t_start;
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
      
      if ($snp_t_start && $snp_t_end && $snp_strand ) {
	my $final_score = $score/$input_length{$q_id};
	print OUT "$q_id\t$t_id\t$snp_t_start\t$snp_t_end\t$snp_strand\t$final_score\n";
	last;
      }
      
      $q_start = $new_q_start;
      $t_start = $new_t_start;
    }
  }
}

	
	
