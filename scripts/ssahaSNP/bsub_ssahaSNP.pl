#! /software/bin/perl
#perl ./bsub_ssahaSNP.pl -reads_dir [reads_dir] -output_dir [out_dir] -target_file [target_dir]/chimp_repeat_masked.fa

use strict;
use Getopt::Long;

our ($reads_dir, $output_dir, $target_file, $ssahaSNP);

GetOptions('reads_dir=s'    => \$reads_dir,
	   'output_dir=s'    => \$output_dir,
	   'target_file=s'   => \$target_file,
	   'ssahaSNP=s'		=> \$ssahaSNP,
	  );

($reads_dir and $output_dir and $target_file) || die "We need reads_dir, output_dir and target_file defined\n";

my $queue = "hugemem -R 'select[mem>6000] rusage[mem=6000]'";

my @tar = split /\//, $target_file;
my $tar_file = $tar[-1];

opendir DIRENTRY, "$reads_dir" || die "Failed to open dir : $!";
my @files = grep(/perlegen-mouse-chip-111/,readdir(DIRENTRY));

#@files = ("perlegen-mouse-chip-1114058134.fastq");
@files = ("test4.fastq");
print "my file is #@files#\n";
my $count = 1;

foreach my $file (@files) {
  open IN, "$reads_dir/$file" or die "can't open in_file:$!";
  my @line = <IN>;
  my @nums = grep (/^\@/, @line);

  my $num = scalar @nums;

#  $num = 10000;
  print "num is $num\n";
  my $n;
  my $size = 20000;

  for ($n = 1;$n<=$num;$n+=$size) {
    my $end = $n+$size-1;
    $end = $num if ($end > $num);

    #my $call = "bsub -q $queue -e $output_dir/ssaha_out\_$file\_$n.err -o $output_dir/ssaha_out\_$file\_$n $ssahaSNP $reads_dir/$file $target_file";# -start $n -end $end -cut 100 -depth 10 -best 1";
    my $call = "bsub -q $queue -J 'ssahaSNP_$n' -e $output_dir/ssaha_out\_$file\_$n.err -o $output_dir/ssaha_out\_$file\_$n\_$tar_file $ssahaSNP $reads_dir/$file  $target_file -start $n -end $end -output cigar -NQS 1 -tags 1 -quality 23 -memory 2000";
    print "$call $count\n";
    system ($call);
    $count++;
    last;
  }
}
