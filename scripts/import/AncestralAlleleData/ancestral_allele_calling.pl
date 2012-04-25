use strict;
use warnings;

use AncestralAlleleCaller;
use Fcntl qw(:flock);
use FindBin qw($Bin);
use Getopt::Long;
use lib($Bin);
my ($fasta_files_dir, $input, $output, $status_file);

GetOptions(
	'fasta_files_dir=s' => \$fasta_files_dir,
	'variation_feature=s' => \$input,
	'save_results=s' => \$output,
	'status_file=s' => \$status_file,
);

my $fasta_db = AncestralAlleleCaller->new($fasta_files_dir); # force to build new index file

open INPUT, "<$input" or die "Can't open $input";
open OUTPUT, ">$output" or die "Can't open $output";
while (<INPUT>) {
	chomp;
	my ($seq_region_name, $seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand, $variation_id) = split /\t/;
	my $ancestral_allele = uc $fasta_db->get_ancestral_allele($seq_region_name, $seq_region_start, $seq_region_end);
	# checks
	# Is insertion - no ancestral allele for an insertion
	if ($seq_region_start > $seq_region_end) {
		$ancestral_allele = '-';
	}
	# Is deletion of size > 50 - no display for large deletions
	if ($seq_region_end - $seq_region_start + 1 > 50) {
		$ancestral_allele = '-';
	}
	# Is on minus strand
	if ($seq_region_strand == -1) {
		if ($ancestral_allele =~ /(A-Z)*/) {
			reverse_comp(\$ancestral_allele);
		} else {
			$ancestral_allele = '-';
		}
	}
	print OUTPUT "$variation_id\t$seq_region_name\t$seq_region_start\t$seq_region_end\t$seq_region_strand\t$ancestral_allele\n";
}
close INPUT;
close OUTPUT;

open STATUS, ">>$status_file";
flock(STATUS, LOCK_EX);
print STATUS $fasta_files_dir, "\n";
print STATUS $input, "\n";
print STATUS $output, "\n";
print STATUS "finished at " . localtime() . "\n";
flock(STATUS, LOCK_UN);
close(STATUS);

1;

