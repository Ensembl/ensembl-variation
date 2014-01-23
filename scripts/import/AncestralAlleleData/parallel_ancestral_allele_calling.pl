# Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

use Fcntl qw(:flock);
use File::Path qw(make_path);
use File::Copy;
use FindBin qw($Bin);
use Getopt::Long;
use Pod::Usage;
use POSIX;

my ($TMP_DIR, $PERLBIN, $fasta_dir, $sub_processes, $help);

GetOptions(
	'TMP_DIR=s'       => \$TMP_DIR,
	'fasta_dir=s'     => \$fasta_dir,
	'sub_processes=s' => \$sub_processes,
	'help|h'          => \$help,
) or die pod2usage(1);
pod2usage(1) if $help;

die "Argument list is not complete, try --help for usage" unless ($TMP_DIR && $fasta_dir && $sub_processes);

$PERLBIN = `which perl`;
my $script = $Bin . '/ancestral_allele_calling.pl';
my $sub_process = 1;
my $bsub_queue_name = 'normal';
my $status_file = "$TMP_DIR/monitor_parallel_processes.txt";

open STATUS, ">>$status_file";
flock(STATUS, LOCK_EX);
while ($sub_process <= $sub_processes) {
	my $pid = fork();
	if ($pid) {
		waitpid($pid, 0);
	} elsif ($pid == 0) {
		my @call_params = (
			"bsub -q $bsub_queue_name",
			"-J test_call_$sub_process",
			"-o $TMP_DIR/test_call_$sub_process.out",
			"-e $TMP_DIR/test_call_$sub_process.err",
			"-M2500000 -R\"select[mem>2500] rusage[mem=2500]\"",
			"$PERLBIN $script",
			"-fasta_files_dir $fasta_dir/$sub_process/",
			"-variation_feature $TMP_DIR/variation_feature.$sub_process.txt",
			"-save_results $TMP_DIR/ancestral_allele_output.$sub_process.txt",
			"-status_file $status_file"
		);
		my $call = join(" ", @call_params);
		system($call);
		exit(0);
	} else {
		die "Couldn't fork: $!\n";
	}
	$sub_process++;
}

print STATUS "finished at " . localtime() . "\n";
flock(STATUS, LOCK_UN);
close(STATUS);

__END__

=head1 NAME

parallel_ancestral_allele_calling.pl

=head1 DESCRIPTION

# steps:
# run different processes
# monitor when they have finished

=head1 SYNOPSIS

parallel_ancestral_allele_calling.pl [arguments]

=head1 ARGUMENTS

=over 4

=item B<--TMP_DIR DIR>

Directory used by pre-processing script to save intermediate results.

=item B<--fasta_dir DIR>

Directory containing sub-directories as divided by pre-processing script.

=item B<--sub_processes NUMBER>

Number of subprocesses computed by load balancing scheme.

=item B<--help>

Show help

=head1

For help with this script address questions to http://lists.ensembl.org/mailman/listinfo/dev

