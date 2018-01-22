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

#perl fastq.pl [fasta_dir] [qual_dir]


use strict;

my ($i, $fasta, $qual, $totalfiles, $totalbases);
my ($ifname, $dname, $fin_unfin, $fname);

my ($fasta_dir, $qual_dir) = @ARGV;

while (<$fasta_dir/*fasta.gz>) {
  my $file = $_;
  $file =~ s/$fasta_dir\///;
  my ($fname) = $file =~ /^(.*)\.gz/; 
  print "fname is $fname\n";
  next if (-e "$fname\.fastq");

  open OUTQ, "> $fasta_dir/$fname\.fastq" or die "cannot open output file: $!";
  process_file($fasta_dir,$qual_dir,$fname);
  close OUTQ;
}
if($totalfiles and $totalfiles > 0) {
	print "total clones $totalfiles for $totalbases bases\n";
}

sub process_file {
    my ($fasta_dir,$qual_dir,$fname) = @_;
	if(-e "$fasta_dir/$fname") {
		open INA,"$fasta_dir/$fname" or return;
	} elsif (-e "$fasta_dir/$fname.Z") {
		open INA,"gunzip -c $fasta_dir/$fname.Z |" or return;
	      } elsif (-e "$fasta_dir/$fname\.gz") {
		open INA,"gunzip -c $fasta_dir/$fname\.gz |" or return;
	}
    $fasta = <INA>;
	if(-e "$qual_dir/$fname.qual") {
		open INQ,"$qual_dir/$fname.qual" or return;
	} elsif (-e "$qual_dir/$fname.qual.Z") {
		open INQ,"gunzip -c $qual_dir/$fname.qual.Z |" or return;
	} elsif (-e "$qual_dir/$fname\.qual.gz") {
		open INQ,"gunzip -c $qual_dir/$fname\.qual.gz |" or return;
	}
    while (defined($qual = <INQ>)) {
          next if $qual =~ /^\s+$/;
	  if ($qual =~ /^>(\S+)/) {
	    my ($header) = $1;
	    if (!($fasta =~ /$header/)) {
		print "out of sync $header\n";
		last;
	    }
	    ##change (f,r) to (p,q)
	    $header =~ s/\_A$/\.p1/;
	    $header =~ s/\_Z$/\.q1/;
	    print OUTQ "\@$header\n";
	    while (defined($fasta = <INA>)) {
	        next if $fasta =~ /^\s+$/;
		last if $fasta =~ /^>/;
		$totalbases += length($fasta)-1;
		print OUTQ $fasta;
	    }
	    print OUTQ "\+$header\n";
	    next;
	  }
	  $qual =~ s/NA/0/g;
	  print OUTQ $qual;
    }
    close INA;
    close INQ;
    $totalfiles++;
}
