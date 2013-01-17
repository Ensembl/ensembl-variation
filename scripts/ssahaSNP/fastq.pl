#!/usr/bin/env perl

=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/legal/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

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
