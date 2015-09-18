#!/usr/bin/env perl

use strict;
use warnings;

my @files = @ARGV;

print qq{
<html>
  <title>VEP documentation</title>
  <body>
};

foreach my $file(@files) {
  open IN, $file or die("ERROR: Could not read from file $file\n");

  my $in_body = 0;

  while(<IN>) {
    if(/\<title\>/) {
      s/title/h1"/g;
      print;
    }

    elsif(/\<\/?body\>/) {
      $in_body = 1 - $in_body;
    }

    elsif($in_body) {
      print;
    }
  }

  close IN;

  print '<div style="clear:both"><h1 style="page-break-after: always">&nbsp</h1></div>';
}

print "</body></html>\n";
