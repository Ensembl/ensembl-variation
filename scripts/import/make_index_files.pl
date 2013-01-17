#!/usr/bin/env perl

use strict;
use warnings;

use Bio::Index::Fastq;

my $index_filename = shift(@ARGV);

my $inx = Bio::Index::Fastq->new(-filename => $index_filename,
				 -write_flag => 1);

$inx->make_index(@ARGV);

exit 0;
