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
use FileHandle;

## configure
my $config = {
  keep_headers  => 'CHROM,fileformat',
  delete_fields => 'QUAL,FILTER,INFO'
};

GetOptions(
  $config,
  'help|h',
  'input_file|i=s',
  'output_file|o=s',
  'keep_headers|k=s',
  'delete_fields|d=s'
);

if($config->{help}) {
  usage();
  exit(0);
}

$config->{$_} = [split(/\,/, ($config->{$_} || ''))] for qw(keep_headers delete_fields);

## input
my $in = *STDIN;

# file provided?
if($config->{input_file}) {
  my $file = $config->{input_file};
  
  die("ERROR: Input file not found\n") unless -e $file;
  
  $in = FileHandle->new();
  
  # gzipped?
  if($file =~ /\.gz$/ || -B $file) {
    $in->open('gzip -dc '.$file.' | ') or die("ERROR: Could not read from input file ".$file."\n");
  }
  
  # uncompressed
  else {
    die("ERROR: Input file looks like a binary (is it gzipped but not named [file].gz?)\n") if -B $file;
    $in->open($file) or die("ERROR: Could not read from input file ".$file."\n");
  }
}

## output
my $out = *STDOUT;

if($config->{output_file}) {
  my $file = $config->{output_file};
  $out = FileHandle->new();
  
  if($file =~ /\.gz/) {
    $out->open('| bgzip -c > '.$file) or die("ERROR: Could not write to output file ".$file."\n");
  }
  else {
    $out->open('> '.$file) or die("ERROR: Could not write to output file ".$file."\n");
  }
}

my $meta = {};
my $headers = {};

while(<$in>) {
  
  # headers
  if(/^\#/) {
      
    # parse column headers
    if(/^\#CHROM/) {
      my $l = $_;
      chomp($l);
      my @data = split(/\t/, $l);
    
      $headers->{$data[$_]} = $_ for (0..$#data);
    }
    
    if(scalar @{$config->{keep_headers}} && $config->{keep_headers}->[0] eq 'ALL') {
      print $out $_;
      next;
    }
    
    foreach my $inc(@{$config->{keep_headers}}) {
      
      if(/^\#+$inc/) {
        print $out $_;
        next;
      }
    }
  }
  
  # data
  else {
    die("ERROR: Failed to parse headers\n") unless scalar keys %$headers;
    
    chomp;
    my @data = split /\t/, $_;
    
    # has genotypes?
    if(!exists($meta->{has_genotypes})) {
      
      # default to no
      $meta->{has_genotypes} = 0;
      
      if($headers->{FORMAT}) {
        my $part_num = 0;
        
        foreach my $part(split(/\:/, $data[$headers->{FORMAT}] || '')) {
          $meta->{has_genotypes} = ++$part_num if $part eq 'GT';
        }
      }
    }
    
    # delete/null fields 
    $data[$headers->{$_}] = '.' for @{$config->{delete_fields}};
    
    # remove non-genotype bits
    if($meta->{has_genotypes}) {
      
      $data[$headers->{FORMAT}] = 'GT';
      
      my $from = $headers->{FORMAT} + 1;
      my $gt_part = $meta->{has_genotypes} - 1;
    
      for my $i($from..$#data) {
        my @parts = split(/\:/, $data[$i]);
        $data[$i] = $parts[$gt_part];
      }
    }
    
    print $out join("\t", @data)."\n";
  }
}

sub usage {
    my $usage =<<END;

Script to strip down a VCF file's size, for use with VCFCollection API.

NB: Stripped VCFs should NOT be shared except for use with the Ensembl Variation API!

Usage:
perl $0 [arguments]

Arguments
=============

-h | --help            Display this message and quit

-i | --input_file      Input file (reads from STDIN if not specified)
                       File may be uncompressed, gzipped or bgzipped.
                       
-o | --output_file     Output file (writes to STDOUT if not specified)
                       If file name ends in ".gz", output will be compressed with bgzip
                       
-k | --keep_headers    Comma separated list of header line keys to keep
                       By default retains the fileformat and CHROM header lines
                       Set to "ALL" to retain all header lines.
                       
-d | --delete_fields   Comma separated list of fields/columns to "blank" i.e. set to "."
                       By default deletes QUAL, FILTER and INFO.
                       NB: future API versions may use information from the INFO field.

END

    print $usage;
}