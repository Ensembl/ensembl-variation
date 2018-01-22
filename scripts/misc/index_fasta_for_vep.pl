#!/usr/bin/env perl

=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

use strict;
use Getopt::Long;
use Bio::DB::Fasta;
use File::Copy;

my $config = {};

GetOptions(
  $config,
  'help|h',                  # displays help message
  
  'dir|d=s',                 
  'fasta|f=s',
  'tmpdir|t=s',
  'species|s=s',
  'version|v=i',
  'target|g=s',
) or die "ERROR: Failed to parse command-line flags\n";

do { usage(); exit(0); } if $config->{help};

$config->{dir}     ||= $ENV{HOME}.'/.vep/';
$config->{tmpdir}  ||= '/tmp';

die("ERROR: no FASTA directory specified (--fasta)\n") unless defined($config->{fasta});

my %match_species = ();
if($config->{species} ne 'all') {
  $match_species{lc($_)} = 1 for split(/\,/, $config->{species});
}

# get all species cache dirs
opendir DIR, $config->{dir};
my @caches = grep {
  !/^\./ &&
  -d $config->{dir}.'/'.$_
} readdir DIR;
closedir DIR;

foreach my $cache(@caches) {
  my $species = $cache;
  $species =~ s/\_merged|\_refseq//;
  my ($fastadir, $tmpdir) = ($config->{fasta}, $config->{tmpdir});
  
  if(scalar keys %match_species) {
    next unless $match_species{lc($species)};
  }
  
  print STDERR ">> Processing $cache\n";
  
  if(!opendir DIR, "$fastadir/$species/dna") {
    warn "WARNING: Could not read from directory $fastadir/$species/dna\n$@\n";
    next;
  }
  my @files = grep {$_ !~ /^\./} readdir DIR;
  closedir DIR;

  # remove repeat/soft-masked files
  @files = grep {!/_(s|r)m\./} @files;

  my ($file) = grep {/primary_assembly.fa.gz$/} @files;
  ($file) = grep {/toplevel.fa.gz$/} @files if !defined($file);

  unless(defined($file)) {
    warn "WARNING: No FASTA file found for $cache\n";
    next;
  }
  
  my $unz = $file;
  $unz =~ s/\.gz$//;
  
  print STDERR " > Unpacking $fastadir/$species/dna/$file\n";
  my $out = `zcat $fastadir/$species/dna/$file > $tmpdir/$unz 2>&1`;
  die("ERROR: failed to unpack FASTA file\n$out") if $out;
  
  print STDERR " > Indexing $tmpdir/$unz\n";
  Bio::DB::Fasta->new("$tmpdir/$unz");
  
  if(!-e "$tmpdir/$unz\.index") {
    warn "WARNING: Failed to create index or index file not found\n";
    next;
  }
  
  # work out where to put the index file
  my $tdir;
  
  if(!defined($config->{target})) {
    my $cdir = $config->{dir};
    opendir DIR, $cdir.'/'.$cache;
    my @dirs = grep {
      !/^\./ && -d $cdir.'/'.$cache.'/'.$_
    } readdir DIR;
    closedir DIR;
  
    my %vs = map {(split('_', $dirs[$_]))[0] => $dirs[$_]} (0..$#dirs);
    $tdir = $vs{(sort {$a <=> $b} keys %vs)[-1]};
    if(defined($config->{version})) {
      $tdir = $vs{$config->{version}};
    }
    unless(defined($tdir)) {
      warn "WARNING: No matching cache directory (version ".($config->{version} || '?').") found in $cache";
      next;
    }
    
    $tdir = "$cdir/$cache/$tdir";
  }
  else {
    $tdir = $config->{target};
  }
  
  print STDERR " > Moving $tmpdir/$unz\.index to $tdir/$unz\.index\n";
  move("$tmpdir/$unz\.index", "$tdir/$unz\.index");
  
  unlink("$tmpdir/$unz");
}

sub usage {
  
  print qq{
Create FASTA index files for VEP caches

Flags: (defaults follow in parentheses)

--dir | -d        Root directory of VEP caches (\$HOME/.vep/)
--fasta | -f      Directory containing FASTA file dirs
--tmpdir | -t     Temporary dir - used to unpack FASTA gzip files (/tmp)
--target | -g     Target dir for index files (per-species cache dir)
--species | -s    Species (all species found)
--version | -v    Version (latest)
};

  return 0;
}
