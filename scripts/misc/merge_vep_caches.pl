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
use File::Copy;
use Storable qw(nstore_fd fd_retrieve);
use Bio::EnsEMBL::Registry;

my $config = {};

GetOptions(
  $config,
  'help|h',                  # displays help message
  
  'dir|d=s',                 
  'version|v=s',
  'species|s=s',
) or die "ERROR: Failed to parse command-line flags\n";

do { usage(); exit(0); } if $config->{help};

$config->{dir}     ||= $ENV{HOME}.'/.vep/';
$config->{version} ||= Bio::EnsEMBL::Registry->software_version;
$config->{species} ||= 'all';
$config->{compress}||= 'zcat';

$config->{cache_region_size} = 1e6;

my %match_species = ();
if($config->{species} ne 'all') {
  $match_species{$_} = 1 for split(/\,/, $config->{species});
}

# get all species with a refseq cache
opendir DIR, $config->{dir};
my @refseq_species = grep {
  !/^\./ &&
  -d $config->{dir}.'/'.$_ &&
#  -d $config->{dir}.'/'.$_.'/'.$config->{version} &&
  -d $config->{dir}.'/'.$_.'_refseq' # &&
#  -d $config->{dir}.'/'.$_.'_refseq/'.$config->{version}
} readdir DIR;
closedir DIR;

# get matched assemblies
my %assemblies = ();

foreach my $species(@refseq_species) {
  opendir ENS, $config->{dir}.'/'.$species;
  opendir REF, $config->{dir}.'/'.$species.'_refseq';

  my $v = $config->{version};
  my %ens_ass_match = map {$_ => 1} grep {$_ =~ /^$v/} readdir ENS;
  my %ref_ass_match = map {$_ => 1} grep {$_ =~ /^$v/} readdir REF;

  foreach my $ass(keys %ens_ass_match) {
    if(defined($ref_ass_match{$ass})) {
      $ass =~ s/^$v\_?//;
      push @{$assemblies{$species}}, $ass;
    }
  }
}

foreach my $species(@refseq_species) {
  next unless $config->{species} eq 'all' || defined($match_species{$species});
  
  print ">>>> Processing $species\n";

  foreach my $ass(@{$assemblies{$species}}) {

    print ">>>> Processing assembly $ass\n" if $ass;
  
    my $ens_root = join("/", $config->{dir}, $species, $config->{version}.($ass ? '_'.$ass : ''));
    my $ref_root = join("/", $config->{dir}, $species.'_refseq', $config->{version}.($ass ? '_'.$ass : ''));
    my $mrg_root = join("/", $config->{dir}, $species.'_merged', $config->{version}.($ass ? '_'.$ass : ''));
  
    # create merged dir
    mkdir($config->{dir}.'/'.$species.'_merged/');
    mkdir($mrg_root);
  
    # copy info
    copy($ens_root.'/info.txt', $mrg_root.'/info.txt');
  
    opendir ENSROOT, $ens_root;
    opendir REFROOT, $ref_root;
  
    # list chromosomes
    my %ens_chrs = map {$_ => 1} grep {-d $ens_root.'/'.$_ && !/^\./} readdir ENSROOT;
    my %ref_chrs = map {$_ => 1} grep {-d $ref_root.'/'.$_ && !/^\./} readdir REFROOT;
  
    closedir ENSROOT;
    closedir REFROOT;
  
    # mkdirs
    my %mrg_chrs = map {$_ => 1} (keys %ens_chrs, keys %ref_chrs);
  
    foreach my $chr(keys %mrg_chrs) {
      print " >>> Processing chromosome $chr\n";
    
      mkdir($mrg_root.'/'.$chr);
    
      # exists in both
      if(-d $ens_root.'/'.$chr && -d $ref_root.'/'.$chr) {
        opendir ENSCHR, $ens_root.'/'.$chr;
        opendir REFCHR, $ref_root.'/'.$chr;
      
        my %ens_files = map {$_ => 1} grep {!/^\./} readdir ENSCHR;
        my %ref_files = map {$_ => 1} grep {!/^\./ && !/var/ && !/reg/} readdir REFCHR;
      
        closedir ENSCHR;
        closedir REFCHR;
      
        # simply copy all var and reg files since they will be the same
        print "  >> Copying var and reg files\n";
      
        foreach my $file(grep {/var|reg/} keys %ens_files) {
          copy($ens_root.'/'.$chr.'/'.$file, $mrg_root.'/'.$chr.'/'.$file);
          delete $ens_files{$file};
        }
      
        print "  >> Merging transcript files\n";
      
        my %mrg_files = map {$_ => 1} (keys %ens_files, keys %ref_files);
      
        foreach my $file(keys %mrg_files) {
        
          # exists in both, need to concatenate
          if(-e $ens_root.'/'.$chr.'/'.$file && -e $ref_root.'/'.$chr.'/'.$file) {
            print "   > Merging $chr $file\n";
          
            # read in Ensembl cache
            open my $ens_fh, $config->{compress}." ".$ens_root.'/'.$chr.'/'.$file." |";
            my $ens_cache;
            $ens_cache = fd_retrieve($ens_fh);
            close $ens_fh;
          
            # add a flag to each transcript indicating which cache it came from
            $_->{_source_cache} = 'Ensembl' for @{$ens_cache->{$chr}};
          
            # do same for RefSeq
            open my $ref_fh, $config->{compress}." ".$ref_root.'/'.$chr.'/'.$file." |";
            my $ref_cache;
            $ref_cache = fd_retrieve($ref_fh);
            close $ref_fh;
            $_->{_source_cache} = 'RefSeq' for @{$ref_cache->{$chr}};
          
            # merge and sort transcript lists
            my $mrg_cache;
            @{$mrg_cache->{$chr}} = sort {$a->{start} <=> $b->{start}} (@{$ens_cache->{$chr}}, @{$ref_cache->{$chr}});
          
            # dump to new file
            open my $mrg_fh, "| gzip -9 -c > ".$mrg_root.'/'.$chr.'/'.$file or die "ERROR: Could not write to dump file";
            nstore_fd($mrg_cache, $mrg_fh);
            close $mrg_fh;
          }
        
          # otherwise simply copy
          else {
            my $root = -e $ens_root.'/'.$chr.'/'.$file ? $ens_root : $ref_root;
            copy($root.'/'.$chr.'/'.$file, $mrg_root.'/'.$chr.'/'.$file);
          }
        }
      }
    
      # only exists in one, simply copy all files
      else {
        print "  >> Copying all files\n";
      
        my $root = -d $ens_root.'/'.$chr ? $ens_root : $ref_root;
      
        opendir CHR, $root.'/'.$chr;
        copy($root.'/'.$chr.'/'.$_, $mrg_root.'/'.$chr.'/'.$_) for grep {!/^\./} readdir CHR;
        closedir CHR;
      }
    }
  }
}

print "Finished\n";

sub usage {
  
  print qq{
Merge Ensembl- and RefSeq-based VEP caches

Flags: (defaults follow in parentheses)

--dir | -d        Root directory of VEP caches (\$HOME/.vep/)
--version | -v    Cache version (current API version)
--species | -s    Species (all species found with both cache types)
};

  return 0;
}
