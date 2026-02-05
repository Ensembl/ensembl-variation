=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2026] EMBL-European Bioinformatics Institute

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

=head1 NAME

Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin

=head1 SYNOPSIS

  package FunkyPlugin;

  use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

  sub new {
    
  }
  
  sub feature_types {
    
  }

  sub get_header_info {
    return {
      # Same name as in parse_data() result hash
      PLUGIN_SCORE => "Description of funky plugin"
    };
  }

  sub run {
    my ($self, $transcript_variation_allele) = @_;

    my @data    = @{$self->get_data(...)} # get data from parse_data()
    my $results = ...                     # do analysis
    
    return {
      FUNKY_PLUGIN => $results
    };
  }
  
  sub parse_data {
    my ($self, $line) = @_;
    # Parse every non-header line based on a delimiter, i.e for `\t`,
    # `CHROM   POS   REF   ALT   SCORE`
    # `1   1234567   A   T   0.95` <<<--- starts here
    my ($chrom, $start, $ref, $alt, $score) = split /\t/, $line;
    
    ... # parse data
    
    # example of hash returned with parsed data
    return {
      ref => $ref,
      alt => $alt,
      start => $start,
      result => {
        # Same name as in get_header_info() description
        PLUGIN_SCORE => $score
      }
    };
  }

  1;

=head1 DESCRIPTION

To make writing VEP plugin modules that load data from Tabix
files easier, get your plugin to inherit from this class,
override (at least) the feature_types, get_header_info, run
and parse_data methods to behave according to the
documentation below, and then run the VEP with your plugin
using the --plugin <module name> command line option.

=cut

package Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);
use Digest::MD5 qw(md5_hex);

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

use strict;
use warnings;

my $DEFAULT_EXPAND_LEFT = 0;
my $DEFAULT_EXPAND_RIGHT = 1e6;
my $DEFAULT_CACHE_SIZE = 1;

my $DEBUG = 0;

our ($CAN_USE_HTS_PM, $CAN_USE_TABIX_PM);

BEGIN {
  if (eval q{ require Bio::DB::HTS::Tabix; 1 }) {
    $CAN_USE_HTS_PM = 1;
  }
  if (eval q{ require Tabix; 1 }) {
    $CAN_USE_TABIX_PM = 1;
  }
  unless($CAN_USE_TABIX_PM || $CAN_USE_HTS_PM) {
    # test tabix
    die "ERROR: tabix does not seem to be in your path\n" unless `which tabix 2>&1` =~ /tabix$/;
  }
}

=head2 new

  Arg [1]    : a VEP configuration hashref
  Arg [>1]   : any parameters passed on the VEP command line, will be
               stored as a listref in $self->{params}
  Description: Creates and returns a new instance of this plugin based on
               the constructor from
               Bio::EnsEMBL::Variation::Utils::BaseVepPlugin
  Returntype : Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin
               instance (most likely a user-defined subclass)
  Status     : Experimental
  
=cut

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  return $self;
}

=head2 get_user_params

  Description: Sets user-defined parameters for expand_left, expand_right
               and cache_size in the plugin instance
  Status     : Experimental
  
=cut

sub get_user_params {
  my $self = shift;

  my $params = $self->params_to_hash;

  foreach my $param_name(qw(expand_left expand_right cache_size)) {
    $self->$param_name($params->{$param_name}) if exists($params->{$param_name});
  }
}

=head2 get_header_info

  Description: Return a hashref with any Extra column keys as keys and a
               description of the data as a value, this will be included
               in the VEP output file header to help explain the results
               of this plugin. Plugins that do not want to include
               anything in the header should return undef.
  Returntype : hashref or undef
  Status     : Experimental

=cut

=head2 files

  Description: Return valid Tabix file paths from plugin parameters
               
               Plugin parameters are potential files if:
                 1. Their keys start with "file" and digits (e.g. file1).
                    The user would pass the parameters like so:
                      "--plugin MyPlugin,file1=ex1,file2=ex2".
                 2. Parameters end in ".gz" or if their filenames exist
                    (this step is skipped if step 1 already found
                    matching parameters).
               
               Finally, check_file() is used to confirm the existence of
               Tabix files and their indexes in the same path. Raises
               error if any potential file does not exist.
  Returntype : arrayref
  Status     : Experimental
  
=cut

sub files {
  my $self = shift;

  if(!exists($self->{_files})) {
    my @files;

    # from hashref of plugin parameters, get files whose keys start with "file"
    # followed by one or more digits
    my $params = $self->params_to_hash;
    foreach my $key(grep {/^file\d+/} keys %$params) {
      push @files, $params->{$key};
    }

    # if no files match the previous criteria, get files from listref of plugin
    # parameters that either end in ".gz" or that are valid filenames
    unless(@files) {
      foreach my $p(@{$self->params}) {
        push @files, $p if $p =~ /\.gz$/ || -e $p;
      }
    }

    $self->check_file($_) for @files;

    $self->{_files} = \@files;
  }

  return $self->{_files};
}

=head2 add_file

  Arg[1]     : file path
  Description: Validate and append Tabix file path to $self->{_files}
  
               check_file() is used to confirm the existence of Tabix
               files and their indexes in the same path. Raises error if
               any tested file does not exist.
  Returntype : boolean
  Status     : Experimental
  
=cut

sub add_file {
  my ($self, $file) = @_;

  $self->check_file($file);

  push @{$self->{_files}}, $file;
}

=head2 check_file

  Arg[1]     : file path
  Description: Validate the existence of Tabix files and their indexes in
               the same path. Raises error if any tested file does not
               exist.
  Returntype : boolean
  Status     : Experimental
  
=cut

sub check_file {
  my ($self, $file) = @_;

  die("ERROR: No file specified\n") unless $file;

  if($file !~ /tp\:\/\//) {
    die "ERROR: Data file $file not found\n" unless -e $file;
    die "ERROR: Tabix index file $file\.tbi not found - perhaps you need to create it first?\n" unless -e $file.'.tbi';
  }

  return 1;
}

=head2 expand_left

  Arg[1]     : positions to expand to the left (optional)
  Description: Store expand left position in $self->{_expand_left}
               Without arguments, value is set to $DEFAULT_EXPAND_LEFT
  Returntype : scalar
  Status     : Experimental
  
=cut

sub expand_left {
  my $self = shift;

  $self->{_expand_left} = shift if @_;
  $self->{_expand_left} = $DEFAULT_EXPAND_LEFT if !exists($self->{_expand_left});

  return $self->{_expand_left};
}

=head2 expand_right

  Arg[1]     : positions to expand to the right (optional)
  Description: Store expand right position in $self->{_expand_right}
               Without arguments, value is set to $DEFAULT_EXPAND_RIGHT
  Returntype : scalar
  Status     : Experimental
  
=cut

sub expand_right {
  my $self = shift;

  $self->{_expand_right} = shift if @_;
  $self->{_expand_right} = $DEFAULT_EXPAND_RIGHT if !exists($self->{_expand_right});

  return $self->{_expand_right};
}

=head2 cache_size

  Arg[1]     : cache size (optional)
  Description: Get cache size and store it in $self->{_cache_size}
               In case of no argument, cache size is set to
               $DEFAULT_CACHE_SIZE
  Returntype : scalar
  Status     : Experimental
  
=cut

sub cache_size {
  my $self = shift;

  $self->{_cache_size} = shift if @_;
  $self->{_cache_size} = $DEFAULT_CACHE_SIZE if !exists($self->{_cache_size});

  return $self->{_cache_size}; 
}

=head2 get_data

  Arg[1]     : chromosome
  Arg[2]     : start
  Arg[3]     : end
  Arg[4]     : file (if not defined, retrieves data from all $self->files)
  Description: Get data from files specified in $self->files for a given
               genomic region
  Returntype : hashref
  Status     : Experimental
  
=cut

sub get_data {
  my ($self, $c, $s, $e, $f) = @_;

  die("ERROR: No chromosome specified\n") unless $c;
  die("ERROR: No start specified\n") unless $s;
  die("ERROR: No end specified\n") unless $e;

  # we use two levels of caching
  # 1) results cache for specific coords
  # 2) region cache to reduce multiple lookups for sequential "close" coords

  my $pos_string = join("_", $c, $s, $e);

  my $cache = $self->cache($f);

  # check results cache first
  if(exists($cache->{results}) && exists($cache->{results}->{$pos_string})) {
    print STDERR "Using results cache\n" if $DEBUG;
    return $cache->{results}->{$pos_string};
  }

  # now check the region cache
  my (@result, @missing_regions, $hit_cache);
  my $regions_used = 0;

  # we only use the region cache if we are allowed to expand left or right
  my ($expand_left, $expand_right) = ($self->expand_left, $self->expand_right);
  my $use_data_cache = ($expand_left || $expand_right) ? 1 : 0;

  if($use_data_cache && $cache->{$c}) {
    my $regions = $cache->{$c}->{regions};

    # iterate through them backwards as most likely to hit last one pushed on to array first
    for(my $i=(scalar @$regions - 1); $i>=0; $i--) {
      my $region = $regions->[$i];
      if(overlap($s, $e, $region->[0], $region->[1])) {

        print STDERR "Using data cache\n" if $DEBUG;

        # flag that we've hit the cache and store how many regions we've used
        $hit_cache = 1;
        
        # for partial overlaps we store the bits we don't have in @missing_regions
        if($s < $region->[0]) {
          push @missing_regions, [$s, $region->[0] - 1];
        }
        if($e > $region->[1]) {
          push @missing_regions, [$region->[1] + 1, $e];
        }

        my $filtered = $self->_filter_by_pos($cache->{$c}->{data}->[$i], $s, $e);
        $regions_used++ if scalar @$filtered;

        push @result, @$filtered;
      }
    }
  }

  # if we hit the cache, we can assume that our original start-end has been covered
  # with any missing pieces added to @missing_regions
  # otherwise we didn't hit the cache, so we add the original start-end
  push @missing_regions, [$s, $e] unless $hit_cache;

  foreach my $region(@missing_regions) {

    my ($r_s, $r_e) = @$region;

    # expand?
    $r_s -= $expand_left;
    $r_e += $expand_right;

    my $tmp_data = $self->_get_data_uncached($c, $r_s, $r_e, $f);

    # cache the data
    $self->_add_data_to_cache($c, $r_s, $r_e, $f, $tmp_data) if $use_data_cache;

    # we don't need to filter it unless we're using the cache
    push @result, $use_data_cache ? @{$self->_filter_by_pos($tmp_data, $s, $e)} : @$tmp_data;
    
    $regions_used++;
  }

  # now unique it, but only if we queried more than one region
  my %seen = ();
  my @uniq_result = ();

  if($regions_used > 1) {
    foreach my $d(@result) {
      push @uniq_result, $d->[1] unless $seen{$d->[0]};
      $seen{$d->[0]} = 1;
    }
  }
  else {
    @uniq_result = map {$_->[1]} @result;
  }

  $self->_add_result_to_cache($f, $pos_string, \@uniq_result);
  
  return \@uniq_result;
}

=head2 cache

  Arg[1]     : File (optional; if not defined, conglomerate data from all files)
  Description: Get cache (empty if no cache available)
  Returntype : hash
  Status     : Experimental
  
=cut

sub cache {
  my ($self, $file) = @_;

  my $cache;
  if (defined $file) {
    $cache = $self->{_tabix_cache}->{$file} ||= {};
  } else {
    $cache = $self->{_tabix_cache} ||= {};
  }
  return $cache;
}

=head2 parse_data

  Args       : A string containing a line from input files
  Description: Data parsing logic (replace this method in a new plugin)
  
               This method runs before the run() method from
               Bio::EnsEMBL::Variation::Utils::BaseVepPlugin and allows
               to customise how data is parsed per line of a Tabix file.
  
               When implementing parse_data() in a new plugin, start by
               splitting the line by tabs (or another delimeter) into
               different variables, for instance:
                 
                 my ($self, $line) = @_;
                 my ($chr, $start, $ref, $alt, $score)=split /\t/, $line;
                 
               Manipulate parsed data as needed and return hash of column
               names and values, such as:
                 
                  return {
                    ref => $ref,
                    alt => $alt,
                    start => $start,
                    result => {
                      PLUGIN_SCORE => $score
                    }
                  };
               
               Parsed data will be accessible in method run() by using:
                 
                 my @data = @{$self->get_data(...)}
                 
  Returntype : hash
  Status     : Experimental
  
=cut

sub parse_data {
  my $self = shift;
  return shift;
}

=head2 identify_data

  Arg[1]     : Line
  Arg[2]     : Parsed line
  Description: Get MD5 digest of a line in hexadecimal format
  Returntype : scalar
  Status     : Experimental
  
=cut

sub identify_data {
  my ($self, $line, $parsed) = @_;
  return md5_hex($line);
}

sub _get_data_uncached {
  my $self = shift;

  if($CAN_USE_HTS_PM) {
    return $self->_get_data_hts(@_);
  }
  elsif($CAN_USE_TABIX_PM) {
    return $self->_get_data_pm(@_);
  }
  else {
    return $self->_get_data_cl(@_);
  }
}

sub _add_data_to_cache {
  my ($self, $c, $s, $e, $f, $data) = @_;

  my $cache = $self->cache($f)->{$c} ||= {};

  push @{$cache->{regions}}, [$s, $e];
  push @{$cache->{data}}, $data;

  # trim
  if(scalar @{$cache->{regions}} > $self->cache_size) {
    shift @{$cache->{regions}};
    shift @{$cache->{data}};
  } 
}

sub _add_result_to_cache {
  my ($self, $file, $pos_string, $result) = @_;

  my $cache = $self->cache($file);

  $cache->{results}->{$pos_string} = $result;
  push @{$cache->{results_order}}, $pos_string;

  if(scalar @{$cache->{results_order}} > $self->cache_size) {
    my $del = shift @{$cache->{results_order}};
    delete $cache->{results}->{$del};
  }
}

sub _filter_by_pos {
  my ($self, $list, $s, $e) = @_;

  my @result;

  foreach my $data(@$list) {
    my $start = $self->get_start($data->[1]);
    push @result, $data if overlap($start, $self->get_end($data->[1]), $s, $e);
    last if $start > $e;
  }

  return \@result;
}

sub _get_data_hts {
  my ($self, $c, $s, $e, $f) = @_;

  my @data;

  my @files = defined $f ? ($f) : @{$self->files};
  foreach my $file(@files) {
    my $hts_obj = $self->_hts_obj($file);
    my $valids = $self->{_valids}->{$file} ||= $hts_obj->seqnames;

    my $iter = $hts_obj->query(
      sprintf(
        "%s:%i-%i",
        $self->_get_source_chr_name($c, $file, $valids),
        $s, $e
      )
    );
    next unless $iter;

    while(my $line = $iter->next) {
      my $parsed = $self->parse_data($line,$file);

      push @data, [$self->identify_data($line), $parsed] if $parsed;
    }
  }

  return \@data;
}

sub _get_data_pm {
  my ($self, $c, $s, $e, $f) = @_;

  my @data;

  my @files = defined $f ? ($f) : @{$self->files};
  foreach my $file(@files) {
    my $tabix_obj = $self->_tabix_obj($file);
    my $valids = $self->{_valids}->{$file} ||= [$tabix_obj->getnames];

    my $iter = $tabix_obj->query($self->_get_source_chr_name($c, $file, $valids), $s, $e);
    next unless $iter && $iter->{_};

    while(my $line = $tabix_obj->read($iter)) {
      my $parsed = $self->parse_data($line,$file);

      push @data, [$self->identify_data($line), $parsed] if $parsed;
    }
  }

  return \@data;
}

sub _get_data_cl {
  my ($self, $c, $s, $e, $f) = @_;
  
  my @data;

  my @files = defined $f ? ($f) : @{$self->files};
  foreach my $file(@files) {
    my $valids = $self->{_valids}->{$file} ||= [split("\n", `tabix -l $file`)];

    open TABIX, sprintf(
      "tabix -f %s %s:%i-%i |",
      $file,
      $self->_get_source_chr_name($c, $file, $valids),
      $s, $e
    );

    while(<TABIX>) {
      chomp;
      s/\r$//g;

      my $parsed = $self->parse_data($_,$file);
      
      push @data, [$self->identify_data($_), $parsed] if $parsed;
    }

    close TABIX;
  }

  return \@data;
}

sub _hts_obj {
  my ($self, $file) = @_;
  return $self->{_tabix_obj}->{$file} ||= Bio::DB::HTS::Tabix->new(filename => $file);
}

sub _tabix_obj {
  my ($self, $file) = @_;
  return $self->{_tabix_obj}->{$file} ||= Tabix->new(-data => $file);
}

sub _get_source_chr_name {
  my ($self, $chr, $set, $valids) = @_;

  $set    ||= 'default';
  $valids ||= [];

  my $chr_name_map = $self->{_chr_name_map}->{$set} ||= {};

  if(!exists($chr_name_map->{$chr})) {
    my $mapped_name = $chr;

    @$valids = @{$self->can('valid_chromosomes') ? $self->valid_chromosomes : []} unless @$valids;
    my %valid = map {$_ => 1} @$valids;

    unless($valid{$chr}) {

      # still haven't got it
      if($mapped_name eq $chr) {

        # try adding/removing "chr"
        if($chr =~ /^chr/i) {
          my $tmp = $chr;
          $tmp =~ s/^chr//i;

          $mapped_name = $tmp if $valid{$tmp};
        }
        elsif($valid{'chr'.$chr}) {
          $mapped_name = 'chr'.$chr;
        }
      }
    }

    $chr_name_map->{$chr} = $mapped_name;
  }

  return $chr_name_map->{$chr};
}

=head2 FREEZE

  Description: Delete Tabix object from plugin instance
  Returntype : Return value of delete()
  Status     : Experimental
  
=cut

sub FREEZE {
  my $self = shift;
  delete $self->{_tabix_obj};
}

1;

