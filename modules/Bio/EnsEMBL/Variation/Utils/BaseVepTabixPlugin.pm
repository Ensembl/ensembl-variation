=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

  sub run {
    my ($self, $transcript_variation_allele) = @_;

    my $results = ... # do analysis
    
    return {
      FUNKY_PLUGIN => $results
    };
  }

  1;

=head1 DESCRIPTION

To make writing plugin modules for the VEP easier, get 
your plugin to inherit from this class, override (at least)
the feature_types, get_header_info and run methods to behave
according to the documentation below, and then run the VEP
with your plugin using the --plugin <module name> command
line option.

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

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  return $self;
}

sub get_user_params {
  my $self = shift;

  my $params = $self->params_to_hash;

  foreach my $param_name(qw(expand_left expand_right cache_size)) {
    $self->$param_name($params->{$param_name}) if exists($params->{$param_name});
  }
}

sub files {
  my $self = shift;

  if(!exists($self->{_files})) {
    my @files;

    my $params = $self->params_to_hash;
    foreach my $key(grep {/^file\d+/} keys %$params) {
      push @files, $params->{$key};
    }

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

sub add_file {
  my ($self, $file) = @_;

  $self->check_file($file);

  push @{$self->{_files}}, $file;
}

sub check_file {
  my ($self, $file) = @_;

  die("ERROR: No file specified\n") unless $file;

  if($file !~ /tp\:\/\//) {
    die "ERROR: Data file $file not found\n" unless -e $file;
    die "ERROR: Tabix index file $file\.tbi not found - perhaps you need to create it first?\n" unless -e $file.'.tbi';
  }

  return 1;
}

sub expand_left {
  my $self = shift;

  $self->{_expand_left} = shift if @_;
  $self->{_expand_left} = $DEFAULT_EXPAND_LEFT if !exists($self->{_expand_left});

  return $self->{_expand_left};
}

sub expand_right {
  my $self = shift;

  $self->{_expand_right} = shift if @_;
  $self->{_expand_right} = $DEFAULT_EXPAND_RIGHT if !exists($self->{_expand_right});

  return $self->{_expand_right};
}

sub cache_size {
  my $self = shift;

  $self->{_cache_size} = shift if @_;
  $self->{_cache_size} = $DEFAULT_CACHE_SIZE if !exists($self->{_cache_size});

  return $self->{_cache_size}; 
}

sub get_data {
  my ($self, $c, $s, $e) = @_;

  die("ERROR: No chromosome specified\n") unless $c;
  die("ERROR: No start specified\n") unless $s;
  die("ERROR: No end specified\n") unless $e;

  # we use two levels of caching
  # 1) results cache for specific coords
  # 2) region cache to reduce multiple lookups for sequential "close" coords

  my $pos_string = join("_", $c, $s, $e);

  my $cache = $self->cache;

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

    my $tmp_data = $self->_get_data_uncached($c, $r_s, $r_e);

    # cache the data
    $self->_add_data_to_cache($c, $r_s, $r_e, $tmp_data) if $use_data_cache;

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

  $self->_add_result_to_cache($pos_string, \@uniq_result);
  
  return \@uniq_result;
}

sub cache {
  my $self = shift;
  my $cache = $self->{_tabix_cache} ||= {};
  return $cache;
}

sub parse_data {
  my $self = shift;
  return shift;
}

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
  my ($self, $c, $s, $e, $data) = @_;

  my $cache = $self->cache->{$c} ||= {};

  push @{$cache->{regions}}, [$s, $e];
  push @{$cache->{data}}, $data;

  # trim
  if(scalar @{$cache->{regions}} > $self->cache_size) {
    shift @{$cache->{regions}};
    shift @{$cache->{data}};
  } 
}

sub _add_result_to_cache {
  my ($self, $pos_string, $result) = @_;

  my $cache = $self->cache;

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
  my ($self, $c, $s, $e) = @_;

  my @data;

  foreach my $file(@{$self->files}) {
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
      my $parsed = $self->parse_data($line);

      push @data, [$self->identify_data($line), $parsed] if $parsed;
    }
  }

  return \@data;
}

sub _get_data_pm {
  my ($self, $c, $s, $e) = @_;

  my @data;

  foreach my $file(@{$self->files}) {
    my $tabix_obj = $self->_tabix_obj($file);
    my $valids = $self->{_valids}->{$file} ||= [$tabix_obj->getnames];

    my $iter = $tabix_obj->query($self->_get_source_chr_name($c, $file, $valids), $s, $e);
    next unless $iter && $iter->{_};

    while(my $line = $tabix_obj->read($iter)) {
      my $parsed = $self->parse_data($line);

      push @data, [$self->identify_data($line), $parsed] if $parsed;
    }
  }

  return \@data;
}

sub _get_data_cl {
  my ($self, $c, $s, $e) = @_;
  
  my @data;

  foreach my $file(@{$self->files}) {
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

      my $parsed = $self->parse_data($_);
      
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

sub FREEZE {
  my $self = shift;
  delete $self->{_tabix_obj};
}

1;

