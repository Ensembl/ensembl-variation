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
 <helpdesk@ensembl.org>.

=cut

#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::ReadCoverageCollectionAdaptor
#
# Copyright (c) 2005 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::ReadCoverageCollectionAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $rca = $reg->get_adaptor("human","variation","readcoveragecollection");
  $sa = $reg->get_adaptor("human","core","slice");
  $pa = $reg->get_adaptor("human","variation","population");

  # get read coverage in a region for a certain population in in a certain level
  $slice = $sa->fetch_by_region('chromosome', 'X', 1e6, 2e6);
  $population = $pa->fetch_by_name("UNKNOWN");
  $level = 1;
  foreach $rc (@{$vfa->fetch_all_by_Slice_Sample($slice,$population,$level)}) {
    print $rc->start(), '-', $rc->end(), "\n";
  }


=head1 DESCRIPTION

This adaptor provides database connectivity for ReadCoverageCollection objects.
Coverage information for reads can be obtained from the database using this
adaptor.  See the base class BaseFeatureAdaptor for more information.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::ReadCoverageCollectionAdaptor;

use Bio::EnsEMBL::Variation::ReadCoverageCollection;
use Bio::EnsEMBL::Mapper::RangeRegistry;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor');
my %Printable = ( "\\"=>'\\', "\r"=>'r', "\n"=>'n', "\t"=>'t', "\""=>'"' );

#my $_pack_type = "n";

=head2 fetch_all_by_Slice_SampleId

    Arg[0]      : Bio::EnsEMBL::Slice $slice
    Arg[1]      : (optional) $sample_id
    Arg[2]      : (optional) int $display_size (default 700)
    Arg[3]      : (optional) int $display_type (one of "AVERAGE" or "MAX","MIN") (default "AVERAGE")
    Arg[4]      : (optional) int $window_size
    Example     : my $reads_coverages = $rcca->fetch_all_by_Slice_SampleId($slice,$sample_id,$display_size,$display_type,$window_size);
               Window_size defines which set of pre-averaged scores to use. 
	       Valid values are 50, 500 or 5000. There is no need to define 
               the window_size because the program will select the most 
               appropriate window_size to use based on the slice_length and the
               display_size.
    Description : Gets all the read coverage collections for a given sample(strain or individual) in a given slice
    ReturnType  : listref of Bio::EnsEMBL::Variation::ReadCoverageCollection
    Exceptions  : thrown on bad arguments
    Caller      : general
    Status      : At Risk

=cut

sub fetch_all_by_Slice_SampleId{
    my $self = shift;
    my ($slice,$sample_id,$display_size,$display_type,$window_size) = @_;

     my $rcc = [];

    if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
      throw('Bio::EnsEMBL::Slice arg expected');
    }
#    if (defined $sample) {
#      if (!ref($sample) || !$sample->isa('Bio::EnsEMBL::Variation::Sample')) {
#	throw('Bio::EnsEMBL::Variation::Sample arg expected');
#      }
#      if (! defined($sample->dbID())) {
#	throw("Sample arg must have defined dbID");
#      }
#    }

   #default display_size is 700
    if (!defined $display_size) {
	$display_size = 700;
    }

    #default display_mode is AVERAGE
    if (!defined $display_type) {
	$display_type = "AVERAGE";
    }

    #set up bucket object for storing bucket_size number of scores 
    my $bucket_size = ($slice->length)/$display_size;
    
    #default window size is the largest bucket that gives at least 
    #display_size values ie get speed but reasonable resolution
    my @window_sizes = (50,500,5000);
    #my $num_windows = 100;

    #check if valid window_size
    my $found = 0;
    if (defined $window_size) {
	foreach my $win_size (@window_sizes) {
	    if ($win_size == $window_size) {
		$found = 1;
		last;
	    }
	}
	if (!$found) {
	    warning("Invalid window_size $window_size");
	    return $rcc;
	}
    }
 
    if (!defined $window_size) {
	#set window_size to be the largest for when for loop fails
	$window_size = $window_sizes[scalar(@window_sizes)-1];
	for (my $i = 1; $i < scalar(@window_sizes); $i++) {
	    if ($bucket_size < $window_sizes[$i]) {
		$window_size = $window_sizes[$i-1];
		last;
	    }
	}
    }
    
    #print "window_size is $window_size\n";

    $rcc = $self->_fetch_all_by_Slice_WindowSize_SampleId($slice,$window_size,$sample_id) if defined $sample_id;
    $rcc = $self->_fetch_all_by_Slice_WindowSize_SampleId($slice,$window_size) if !$sample_id;


    if (scalar(@$rcc ==0)) {
      return $rcc;
    }

   #Find the min and max scores for y axis scaling. Save in first
   #read_coverage_collection object
    my ($min_y_axis, $max_y_axis) =  _find_min_max_score($rcc);

    #add min and max scores to the first read_coverage_collection object
    if ((scalar @$rcc) > 0) {
	$rcc->[0]->y_axis_min($min_y_axis);
	$rcc->[0]->y_axis_max($max_y_axis);
    }
   return ($rcc);

}

sub _fetch_all_by_Slice_WindowSize_SampleId {

  my ($self,$slice,$window_size,$sample_id) = @_;

  my $rcc = [];
  my $rcs = [];
  my $reads_string_avg=[];
  my $reads_string_min=[];
  my $reads_string_max=[];

  my ($window_start,$window_end);

  my $last_sql = (defined $sample_id) ? "AND sample_id = ?" : '';
  my $sql = qq{SELECT window_size,window_start,window_end,read_coverage_string_avg,read_coverage_string_min,read_coverage_string_max,sample_id
               FROM read_coverage_collection
               WHERE seq_region_id = ?
               AND window_size = ?
               AND window_end > ?
               AND window_start < ?
               $last_sql
              };

  my $sth = $self->prepare($sql);
  $sth->execute($slice->get_seq_region_id,$window_size,$slice->start,$slice->end,$sample_id) if defined $sample_id;
  $sth->execute($slice->get_seq_region_id,$window_size,$slice->start,$slice->end) if !defined $sample_id;
  #these bind_param are not working here
  #$sth->bind_param(1,$slice->get_seq_region_id,SQL_INTEGER);
  #$sth->bind_param(2,$window_size,SQL_INTEGER);
  #$sth->bind_param(3,$slice->start,SQL_INTEGER);
  #$sth->bind_param(4,$slice->end,SQL_INTEGER);

  while (my @arrays = $sth->fetchrow_array()) { 
    $window_size = $arrays[0];
    $window_start = $arrays[1];
    $window_end = $arrays[2];
    $reads_string_avg = _unpack_strings($arrays[3]);
    $reads_string_min = _unpack_strings($arrays[4]);
    $reads_string_max = _unpack_strings($arrays[5]);
    $sample_id = $arrays[6] if $arrays[6];
    #need to find the start_offset for slice_start and end_offset for slice_end
    #in the read coverage collection row

    #print "window_size is $window_size,window_start is $window_start,window_end is $window_end,avg is @$reads_string_avg,min is @$reads_string_min,max is @$reads_string_max\n";
    for (my $j = 0; $j < @$reads_string_avg; $j++) {
      my $rc = {};
      my $seq_region_start = $window_start + $window_size * $j ;
      my $seq_region_end = $seq_region_start + $window_size -1;
      next if ($seq_region_start > $slice->end or $seq_region_end < $slice->start);
      #covert seq_region_start and seq_region_end tp a start and end relative to the slice
      my ($start,$end);
      if ($slice->strand == -1) {
	$start = $slice->end - $seq_region_end + 1;
	$end = $slice->end - $seq_region_start + 1;
      }
      else {
	$start = $seq_region_start - $slice->start + 1;
	$end = $seq_region_end - $slice->start +1;
      }
      $rc->{'start'} = $start;
      $rc->{'end'} = $end;
      $rc->{'seq_region_start'} = $seq_region_start;
      $rc->{'seq_region_end'} = $seq_region_end;
      $rc->{'strand'} = $slice->strand;
      $rc->{'read_coverage_avg'} = $reads_string_avg->[$j];
      $rc->{'read_coverage_min'} = $reads_string_min->[$j];
      $rc->{'read_coverage_max'} = $reads_string_max->[$j];
      push @{$rcs}, $rc;
    }
  }#end while loop

  foreach my $rc (@$rcs) {
    push @$rcc,$self-> _create_feature_fast('Bio::EnsEMBL::Variation::ReadCoverageCollection',
					    {'adaptor' => $self,
					     'slice' => $slice,
					     'start' => $rc->{'start'},
					     'end'   => $rc->{'end'},
					     'strand' => $rc->{'strand'},
					     'window_size' => $window_size,
					     'window_start' => $window_start,
					     'window_end' => $window_end,
					     'seq_region_start' => $rc->{'seq_region_start'},
					     'seq_region_end' => $rc->{'seq_region_end'},
					     'seq_region_strand' => 1,
					     'sample_id' => (defined $sample_id) ? $sample_id : '',
					     'read_coverage_avg' => $rc->{'read_coverage_avg'},
					     'read_coverage_min' => $rc->{'read_coverage_min'},
					     'read_coverage_max' => $rc->{'read_coverage_max'},
					    });
  }

  #sort into numerical order based on position
  my @sorted_features = sort {$a->{start} <=> $b->{start}} @$rcc;
  return\@sorted_features;
}

sub _unpack_strings {

  my $string = shift;

  my $pack_type = "n" x (length($string)/2);
  my @arrays = unpack($pack_type,$string);
  return \@arrays;
}

sub escape ($) {
  local $_ = ( defined $_[0] ? $_[0] : '' );
  s/([\r\n\t\\\"])/\\$Printable{$1}/sg;
  return $_;
}


sub _find_min_max_score {

    my ($scores) = @_;
    my $min; 
    my $max;

    foreach my $score (@$scores) {
	#find min and max of diff scores
	if (defined $score->read_coverage_min) {
	    #if min hasn't been defined yet, then define min and max
	    unless (defined $min and defined $max) {
		$min = $score->read_coverage_min;
		$max = $score->read_coverage_max;
	    }
	    if ($min > $score->read_coverage_min) {
		$min = $score->read_coverage_min;
	    }
	    if ($max < $score->read_coverage_max) {
		$max = $score->read_coverage_max;
	    }
	}
    }

    return ($min, $max);
}


1;
