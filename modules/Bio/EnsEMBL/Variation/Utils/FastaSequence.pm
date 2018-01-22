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

Bio::EnsEMBL::Variation::Utils::FastaSequence

=head1 SYNOPSIS

  use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);

  my $fasta_db = setup_fasta(
    -FASTA => 'genome.fa',
    -ASSEMBLY => 'GRCh38'
  );

  This module overrides the sequence fetching method in
  Bio::EnsEMBL::Slice with calls to either Bio::DB::Fasta or Bio::DB::HTS::Faidx

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::Utils::FastaSequence;

use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);
use Bio::EnsEMBL::Slice;
use Bio::DB::Fasta;

use vars qw(@ISA @EXPORT_OK);
@ISA = qw(Exporter);

@EXPORT_OK = qw(
  &setup_fasta
  &revert_fasta
  &clear_fasta_cache
);

our $DEBUG = 0;

our $CAN_USE_FAIDX;

# predefine PAR regions for human assemblies
our %PARS = (
  GRCh37 => [
    {
      start => 10001, # start of region on Y
      end => 2649520, # end of region on Y
      adj => 50000    # how much to adjust Y coords by to get X coords
    }
  ],
  GRCh38 => [
    {
      start => 10001,   # start of region on Y
      end   => 2781479, # end of region on Y
      adj   => 0        # how much to adjust Y coords by to get X coords
    },
    {
      start => 56887903, # start of region on Y
      end   => 57217415, # end of region on Y
      adj   => 98813480  # how much to adjust Y coords by to get X coords
    }
  ],
);

BEGIN {
  if (eval { require Bio::DB::HTS::Faidx; 1 }) {
    $CAN_USE_FAIDX = 1;
  }
  else {
    $CAN_USE_FAIDX = 0;
  }
}

=head2 setup_fasta

  Arg [-FASTA]    : string - the path to the FASTA file. May be bgzipped
                    if Faidx is installed

  Arg [-TYPE]     : (optional) string - one of Faidx or Bio::DB::Fasta.
                    Defaults to Faidx if installed.

  Arg [-ASSEMBLY] : (optional) string - assembly name; used for Homo sapiens
                    to assign sequence from X PAR regions to corresponding Y
                    regions. These are otherwise "N"s in the FASTA file.

  Arg [-PARS]     : (optional) arrayref - manually define PARs as follows:
    $pars = [
      {
        start => 10001,   # start of region on Y
        end   => 2781479, # end of region on Y
        adj   => 0        # how much to adjust Y coords by to get X coords
      }
    ];

  Arg [-OFFLINE]  : (optional) boolean - set to a true value to prevent
                    falling back to using the original
                    Bio::EnsEMBL::Slice->seq() method if a sequence name
                    is not found in the FASTA file.

  Arg [-SYNONYMS] : (optional) hashref - two-level hashref of chromosome
                    name synonyms e.g.

                    {
                      foo => {
                        bar => 1,
                        car => 1,
                      }
                      goo => {
                        far => 1,
                        star => 1,
                      }
                    }

  Example         :
    $fasta_db = setup_fasta(
      -fasta    => '/path/to/homo_sapiens.fa.gz',
      -assembly => 'GRCh38',
      -type     => 'Bio::DB::HTS::Faidx',
    );

  Description: Set up a FASTA DB to override the seq() method in Bio::EnsEMBL::Slice
  Returntype : Bio::DB::Fasta or Bio::DB::HTS::Faidx
  Exceptions : no file, wrong path, wrong file type
  Caller     : VEP, VariationEffect pipeline
  Status     : Stable

=cut

sub setup_fasta {
  my (
    $fasta,
    $assembly,
    $pars,
    $offline,
    $type,
    $synonyms
  ) = rearrange([qw(
    FASTA
    ASSEMBLY
    PARS
    OFFLINE
    TYPE
    SYNONYMS
  )], @_);

  throw("ERROR: No FASTA file specified\n") unless $fasta;
  throw("ERROR: Specified FASTA file/directory $fasta not found\n") unless -e $fasta;

  # work out type
  if($type) {
    if($type eq 'Bio::DB::HTS::Faidx' && !$CAN_USE_FAIDX) {
      throw("ERROR: Bio::DB::HTS::Faidx not installed\n");
    }
    elsif($type ne 'Bio::DB::Fasta') {
      throw("ERROR: Unrecognised index type $type\n");
    }
  }

  # if given FASTA file is gzipped, have a quick look for the uncompressed file
  if($fasta =~ /\.gz$/ && !$CAN_USE_FAIDX) {

    # but first check that the unpacked file doesn't exist
    my $unpacked_fa = $fasta;
    $unpacked_fa =~ s/\.gz$//;

    # if it does, we can use it instead
    if(-e $unpacked_fa) {
      $fasta = $unpacked_fa;
      warning("Switching to using $unpacked_fa\n");
    }
  }

  $type ||= $CAN_USE_FAIDX ? 'Bio::DB::HTS::Faidx' : 'Bio::DB::Fasta';

  throw("ERROR: Cannot index bgzipped FASTA file with Bio::DB::Fasta\n") if $type eq 'Bio::DB::Fasta' && $fasta =~ /\.gz$/;

  my $fasta_db = _get_fasta_db($fasta, $type);

  # dont re-redefine
  unless($Bio::EnsEMBL::Slice::_fasta_redefined) {

    # "back up" the original seq and subseq method
    # that way we can fall back to it later
    no warnings 'redefine';
    *Bio::EnsEMBL::Slice::_fasta_old_db_seq = \&Bio::EnsEMBL::Slice::seq unless $offline;

    no warnings 'redefine';
    *Bio::EnsEMBL::Slice::_fasta_old_db_subseq = \&Bio::EnsEMBL::Slice::subseq unless $offline;

    # redefine Slice's seq() and subseq() methods
    # to a new method that uses the FASTA sequence
    # and some nice caching
    no warnings 'redefine';
    *Bio::EnsEMBL::Slice::seq = _new_slice_seq();

    no warnings 'redefine';
    *Bio::EnsEMBL::Slice::subseq = _new_slice_seq();

    # this method does the actual fetching
    # separate so it can easily use Bio::DB::Fasta or Faidx
    *Bio::EnsEMBL::Slice::_raw_seq = _raw_seq();

    # store a flag to indicate we have redefined
    # saves us setting up again
    $Bio::EnsEMBL::Slice::_fasta_redefined = 1;
  }

  # we need to tell Slice about the fasta DB
  # may as well use a package variable
  no warnings 'once';
  $Bio::EnsEMBL::Slice::fasta_db = $fasta_db;

  # we also need to tell it about PARs, ugh hacky
  $pars ||= $PARS{$assembly} if $assembly && defined($PARS{$assembly});

  $Bio::EnsEMBL::Slice::_fasta_PARs = $pars if $pars;
  $Bio::EnsEMBL::Slice::_fasta_synonyms = $synonyms || {};

  return $fasta_db;
}


=head2 revert_fasta
  Description: Revert Bio::EnsEMBL::Slice->seq to original behaviour
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub revert_fasta {
  no warnings 'redefine';
  *Bio::EnsEMBL::Slice::seq = \&Bio::EnsEMBL::Slice::_fasta_old_db_seq;

  no warnings 'redefine';
  *Bio::EnsEMBL::Slice::subseq = \&Bio::EnsEMBL::Slice::_fasta_old_db_subseq;

  undef $Bio::EnsEMBL::Slice::fasta_db;
  undef $Bio::EnsEMBL::Slice::_fasta_PARs;
  undef $Bio::EnsEMBL::Slice::_fasta_synonyms;
  clear_fasta_cache();

  $Bio::EnsEMBL::Slice::_fasta_redefined = 0;
}


=head2 clear_fasta_cache
  Description: Manually clear the cache used by the sequence fetching
               code. Not necessary in normal use as the cache clears
               itself periodically.
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub clear_fasta_cache {
  undef $Bio::EnsEMBL::Slice::_fasta_sequence_cache;
}


# creates the FASTA DB object
# returns either a Bio::DB::Fasta or a Bio::DB::HTS::Faidx
# uses a lock file to avoid partial indexes getting used
sub _get_fasta_db {
  my $fasta = shift;
  my $type = shift;

  # check lock file
  my $lock_file = $fasta;
  $lock_file .= -d $fasta ? '/.vep.lock' : '.vep.lock';

  # lock file exists, indexing failed
  if(-e $lock_file) {
    print STDERR "Lock file found, removing to allow re-indexing\n" if $DEBUG;
    for(qw(.fai .index .gzi /directory.index /directory.fai .vep.lock)) {
      unlink($fasta.$_) if -e $fasta.$_;
    }
  }

  my $index_exists = 0;

  for my $fn(map {$fasta.$_} qw(.fai .index .gzi /directory.index /directory.fai)) {
    if(-e $fn) {
      $index_exists = 1;
      last;
    }
  }

  # create lock file
  unless($index_exists) {
    open LOCK, ">$lock_file" or throw("ERROR: Could not write to FASTA lock file $lock_file\n");
    print LOCK "1\n";
    close LOCK;
    print STDERR "Indexing $fasta\n" if $DEBUG;
  }

  # run indexing
  my $fasta_db;
  if($type eq 'Bio::DB::HTS::Faidx' && $CAN_USE_FAIDX) {
    $fasta_db = Bio::DB::HTS::Faidx->new($fasta);
  }
  elsif(!$type || $type eq 'Bio::DB::Fasta') {
    $fasta_db = Bio::DB::Fasta->new($fasta);
  }
  else {
    throw("ERROR: Don't know how to index with type $type\n");
  }
  print STDERR "Finished indexing\n" if $DEBUG;

  throw("ERROR: Unable to create FASTA DB\n") unless $fasta_db;

  # remove lock file
  unlink($lock_file) unless $index_exists;

  return $fasta_db;
}

# the method called in place of Bio::EnsEMBL::Slice->seq()
# most of the logic in this method is for caching sequence
# sequence is cached in 1MB range around the current slice
sub _new_slice_seq {
  return sub {
    my ($self, $start, $end, $strand) = @_;
    my ($seq, $length) = ('', 0);

    $start = $start ? ($self->start + $start) - 1 : $self->start;
    $end   = $end   ? ($self->start + $end) - 1   : $self->end;
    $strand = defined($strand) ? $strand * $self->strand : $self->strand;
    my $sr_name = $self->seq_region_name;

    # indels
    return "" if $start > $end;

    my $fasta_db = $Bio::EnsEMBL::Slice::fasta_db;
    my $fa_length = $fasta_db->length($sr_name);
    unless($fa_length && $fa_length > 0) {
      foreach my $alt(
        $sr_name =~ /^chr/i ?
          (substr($sr_name, 3), keys %{$Bio::EnsEMBL::Slice::_fasta_synonyms->{substr($sr_name, 3)} || {}}) :
          (
            'chr'.$sr_name,
            'CHR'.$sr_name,
            keys %{$Bio::EnsEMBL::Slice::_fasta_synonyms->{'chr'.$sr_name} || {}},
            keys %{$Bio::EnsEMBL::Slice::_fasta_synonyms->{'CHR'.$sr_name} || {}}
          ),
        keys %{$Bio::EnsEMBL::Slice::_fasta_synonyms->{$sr_name} || {}}
      ) {
        $fa_length = $fasta_db->length($alt);

        if($fa_length && $fa_length > 0) {
          print STDERR "USING SYNOYM $alt FOR $sr_name\n" if $DEBUG;
          $sr_name = $alt;
          last;
        }
      }

      if(!($fa_length && $fa_length > 0) && $self->can('_fasta_old_db_seq')) {
        print STDERR "USING DATABASE\n" if $DEBUG;
        return
          scalar(@_) > 1 ?
          $self->_fasta_old_db_subseq($start, $end, $strand) :
          $self->_fasta_old_db_seq();
      }
    }

    my $cache = $Bio::EnsEMBL::Slice::_fasta_sequence_cache->{$sr_name} ||= [];

    # find matching regions, if any
    my $region;
    for(my $i=0; $i<scalar(@$cache); $i++) {
      my $tmp_region = $cache->[$i];
      if(overlap($start, $end, $tmp_region->{start}, $tmp_region->{end})) {
        $region = $tmp_region;
        last;
      }
    }

    if($region) {
      my ($region_start, $region_end, $region_strand) = ($region->{start}, $region->{end}, $region->{strand});

      print STDERR "USING CACHE\n" if $DEBUG;

      my $updated = 0;

      ## Get any extra sequence we need
      ## What we do is append sequence to the current region
      ## This means we fetch it on the same strand as the original region
      ## even if the requested sequence is for the opposite strand,
      ## the idea being that most times you will be requesting sequence for the same strand in a given region
      ## revcomp is slow so trying to reduce the number of calls to it
      ## We then revcomp later if you are actually asking for the opposite

      # overhang 3'
      if($start < $region_start) {

        # get missing sequence
        my $missing_seq = $self->_raw_seq($sr_name, $start, $region_start - 1, $region_strand);

        # reverse strand: append to current seq
        if($region_strand < 0) {
          $region->{seq} .= $missing_seq;
        }

        # forward strand: prepend to current seq
        else {
          $region->{seq} = $missing_seq.$region->{seq};
        }

        # update region start
        $region->{start} = $start;
        $updated = 1;
      }

      # overhang 5'
      if($end > $region_end) {

        # get missing sequence
        my $missing_seq = $self->_raw_seq($sr_name, $region_end + 1, $end, $region_strand);

        # reverse strand: prepend to current seq
        if($region_strand < 0) {
          $region->{seq} = $missing_seq.$region->{seq};
        }

        # forward strand: append to current seq
        else {
          $region->{seq} .= $missing_seq;
        }

        # update region end
        $region->{end} = $end;
        $updated = 1;
      }

      # delete revcomped seq if coords updated
      delete $region->{revseq} if $updated;

      # now reverse comp if requested strand opposite
      if($strand != $region_strand) {

        # get and cache revcomped seq
        if(!$region->{revseq}) {
          my $revseq = $region->{seq};
          reverse_comp(\$revseq);
          $region->{revseq} = $revseq;
        }

        my $substr_start = $region_strand < 0 ? $start - $region->{start} : $region->{end} - $end;
        $seq = substr($region->{revseq}, $substr_start, ($end - $start) + 1);
      }

      else {
        my $substr_start = $region_strand < 0 ? $region->{end} - $end : $start - $region->{start};
        $seq = substr($region->{seq}, $substr_start, ($end - $start) + 1);
      }
    }

    # no sequence in cache
    # we need to fetch
    if(!$seq) {

      # do raw fetch
      $seq = $self->_raw_seq($sr_name, $start, $end, $strand);

      # prune the cache
      # do this before we add the new one to do one fewer greps
      @$cache = grep {overlap($_->{start}, $_->{end}, $start - 1e6, $end + 1e6)} @$cache;

      # use unshift on the cache array
      # this way we'll likely find this sooner next time
      unshift @$cache, {
        seq    => $seq,
        start  => $start,
        end    => $end,
        strand => $strand
      };
    }

    return $seq;
  };
}

# does the actual sequence fetch
sub _raw_seq {
  return sub {
    my ($self, $sr_name, $start, $end, $strand) = @_;

    # indels
    return "" if $start > $end;

    # get fasta DB from package variable
    my $fasta_db = $Bio::EnsEMBL::Slice::fasta_db;

    throw("ERROR: No FASTA DB defined\n") unless $fasta_db;

    my $requested_length = ($end - $start) + 1;

    ## handle PARS
    my ($pre, $post) = ('', '');

    print STDERR "USING FASTA\n" if $DEBUG;

    # There's an assumption here that we won't get asked for a sequence that overlaps more than one PAR
    # Let's face it though that would be >54Mbp so unlikely!
    if(
      $sr_name eq 'Y' &&
      $Bio::EnsEMBL::Slice::_fasta_PARs &&
      (my ($par) = grep {overlap($_->{start}, $_->{end}, $start, $end)} @{$Bio::EnsEMBL::Slice::_fasta_PARs})
    ) {

      # check for partial overlap at 5' end
      if($start < $par->{start}) {

        # get chunk of non-PAR sequence
        my $tmp_seq = $self->_raw_seq('Y', $start, $par->{start} - 1, $strand);

        # then make it a prefix (+ve) or a suffix (-ve) depending on requested strand
        if($strand > 0) {
          $pre = $tmp_seq;
        }
        else {
          $post = $tmp_seq;
        }

        # adjust start for the remaining sequence to be requested
        $start = $par->{start};
      }

      # check for partial overlap at 3' end
      if($end > $par->{end}) {

        # get chunk of non-PAR sequence
        my $tmp_seq = $self->_raw_seq('Y', $par->{end} + 1, $end, $strand);

        # then make it a suffix (+ve) or a prefix (-ve) depending on requested strand
        if($strand > 0) {
          $post = $tmp_seq;
        }
        else {
          $pre = $tmp_seq;
        }

        # adjust end for the remaining sequence to be requested
        $end = $par->{end};
      }

      # rest of seq is fetched from X on adjusted coords
      $sr_name = 'X';
      $start  += $par->{adj};
      $end    += $par->{adj};
    }

    # FASTA file does not contain this chrom
    my $fa_length = $fasta_db->length($sr_name);
    unless($fa_length && $fa_length > 0) {
      warning("WARNING: FASTA does not contain sequence for ".$sr_name."\n");
      return 'N' x (($end - $start) + 1);
    }

    # requested slice beyond limits?
    if($start < 1) {
      print STDERR "SLICE HAS START < 1\n" if $DEBUG;

      my $tmp_seq = 'N' x (1 - $start);

      if($strand > 0) {
        $pre = $tmp_seq.$pre;
      }
      else {
        $post .= $tmp_seq;
      }

      $start = 1;
    }

    if($end > $fa_length) {
      print STDERR "SLICE HAS END > SEQ LENGTH\n" if $DEBUG;

      my $tmp_seq = 'N' x ($end - $fa_length);

      if($strand > 0) {
        $post .= $tmp_seq;
      }
      else {
        $pre = $tmp_seq.$post;
      }

      $end = $fa_length;
    }

    my ($seq, $length);

    # different modules have different calls
    if($fasta_db->isa('Bio::DB::HTS::Faidx')) {
      my $location_string = $sr_name.":".$start."-".$end ;
      ($seq, $length) = $fasta_db->get_sequence($location_string);
      $seq = uc($seq);
    }
    elsif($fasta_db->isa('Bio::DB::Fasta')) {
      $seq = uc($fasta_db->seq($sr_name, $start => $end));
    }
    else {
      throw("ERROR: Don't know how to fetch sequence from a ".ref($fasta_db)."\n");
    }

    reverse_comp(\$seq) if $seq && defined($strand) && $strand < 0;

    # add on PAR overlapped chunks
    $seq = $pre.$seq.$post;

    # default to a string of Ns if we couldn't get sequence
    if(!$seq) {
      warning(sprintf("WARNING: Could not obtain sequence for %s:%i-%i:%i\n", $sr_name, $start, $end, $strand));
      $seq = 'N' x (($self->end - $self->start) + 1);
    }

    print STDERR "GOT SEQ LENGTH ".length($seq)." REQUESTED ".$requested_length."\n" if $DEBUG;

    return $seq;
  }
}

1;
