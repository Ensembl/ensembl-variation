use strict;
use warnings;

package LRGMapping;

#Êuse lib '/nfs/acari/dr2/projects/src/ensembl/ensembl/modules';

use Bio::EnsEMBL::LRGSlice;
use Bio::EnsEMBL::MappedSliceContainer;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use LRG::LRG;
use LRG::LRGImport;
use List::Util qw (min max);


# global variables
#our $input_dir = '/tmp';
#our $output_dir = '/tmp';

our $input_dir = '/lustre/scratch103/ensembl/pl4/tmp';
our $output_dir = '/lustre/scratch103/ensembl/pl4/tmp';

our $target_dir = '/lustre/work1/ensembl/yuan/SARA/human/ref_seq_hash';
our %rec_seq;

#ÊKeep different db connections for a core db with write access to and one with read-only access
our $dbCore_rw;
our $dbCore_ro;

our $lrg_name;
our $registry_file;
our $current_assembly;

our $mapping_num = 1;

sub mapping {
  my $lrg_seq = shift;
  
  my $maps = map_to_genome($lrg_seq);
  my $mainmap = join_mappings($maps);
  
  return $mainmap;
}

sub map_to_genome {
  my $q_seq = shift;
  my $lrg_offset = shift;
  
  my $SPLIT_SIZE = 25000;
  my $SPLIT_OVERLAP = 1000;
  my $MIN_QUERY_LENGTH = 15000;
  
  my $q_len = length($q_seq);
  $lrg_offset ||= 0;
  
  my @maps;
  
  # If the remainder of the split sequence will not be smaller than the minimum query length, it is ok to split the sequence
  if ($q_len + $SPLIT_OVERLAP > $SPLIT_SIZE + $MIN_QUERY_LENGTH) {
    # Chop off a split at the 5' end 
    my $q_split = substr($q_seq,0,$SPLIT_SIZE);
    #ÊAppend a small overlap to the remainder
    my $q_rest = substr($q_split,-$SPLIT_OVERLAP) . substr($q_seq,$SPLIT_SIZE);
    
    # Map the split to the genome and add it to the array
    @maps = (@maps,@{map_to_genome($q_split,$lrg_offset)});
    #ÊSend the rest of the sequence back to this method to do a recursion
    @maps = (@maps,@{map_to_genome($q_rest,$lrg_offset+($SPLIT_SIZE-$SPLIT_OVERLAP))});
    
    return \@maps;
  }
  # Else, we map the supplied piece
  @maps = @{ssaha_mapping($q_seq)};
  
  # For all mappings, add the LRG offset to the LRG coordinates
  foreach my $map (@maps) {
    $map->{'lrg_start'} += $lrg_offset;
    $map->{'lrg_end'} += $lrg_offset;
    # Do this also for all pairs in the pairs array
    foreach my $pair (@{$map->{'pairs'}}) {
      $pair->[1] += $lrg_offset;
      $pair->[2] += $lrg_offset;
    }
  }
  
  return \@maps;
}

# Join together overlapping mappings to one main map
sub join_mappings {
  my $mapref = shift;
  
  # Order the mappings in ascending LRG start coords order
  my @mappings = sort {$a->{'lrg_start'} <=> $b->{'lrg_start'}} @{$mapref};
  
  my %mainmap;
  while (my $submap = shift(@mappings)) {
    #ÊInitialize the mainmap by setting it to the first mapped piece
    if (!defined($mainmap{'lrg_id'})) {
      %mainmap = %{$submap};
      next;
    }
    #ÊCheck that this piece is overlapping the previous one (or possibly lying back-to-back in case of no overlap)
    die("Mapped pieces are not overlapping!") unless ($submap->{'lrg_end'} > $mainmap{'lrg_end'} && $submap->{'lrg_start'} - $mainmap{'lrg_end'} <= 1);
    # Do a check to make sure that the mappings are not on opposite strands
    die("Mapped pieces are on opposite strands!") unless ($submap->{'strand'} == $mainmap{'strand'});
    
    #ÊSort the mainmap's pairs in descending LRG coord order
    my @mainpairs = sort {$b->[2] <=> $a->[2]} @{$mainmap{'pairs'}};
    
    #ÊAdd pairs that fall on or beyond the mainmap's border to the mainmap's pair array
    while (my $subpair = shift(@{$submap->{'pairs'}})) {
      # Go to the next if the pair ends within the overlap with the main map
      next if ($subpair->[2] <= $mainmap{'lrg_end'});
      #ÊIf the pair starts within the overlap and extends beyond, it should be merged with an already existing pair
      # In the special case of no overlap between mappings, 'DNA' pairs that lie back-to-back and have no gap in between should be joined
      if (
	  ($subpair->[1] <= $mainmap{'lrg_end'} && $subpair->[2] > $mainmap{'lrg_end'}) ||
	  ($submap->{'lrg_start'} == ($mainmap{'lrg_end'} + 1) && $subpair->[0] eq 'DNA' && $subpair->[1] == $submap->{'lrg_start'})
	 ) {
	foreach my $mainpair (@mainpairs) {
	  if (
	      $subpair->[0] eq $mainpair->[0] &&
	      (
	       ($subpair->[1] <= $mainmap{'lrg_end'}) ||
	       ($submap->{'lrg_start'} == ($mainmap{'lrg_end'} + 1) && $subpair->[1] == $submap->{'lrg_start'} && $mainpair->[2] == $mainmap{'lrg_end'})
	      )
	     ) {
	    $mainpair->[2] = $subpair->[2];
	    $mainpair->[3] = min($mainpair->[3],$subpair->[3]);
	    $mainpair->[4] = max($mainpair->[4],$subpair->[4]);
	    last;
	  }
	  # If the mainpair ends short of the end of the mainmap, we don't need to loop anymore (since they are sorted)
	  last if ($mainpair->[2] < $mainmap{'lrg_end'});
	}
      }
      # If the pair starts beyond the mainmap, add it to the mainmap's pair array
      push(@mainpairs,$subpair) if ($subpair->[1] > $mainmap{'lrg_end'} )
    }
    
    # Update the mainmap's coordinates
    $mainmap{'lrg_end'} = $submap->{'lrg_end'};
    $mainmap{'chr_start'} = min($mainmap{'chr_start'},$submap->{'chr_start'});
    $mainmap{'chr_end'} = max($mainmap{'chr_end'},$submap->{'chr_end'});
    $mainmap{'score'} += $submap->{'score'};
    
    # Replace the mainmap's pair array
    $mainmap{'pairs'} = \@mainpairs;
  }
  
  return \%mainmap;
}

# Map the supplied sequence to the genome using ssaha2 and return an arrayref with hashes containing mapping details as well as a pairs array
sub ssaha_mapping {
  my $sequence = shift;
  
  my $name = "temp$$".$mapping_num++;
  my $input_file_name = $name.'.fa';
  	
  open OUT, ">$input_dir/$name\.fa" or die $!;
  print OUT '>', $name, "\n", $sequence;
  close OUT;
  	
  my $queue = "-q normal -R'select[mem>4000] rusage[mem=4000]' -M4000000";
  	
  my $output_file_name = "ssaha2_output\_$input_file_name";
  my $input_file = "$input_dir/$input_file_name";
  my $output_file ="$output_dir/$output_file_name";
  	
  my ($subject, %rec_find, %input_length, %done);
  	
  $subject = "$target_dir/ref";
  $rec_seq{$name} = $sequence;
  my $seqobj = Bio::PrimarySeq->new(-id => $name, -seq => $rec_seq{$name});
  $rec_seq{$name} = $seqobj;
  	
  bsub_ssaha_job($queue,$input_file,$output_file,$subject);
  
  my $call = "bsub -o $input_dir/$$.waiting.out -P ensembl-variation -K -w 'done($input_file\_ssaha_job)' -J waiting_process sleep 1"; #waits until all ssaha jobs have finished to continue
  system($call);
  
  # check the output file
  open IN, $output_file or die("Could not read from output file $output_file\n");
  
  my $ok = 0;
  
  while(<IN>) {
    if(/^Query/) {
  	$ok = 1;
  	last;
    }
  }
  close IN;
  
  die "Error in mapping process - no data found in ssaha output file $output_file\n" unless $ok;
  
  my $mapping = parse_ssaha2_out($output_file);
  
  return $mapping;
}

# subroutine to submit a job to the farm
sub bsub_ssaha_job {
  my ($queue, $input_file, $output_file, $subject) = @_;
  
  my $host = `hostname`;
  chomp $host;
  
  my $ssaha_command = "/nfs/users/nfs_y/yuan/ensembl/src/ensembl-variation/scripts/ssahaSNP/ssaha2/ssaha2_v1.0.9_x86_64/ssaha2";
  $ssaha_command .= " -align 1 -kmer 12 -seeds 4 -cut 1000 -output vulgar -depth 10 -best 1 -save $subject $input_file";
#  my $call = "echo '$ssaha_command; scp $output_file $host:$input_dir/' | bsub -E 'scp $host:$input_file $input_dir/' -J $input_file\_ssaha_job -P ensembl-variation $queue -e $output_dir/error_ssaha -o $output_file";
  my $call = "echo '$ssaha_command' | bsub -J $input_file\_ssaha_job -P ensembl-variation $queue -e $output_dir/error_ssaha -o $output_file";
	
  system ($call);
#  print $call, "\n";
}

# subroutine to parse the output from ssaha2
sub parse_ssaha2_out {

#   open OUT, ">$output_dir/mapping_file_$input_file_name" or die "can't open output_file : $!\n";
  my $output_file = shift;

  open SSAHA, "$output_file" or die "can't open $output_file : $!\n";
  my @rec_find;
  
  my $q_seq;
  my $t_seq;
  my $vulgar;
  my $data;
  
  while (<SSAHA>) {
    chomp;
    my $line = $_;
    
    # If the vulgar annotation line is encountered, store it. If one already is stored, parse it with the collected sequences  
    if ($line =~ m/^vulgar\:/) {
      if (defined($vulgar)) {
	$data = parse_vulgar_string($vulgar,$q_seq,$t_seq);
	if (defined($data)) {
	  push(@rec_find,$data);
	}
	$q_seq = undef;
	$t_seq = undef;
      }
      $vulgar = $line;
    }
    elsif ($line =~ m/^Query|Sbjct/) {
      my ($id,$seq) = $line =~ m/^(Query|Sbjct)\s+\d+\s+([^\s]+)\s+\d+$/;
      
      # Remove gap characters
      $seq =~ s/\-//g;
      if ($id eq 'Query') {
	$q_seq .= $seq;
      }
      else {
	$t_seq .= $seq;
      }
    }
  }
  close SSAHA;

  #ÊParse the stored vulgar string and sequences
  if (defined($vulgar)) {  
    $data = parse_vulgar_string($vulgar,$q_seq,$t_seq);
    if (defined($data)) {
      push(@rec_find,$data);
    }
  }
  
  # For now, die if we could not get exactly one mapping
  die("The LRG sequence was not uniquely mapped!") unless (scalar(@rec_find) == 1);
  
  return \@rec_find;
  
    #print ;
=head    
    next if ($_ !~ /^vulgar/);
    if (/^vulgar\:\s+(\S+)\s+(\d+)\s+(\d+)\s+(\+|\-)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\+|\-)\s+(\d+)\s+(.*)$/) {
      my $h = undef;
      $h->{'q_id'} = $1;
      $h->{'q_start'} = $2;
      $h->{'q_end'} = $3;
      $h->{'t_id'} = (split /\-/, $5)[0];
      $h->{'t_start'} = $6;
      $h->{'t_end'} = $7;
      $h->{'score'} = $9;
      $h->{'match'} = $10;
      my $q_strand = $4;
      my $t_strand = $8;
      $q_strand = ($q_strand =~ /\-/) ? -1 : 1;
      $t_strand = ($t_strand =~ /\-/) ? -1 : 1;
      $h->{'q_strand'} = $q_strand;
      $h->{'t_strand'} = $t_strand;
      push @rec_find, $h if $h->{'q_id'};
    }
  }
  close SSAHA;
  
  #unlink $output_file;
  
  # if we have more than one hit
  if(scalar @rec_find > 1 && $rec_find[0]->{'score'} == $rec_find[1]->{'score'}) {
    die "Mapped ",scalar @rec_find, "times\n";
  }
  my $mapping = make_feature_pair($rec_find[0]);
  return $mapping;
=cut
}

# Attempts to parse a vulgar string as outputted from ssaha2. Will return a reference to a hash with the following fields:
#  lrg_id
#  chr_name
#  lrg_start
#  lrg_end
#Ê chr_start
#Ê chr_end
#  strand
#  score
#  pairs
sub parse_vulgar_string {
  my $vulgar = shift;
  my $q_seq = shift;
  my $t_seq = shift;
  
  my ($q_id,$q_start,$q_end,$q_strand,$t_id,$t_start,$t_end,$t_strand,$score,$match) = $vulgar =~ m/^vulgar\:\s+([^\s]+)\s+(\d+)\s+(\d+)\s+([\+\-])\s+([^\s]+)\s+(\d+)\s+(\d+)\s+([\+\-])\s+(\d+)\s+(.+)$/;
  my %data;
  
  unless(
    defined($q_id) &&
    defined($q_start) &&
    defined($q_end) &&
    defined($q_strand) &&
    defined($t_id) &&
    defined($t_start) &&
    defined($t_end) &&
    defined($t_strand) &&
    defined($score) &&
    defined($match)
  ) {
    return undef;
  }
  
  $q_strand = ($q_strand eq '+' ? 1 : -1);
  $t_strand = ($t_strand eq '+' ? 1 : -1);
  
  # If the match is on the minus strand relative to the LRG, flip everything around
  if ($q_strand < 0) {
    ($q_start,$q_end) = ($q_end,$q_start);
    $t_strand *= -1;
    $match = join(" ",reverse(split(/\s+/,$match)));
    reverse_comp(\$q_seq);
    reverse_comp(\$t_seq);
  }
  
  # Get the chromosome name
  my ($t_name) = $t_id =~ m/^([^\-]+)\-/;
  
  $data{'lrg_id'} = $q_id;
  $data{'chr_name'} = $t_name;
  $data{'lrg_start'} = $q_start;
  $data{'lrg_end'} = $q_end;
  $data{'chr_start'} = $t_start;
  $data{'chr_end'} = $t_end;
  $data{'strand'} = $q_strand;
  $data{'score'} = $score;
  
  my @pairs;
  
  #ÊA vulgar label 'M' means that the match is of equal length* while 'G' means an indel
  # *At least this is what I think. I do an extra step to check this and if it's not true, the script will die.
  my $q_offset = 0;
  my $t_offset = 0;
  while ($match =~ m/([MG\d]+)\s(\d+)\s([MG\d]+)/g) {
    my ($lbl,$q_len,$t_len);
    if ($q_strand < 0) {
      ($lbl,$q_len,$t_len) = ($3,$2,$1);
    }
    else {
      ($lbl,$q_len,$t_len) = ($1,$2,$3);
    }
    
    my $pair;
    if ($lbl eq 'M') {
      die('The length of query and target in match pair is different!') unless ($q_len == $t_len);
      
      # Store a coordinate pair with label 'DNA'. These can be used as mapping spans and chunks in the assembly mapping table
      $pair = [
		'DNA',
		$q_start + $q_offset,
		$q_start + $q_offset + $q_len - 1,
		($t_strand > 0 ? $t_start + $t_offset : $t_end - $t_offset - $t_len + 1),
		($t_strand > 0 ? $t_start + $t_offset + $t_len - 1 : $t_end - $t_offset),
		$q_strand
	      ];
      push(@pairs,$pair);
      
      # Go through each nucleotide in the query and target and store SNPs
      my $i = 0;
      while ($i < $q_len) {
	if (substr($q_seq,$q_start + $q_offset + $i - 1,1) ne substr($t_seq,1 + $t_offset + $i - 1,1)) {
	  $pair = [
	    $lbl,
	    $q_start + $q_offset + $i,
	    $q_start + $q_offset + $i,
	    ($t_strand > 0 ? $t_start + $t_offset + $i : $t_end - $t_offset - $i),
	    ($t_strand > 0 ? $t_start + $t_offset + $i : $t_end - $t_offset - $i),
	    $q_strand,
	    substr($q_seq,$q_start + $q_offset + $i - 1,1),
	    substr($t_seq,1 + $t_offset + $i - 1,1)
	  ];
	  push(@pairs,$pair);
	}
	$i++;
      }
    }
    else {
      my $gap_len = max($q_len,$t_len);
      my $gap_string = '-' x $gap_len;
      $pair = [
	$lbl,
	$q_start + $q_offset,
	$q_start + $q_offset + $q_len - 1,
	($t_strand > 0 ? $t_start + $t_offset : $t_end - $t_offset - $t_len + 1),
	($t_strand > 0 ? $t_start + $t_offset + $t_len - 1 : $t_end - $t_offset),
	$q_strand,
	($q_len == 0 ? $gap_string : substr($q_seq,$q_start + $q_offset -1,$q_len)),
	($t_len == 0 ? $gap_string : substr($t_seq,1 + $t_offset -1,$t_len))
      ];
      push(@pairs,$pair);
    }
    $q_offset += $q_len;
    $t_offset += $t_len;
  }
  
  $data{'pairs'} = \@pairs;
  
  return \%data;
}

sub seqname {

  my $self = shift;
  $self->{'seqname'} = shift if(@_);
  return $self->{seqname};

}

sub type {

  my $self = shift;
  $self->{'type'} = shift if(@_);
  return $self->{type};

}

sub identical_matches {

  my $self = shift;
  $self->{'identical_matches'} = shift if(@_);
  return $self->{'identical_matches'};

}

sub get_annotations {
  my $lrg_name = shift;
  my $lrg_coord_system_name = shift;
  
  # These parameters need only be defined if we need to add an entry to core db
  my $chr_name = shift;
  my $lrg_len = shift;
  my $pairs = shift;
  
  #ÊTry to fetch a slice for the LRG (from the read-write db), this should be possible if an entry exists in the core db
  my $sa_rw = $dbCore_rw->get_SliceAdaptor();
  my $lrg_slice = $sa_rw->fetch_by_region($lrg_coord_system_name,$lrg_name);
 
  # If it failed, insert mapping data to the core db
  if (!defined($lrg_slice)) {
  
    # Add a mapping between the LRG and chromosome to the core db
    # For consistency in the annotations, do transfer between LRG and chromosome coord systemsÊeven if there is a perfect match
    $LRGImport::dbCore = $dbCore_rw;
    LRGImport::add_mapping($lrg_name,$lrg_coord_system_name,$lrg_len,$current_assembly,$chr_name,$pairs);
    
    # Reload the db connections to get rid of any caching (is this necessary?)
    Bio::EnsEMBL::Registry->clear();
    Bio::EnsEMBL::Registry->load_all( $registry_file );
    $dbCore_rw = Bio::EnsEMBL::Registry->get_DBAdaptor('human','core_rw');
    
    # Now, try again to fetch the LRG slice
    $sa_rw = $dbCore_rw->get_SliceAdaptor();
    $lrg_slice = $sa_rw->fetch_by_region($lrg_coord_system_name,$lrg_name);
    
    die('Could not get LRG slice despite adding it to core db!') if (!defined($lrg_slice));
  }  
  
  my $genes = get_overlapping_genes($lrg_slice);
  my @gene_features;
  
  foreach my $gene (@{$genes}) {
    my $gene_feature = gene_2_feature($gene,$lrg_slice);
    push(@gene_features,$gene_feature);
  }
  
  return \@gene_features;
}

# Get the outer limits of a feature when projected onto another slice or undef if the projection failed
sub get_feature_limits {
  my $feature = shift;
  my $slice = shift;
  
  my %limits;
  
  #ÊGet the start position of the feature on the slice
  my $lrg_position = lift_coordinate($feature->start(),$feature->slice(),$slice);
  if (defined($lrg_position->{'position'})) {
    $limits{'start'} = $lrg_position->{'position'};
    $limits{'strand'} = $lrg_position->{'strand'};
  }
  else {
    #ÊIf start position lies upstream of LRG slice, indicate that the feature is partial in its 3'-end or 5'-end
    if ($lrg_position->{'reason'} eq 'upstream') {
      $limits{'start'} = 1;
      if ($feature->strand() > 0) {
	$limits{'partial_5'} = 1;
      }
      else {
	$limits{'partial_3'} = 1;
      }
    }
    # If start position lies downstream of LRG, the entire feature must be downstream, so indicate that
    elsif ($lrg_position->{'reason'} eq 'downstream') {
      $limits{'downstream'} = 1;
      $limits{'start'} = $slice->length() + 1;
    }
  }
  
  #ÊGet the end position of the feature on the slice
  $lrg_position = lift_coordinate($feature->end(),$feature->slice(),$slice);
  if (defined($lrg_position->{'position'})) {
    $limits{'end'} = $lrg_position->{'position'};
    $limits{'strand'} = $lrg_position->{'strand'};
  }
  else {
    #ÊIf end position lies downstream of LRG slice, indicate that the feature is partial in its 5'-end or 3'-end
    if ($lrg_position->{'reason'} eq 'downstream') {
      if ($feature->strand() > 0) {
	$limits{'partial_3'} = 1;
      }
      else {
	$limits{'partial_5'} = 1;
      }
      $limits{'end'} = $slice->length();
    }
    # If end position lies upstream of LRG, the entire feature must be upstream, so indicate that
    elsif ($lrg_position->{'reason'} eq 'upstream') {
      $limits{'upstream'} = 1;
      $limits{'end'} = -1;
    }
    
  }
  
  return \%limits;
}

#ÊLift a specific coordinate from a source slice to another slice. Returns a hash with the new coordinate and a field for storing the reason why mapping failed
sub lift_coordinate {
  my $position = shift;
  my $input_slice = shift;
  my $target_slice = shift;
  
  my %mapping;
  
  my $source_slice;
  #ÊIf the source slice is fetched from the read-only database, get the corresponding slice from the read-write db instead
  if ($input_slice->adaptor->dbc->host() eq $dbCore_rw->dbc->host() && $input_slice->adaptor->dbc->dbname() eq $dbCore_rw->dbc->dbname()) {
    $source_slice = $input_slice;
  }
  else {
    my $sa = $dbCore_rw->get_SliceAdaptor();
    $source_slice = $sa->fetch_by_region(
      $input_slice->coord_system_name(),
      $input_slice->seq_region_name(),
      $input_slice->start(),
      $input_slice->end(),
      $input_slice->strand()
    );
  }
  
  #ÊProject the feature onto the target slice
  my $projections = $source_slice->project_to_slice($target_slice);
  
  #ÊIf feature fails to project (lies entirely outside of the supplied slice), return with reason 'failed'
  if (!defined($projections) || scalar(@{$projections}) == 0) {
    $mapping{'reason'} = 'failed';
    return \%mapping;
  }
  
  # If position is 5' of what's been mapped
  $mapping{'reason'} = 'upstream' if ($position < $projections->[0][0]);
  #ÊIf position is 3' of what's been mapped
  $mapping{'reason'} = 'downstream' if ($position > $projections->[-1][1]);
  
  # Return if position was outside of lrg_slice
  return \%mapping if (defined($mapping{'reason'}));
  
  my $new_position = undef;
  my $new_strand = undef;
  # Go over the projected segments and find the coordinate liftover value
  foreach my $segment (@{$projections}) {
    
    # If position falls outside of the segment, skip to the next
    next if ($position < $segment->[0] || $position > $segment->[1]);
    
    # The offset of the coordinate within the segment
    my $offset = $position - $segment->[0];
    
    # The liftover position
    $new_position = $segment->[2]->start() + $offset;
    $new_strand = $segment->[2]->strand();
    
    last;
  }
  
  # If new coordinate still couldn't be determined, this must be because it is deleted in the target slice
  $mapping{'reason'} = 'deleted' if (!defined($new_position));
  
  $mapping{'position'} = $new_position;
  $mapping{'strand'} = $new_strand;
  
  return \%mapping;
}

sub gene_2_feature {

	my $gene = shift;
	my $lrg_slice = shift;
	my $slice_end = $lrg_slice->end;
	my ($current, $entry);
	
	# Check if the boundaries of the gene extends beyond the mapped region, in which case the end points should be set to those of the mapped region and a flag indicating partial overlap should be set
	my $limits = get_feature_limits($gene,$lrg_slice);
	
	my $partial_5 = $limits->{'partial_5'};
	my $partial_3 = $limits->{'partial_3'};
	
	
	my $feat_slice = $gene->feature_Slice();
	
	my $feat_start = $limits->{'start'};
	my $feat_end = $limits->{'end'};
	my $feat_strand = $limits->{'strand'};
	
#	if ($feat_strand < 0) {
#	  ($partial_5,$partial_3) = ($partial_3,$partial_5);
#	}
	
	# create the gene node
	my $gene_node = LRG::Node->new("gene", undef, {'symbol' => $gene->external_name, 'start' => $feat_start, 'end' => $feat_end, 'strand' => $gene->strand});
 
	# If the gene partially overlaps, this should be indicated
	if ($partial_5) {
	  $gene_node->addNode('partial')->content('5-prime');
	}
	if ($partial_3) {
	  $gene_node->addNode('partial')->content('3-prime');
	}
	
	# get xrefs for the gene
	my $entries = $gene->get_all_DBEntries();

	#print "\n\n>>> ", $gene->stable_id, " ", $gene->biotype, " ", $gene->status, "\n";
	my $long_name;
	while($entry = shift @$entries) {
		
		next unless $entry->dbname =~ /GI|RefSeq|HGNC$/;
		next if $entry->dbname =~ /RefSeq_peptide/;
		
		# get synonyms from HGNC entry
		if($entry->dbname eq 'HGNC') {
		        $gene_node->addNode('lrg_gene_name',{'source' => $entry->dbname})->content($entry->display_id);
			foreach my $synonym(@{$entry->get_all_synonyms}) {
				$gene_node->addNode('synonym')->content($synonym);
			}
			if (length($entry->description) > 1) {
			  $long_name = LRG::Node->new('long_name');
			  $long_name->content($entry->description);
			  $gene_node->addExisting($long_name) unless ($gene_node->nodeExists($long_name));
			}
		}
		$gene_node->addExisting(xref($entry));
		
		#print $entry->dbname, " ", $entry->primary_id, ".", $entry->version, " ", $entry->description, "\n";
	}
	# If no HGNC reference was found, add the description from Ensembl as long_name
	if (!defined($long_name) && defined($gene->description()) && length($gene->description()) > 0) {
	  $long_name = LRG::Node->new('long_name');
	  $long_name->content($gene->description());
	  $gene_node->addExisting($long_name);
	}
	
	$gene_node->addEmptyNode('db_xref', {'source' => 'Ensembl', 'accession' => $gene->stable_id});

	# get the xrefs for the transcript
	my %ext = ();
	my %extdesc = ();
	my $note_node;
	
	foreach my $trans(@{$gene->get_all_Transcripts}) {
		
# Check for partial overlaps
	  $limits = get_feature_limits($trans,$lrg_slice);
# If the transcript falls entirely outside of the LRG it should be skipped
	  next if (!defined($limits) || $limits->{'upstream'} || $limits->{'downstream'});
	  
	  $feat_start = $limits->{'start'};
	  $feat_end = $limits->{'end'};
	  $feat_strand = $limits->{'strand'};
	  $partial_5 = $limits->{'partial_5'};
	  $partial_3 = $limits->{'partial_3'}; 	
		
#ÊIf the transcript only partially overlaps the LRG, we need to set the end points of the transcript to the end point of the last exon that falls within the LRG record
	  my $exon_overlap = 1;
	  
	  if ($partial_5 || $partial_3) {
	    my $exons = $trans->get_all_Exons();
	    $exon_overlap = 0;
	    if (scalar(@{$exons}) > 0) {
	      $feat_start = $slice_end;
	      $feat_end = 1;
	      foreach my $exon (@{$exons}) {
#	        my $lrg_exon = $exon->transfer($lrg_slice);
#	        next if (!defined($lrg_exon));
	        $limits = get_feature_limits($exon,$lrg_slice);
	        next if (!defined($limits) || $limits->{'upstream'} || $limits->{'downstream'});
	        
# Set a flag to indicate that there is in fact an exon in this transcript within the scope of the LRG
		$exon_overlap = 1;
		if ($limits->{'end'} > 1) {
		  $feat_start = min($feat_start,$limits->{'start'});
		}
		if ($limits->{'start'} < $slice_end) {
		  $feat_end = max($feat_end,$limits->{'end'});
		}
	      }
	      if ($feat_start > $feat_end) {
	        ($feat_start,$feat_end) = ($feat_end,$feat_start);
	      }
	    }
	  }
		
	  my $cds_node = LRG::Node->new("transcript", undef, {'source' => 'Ensembl', 'start' => $feat_start, 'end' => $feat_end, 'transcript_id' => $trans->stable_id});
		
# If transcript partially overlaps, it should be indicated
	  if ($partial_5) {
	    $cds_node->addNode('partial')->content('5-prime');
	  }
	  if ($partial_3) {
	    $cds_node->addNode('partial')->content('3-prime');
	  }
	  
	  $entries = $trans->get_all_DBEntries();
	  
	  while($entry = shift @$entries) {
	  	
	    next unless $entry->dbname =~ /GI|RefSeq|MIM_GENE|Entrez|CCDS|RFAM|miRBase|pseudogene.org/;
	    next if $entry->dbname =~ /RefSeq_peptide/;
	    
	    if($entry->dbname eq 'RefSeq_dna') {
	      $note_node = LRG::Node->new('long_name');
	      $note_node->content($entry->description);
	      $cds_node->addExisting($note_node) unless ($cds_node->nodeExists($note_node));
	    }
	    
	    if ($entry->dbname !~ /MIM_GENE/) {
	      $cds_node->addExisting(xref($entry));
	    }
	    else {
	      $gene_node->addExisting(xref($entry));
	    }
	    
	    #print $entry->dbname, " ", $entry->primary_id, ".", $entry->version, " ", $entry->description, "\n";
	    
	    $ext{$entry->dbname} = $entry->primary_id;
	    $extdesc{$entry->dbname} = $entry->description;
	  }
	  
	  my $protein = $trans->translation;
	  
	  if($protein && length($protein) > 0) {
	    my $cds_start;
	    my $cds_end;
	    undef($partial_5);
	    undef($partial_3);
	    my $lift = lift_coordinate($trans->coding_region_start,$trans->slice(),$lrg_slice);
	    # In case the position couldn't be found on the LRG slice, investigate why
	    if (!defined($lift->{'position'})) {
	      #ÊIf coding start lies downstream of LRG slice, the protein should be skipped
	      if ($lift->{'reason'} eq 'downstream') {
		$cds_start = $slice_end + 1;
	      }
	      #ÊIf coding start lies upstream of LRG slice, set start of protein to the closest exon start
	      elsif ($lift->{'reason'} eq 'upstream') {
		$cds_start = $feat_start;
		$partial_3 = 1;
	      }
	      #ÊIf it is because of a deletion, don't know what to do... (set to closest, non-deleted position?)
	      elsif ($lift->{'reason'} eq 'deleted') {}
	      #ÊElse if it failed to project, this is likely to be because the transcript is entirely outside of LRG. Shouldn't have gotten this far in that case but still, skip protein in that case
	      elsif ($lift->{'reason'} eq 'failed') {
		$cds_start = $slice_end + 1;
	      }
	    }
	    else {
	      $cds_start = $lift->{'position'};
	    }
	    
	    $lift = lift_coordinate($trans->coding_region_end,$trans->slice(),$lrg_slice);
	    # In case the position couldn't be found on the LRG slice, investigate why
	    if (!defined($lift->{'position'})) {
	      #ÊIf coding end lies upstream of LRG slice, the protein should be skipped
	      if ($lift->{'reason'} eq 'upstream') {
		$cds_end = -1;
	      }
	      #ÊIf coding end lies downstream of LRG slice, set end of protein to the closest exon end
	      elsif ($lift->{'reason'} eq 'downstream') {
		$cds_end = $feat_end;
		$partial_5 = 1;
	      }
	      #ÊIf it is because of a deletion, don't know what to do... (set to closest, non-deleted position?)
	      elsif ($lift->{'reason'} eq 'deleted') {}
	      #ÊElse if it failed to project, this is likely to be because the transcript is entirely outside of LRG. Shouldn't have gotten this far in that case but still, skip protein in that case
	      elsif ($lift->{'reason'} eq 'failed') {
		$cds_end = -1;
	      }
	    }
	    else {
	      $cds_end = $lift->{'position'};
	    }

# Some gene-level data (e.g. NCBI GeneID) is buried in the translation object. So before determining if we skip the transcript, parse the xrefs
	    my $prot_node = LRG::Node->new('protein_product');
	    $prot_node->data(
	      {
	        'source' => 'Ensembl',
	        'accession' => $protein->stable_id,
	        'cds_start' => $cds_start,
	        'cds_end' => $cds_end
	      }
	    );
	    	    
	    $entries = $protein->get_all_DBEntries();
	      
	    while($entry = shift @$entries) {
	    	  
	      next unless $entry->dbname =~ /RefSeq|Uniprot|CCDS|MIM_GENE|Entrez/;
	    	  
	      if($entry->dbname eq 'RefSeq_peptide') {
	        $note_node = LRG::Node->new('long_name');
	        $note_node->content($entry->description);
	        $prot_node->addExisting($note_node) unless ($prot_node->nodeExists($note_node));
	      }
	      
	      if ($entry->dbname !~ /MIM_GENE|Entrez/) {
	        $prot_node->addExisting(xref($entry));
	      }
	      else {
	        my $xref = xref($entry);
	        $gene_node->addExisting($xref) unless ($gene_node->nodeExists($xref));
	      }
	    }

# If the transcript falls entirely outside of the LRG it should be skipped
	    unless ($cds_end < 1 || $cds_start > $slice_end) {
	      
	      if ($partial_5) {
		$prot_node->addNode('partial')->content('5-prime');
	      }
	      if ($partial_3) {
		$prot_node->addNode('partial')->content('3-prime');
	      }
	      
	      $cds_node->addExisting($prot_node);
		      
	      #print $entry->dbname, " ", $entry->primary_id, ".", $entry->version, " ", $entry->description, "\n";
	    }
	  }
	  $gene_node->addExisting($cds_node) unless ($gene->biotype eq 'protein_coding' && (!$protein || !$exon_overlap));
	}
	
	# finish the gene with xrefs
	#$gene_node->addEmptyNode('db_xref', {'source' => 'GeneID', 'accession' => $ext{'EntrezGene'}}) if defined $ext{'EntrezGene'};
	#$gene_node->addEmptyNode('db_xref', {'source' => 'HGNC', 'accession' => $ext{'HGNC'}}) if defined $ext{'HGNC'};
	#$gene_node->addEmptyNode('db_xref', {'source' => 'MIM', 'accession' => $ext{'MIM_GENE'}}) if defined $ext{'MIM_GENE'};
	
	return $gene_node;
}

sub xref {
	my $entry = shift;
	
	my $db = $entry->dbname;
	
	$db =~ s/EntrezGene.*/GeneID/;
	$db =~ s/Uniprot.*/UniProtKB/;
	$db =~ s/RefSeq.*/RefSeq/;
	$db =~ s/MIM.*/MIM/;
	
	my $id = $entry->primary_id.($entry->version > 0 ? ".".$entry->version : "");
	
	my $node = LRG::Node->new('db_xref', undef, {'source' => $db, 'accession' => $id});
	$node->empty(1);
	
	return $node;
}

#ÊCreate an array of FeaturePair objects from a mapping tag, one for each mapping span
sub mapping_2_feature_pair {
	my $mapping = shift;
	my $slice_adaptor = shift;
	
	my $hseqname = $mapping->data->{'chr_name'};
	my $slice = $slice_adaptor->fetch_by_region('chromosome',$hseqname);
	my @feature_pairs;
	
	foreach my $span (@{$mapping->findNodeArray('mapping_span')}) {
		my $hstart = $span->data->{'start'};
		my $hend = $span->data->{'end'};
		my $hstrand = 1;
		my $start = $span->data->{'lrg_start'};
		my $end = $span->data->{'lrg_end'};
		my $strand = $span->data->{'strand'};
		my $identical_matches = 1;
		
		my @type;
		push (@type,['DNA',$start,$end,$hstart,$hend,$strand]);
		my ($s,$e,$hs,$he);
		my $last_end = 0;
		my $last_hend = ($strand > 0 ? ($hstart-1) : ($hend+1));
		
		my $diffs = $span->findNodeArray('diff');
		if (defined($diffs)) {
			foreach my $diff (@{$diffs}) {
				$s = $last_end+1;
				$e = $diff->data->{'lrg_start'}-1;
				$hs = ($strand > 0 ? ($last_hend+1) : ($diff->data->{'end'}+1));
				$he = ($strand > 0 ? ($diff->data->{'start'}-1) : ($last_hend-1));
				push(@type,['M',$s,$e,$hs,$he,$strand]);
				
				$s = $diff->data->{'lrg_start'};
				$e = $diff->data->{'lrg_end'};
				$hs = $diff->data->{'start'};
				$he = $diff->data->{'end'};
				my $seq = $diff->data->{'lrg_sequence'};
				my $hseq = $diff->data->{'genomic_sequence'};
				push(@type,['M',$s,$e,$hs,$he,$strand,$seq,$hseq]);
				
				$last_end = $e;
				$last_hend = ($strand > 0 ? $he : $hs);
				$identical_matches = 0;
			}
		}
		$s = $last_end+1;
		$e = $end;
		$hs = ($strand > 0 ? ($last_hend+1) : $hstart);
		$he = ($strand > 0 ? $hend : ($last_hend-1));
		push(@type,['M',$s,$e,$hs,$he,$strand]);

		my $fp = Bio::EnsEMBL::FeaturePair->new(
			-start    => $start,
			-end      => $end,
			-strand   => $strand,
			-display_id => $hseqname,
			-hstart   => $hstart,
			-hend     => $hend,
			-hstrand  => $hstrand,
			-hseqname => $hseqname,
			-slice    => $slice,
		);
		$fp->identical_matches($identical_matches);
		$fp->type(\@type);
		
		push(@feature_pairs,$fp);		
	}
	return \@feature_pairs;
}

# Convert an XML mapping structure into a pair array. Needs all sequences to be in forward direction (will flip them on its own)
sub mapping_2_pairs {
  my $mapping_node = shift;
  my $lrg_seq = shift;
  my $chr_seq = shift;
  
  # Loop over the mapping spans and create the pairs array
  my @pairs;
  my $mapping_spans = $mapping_node->findNodeArray('mapping_span');
  my $last_lrg_start;
  my $last_lrg_end;
  my $last_chr_start;
  my $last_chr_end;
  my $last_strand = 1;
  my $chr_name = $mapping_node->data->{'chr_name'};
  my $chr_offset_start = $mapping_node->data->{'chr_start'};
  my $chr_offset_end = $mapping_node->data->{'chr_end'};
  my $lrg_offset_start = $mapping_node->data->{'lrg_start'};
  my $lrg_offset_end = $mapping_node->data->{'lrg_end'};
  
  foreach my $span (@{$mapping_spans}) {
    my $lrg_start = $span->data->{'lrg_start'};
    my $lrg_end = $span->data->{'lrg_end'};
    my $chr_start = $span->data->{'start'};
    my $chr_end = $span->data->{'end'};
    my $strand = $span->data->{'strand'};
    
    # If we have changed orientation since last span (or since start when sequence was forward), flip the sequence
    if ($strand != $last_strand) {
      reverse_comp(\$chr_seq);
      $last_strand = $strand;
    }
    
    #ÊTwo mapping spans are separated by a gap, so check if we need to add that
    if (defined($last_lrg_start)) {
      my $lrgstring;
      my $chrstring;
      my $gaplen = max($lrg_start - $last_lrg_end,($strand > 0 ? $chr_start - $last_chr_end : $last_chr_start - $chr_end)) - 1;
      # If the LRG has an insertion relative to the reference, get the string
      if ($lrg_start - $last_lrg_end > 1 && defined($lrg_seq)) {
	$lrgstring = substr($lrg_seq,$last_lrg_end,$lrg_start - $last_lrg_end - 1);
      }
      # Else, if the LRG has a deletion, create a gap string of the same length as the deletion
      elsif ($lrg_start - $last_lrg_end == 1) {
	$lrgstring = '-' x $gaplen;
      }
      
      # If the reference has an insertion
      if (
	  defined($chr_seq) && (
	    ($strand > 0 && ($chr_start - $last_chr_end) > 1) ||
	    ($strand < 0 && ($last_chr_start - $chr_end) > 1)
	  )
	 ) {
	if ($strand > 0) {
	  $chrstring = substr($chr_seq,$last_chr_end-$chr_offset_start+1,$chr_start - $last_chr_end - 1);
	}
	else {
	  $chrstring = substr($chr_seq,$chr_offset_end-$chr_end-1,$last_chr_start - $chr_end - 1);
	}
      }
      #ÊElse, create a gap string
      elsif (
	      ($strand > 0 && ($chr_start - $last_chr_end) == 1) ||
	      ($strand < 0 && ($last_chr_start - $chr_end) == 1)
	    ) {
	$chrstring = '-' x $gaplen;
      }
      
      my $gap_pair = [
	'G',
	$last_lrg_end + 1,
	$lrg_start - 1,
	($strand > 0 ? $last_chr_end + 1 : $chr_end + 1),
	($strand > 0 ? $chr_start - 1 : $last_chr_start - 1),
	$strand,
	$lrgstring,
	$chrstring
      ];
      
      push(@pairs,$gap_pair);
    }
    
    my $pair = [
      'DNA',
      $lrg_start,
      $lrg_end,
      $chr_start,
      $chr_end,
      $strand
    ];
    push(@pairs,$pair);
    
    # Store these coordinates in the last_* variables
    $last_lrg_start = $lrg_start;
    $last_lrg_end = $lrg_end;
    $last_chr_start = $chr_start;
    $last_chr_end = $chr_end;
    
    # Loop over any mismatch information and create pairs
    my $diffs = $span->findNodeArray('diff');
    foreach my $diff (@{$diffs}) {
      # So far, can only handle mismatches, die with an error if something else is encountered
      die("Can not handle diffs of type '" . $diff->data->{'type'} . "'!") unless ($diff->data->{'type'} eq 'mismatch');
      
      $lrg_start = $diff->data->{'lrg_start'};
      $lrg_end = $diff->data->{'lrg_end'};
      $chr_start = $diff->data->{'start'};
      $chr_end = $diff->data->{'end'};
      my $lrg_string = $diff->data->{'lrg_sequence'};
      my $chr_string = $diff->data->{'genomic_sequence'};
      
      $pair = [
	'M',
	$lrg_start,
	$lrg_end,
	$chr_start,
	$chr_end,
	$strand,
	$lrg_string,
	$chr_string
      ];
      push(@pairs,$pair);
    }
  }
  
  my %mapping = (
    'lrg_start' => $lrg_offset_start,
    'lrg_end' => $lrg_offset_end,
    'chr_start' => $chr_offset_start,
    'chr_end' => $chr_offset_end,
    'chr_name' => $chr_name,
    'strand' => $last_strand,
    'pairs' => \@pairs
  );
  return \%mapping;
}

# Convert a pair array structure into a XML structure
sub pairs_2_mapping {
  my $pairs = shift;
  my $chr_assembly = shift;
  my $chr_name = shift;
  my $most_recent = shift;
  my $chr_id = shift;
  
  # Create the top mapping node
  my $mapping_node = LRG::Node->new('mapping');
  $mapping_node->addData({'assembly' => $chr_assembly,'chr_name' => $chr_name});
  $mapping_node->addData({'chr_id' => $chr_id}) if (defined($chr_id));
  $mapping_node->addData({'most_recent' => $most_recent}) if (defined($most_recent));
  
  #ÊOrder the pairs array according to ascending LRG coordinates
  my @sorted_pairs = sort {$a->[1] <=> $b->[1]} @{$pairs};
  
  # Separate the DNA chunks from mismatches in the pairs array. We will create a mapping span for each DNA chunk and include all mismatch information in this span
  # Gaps can be discarded at this point since they are implicit between DNA chunks
  my @dna_chunks;
  my @mismatches;
  while (my $pair = shift(@sorted_pairs)) {
    push(@dna_chunks,$pair) if ($pair->[0] eq 'DNA');
    push(@mismatches,$pair) if ($pair->[0] eq 'M');
  }
  
  # Loop over the DNA chunk array and create mapping spans and diffs as necessary 
  my $chr_map_start = 1e11;
  my $chr_map_end = -1;
  while (my $chunk = shift(@dna_chunks)) {
    my $lrg_start = $chunk->[1];
    my $lrg_end = $chunk->[2];
    my $chr_start = $chunk->[3];
    my $chr_end = $chunk->[4];
    my $strand = $chunk->[5];
    
    my $span_node = LRG::Node->new('mapping_span');
    $span_node->addData(
      {
	'start' => $chr_start,
	'end' => $chr_end,
	'lrg_start' => $lrg_start,
	'lrg_end' => $lrg_end,
	'strand' => $strand
      }
    );
    
    # Loop over the mismatches and create diffs as necessary
    foreach my $mismatch (@mismatches) {
      next unless ($mismatch->[1] >= $lrg_start && $mismatch->[2] <= $lrg_end);
      
      my $diff_node = LRG::Node->newEmpty('diff');
      $diff_node->addData(
	{
	  'type' => 'mismatch',
	  'lrg_start' => $mismatch->[1],
	  'lrg_end' => $mismatch->[2],
	  'start' => $mismatch->[3],
	  'end' => $mismatch->[4],
	  'lrg_sequence' => $mismatch->[6],
	  'genomic_sequence' => $mismatch->[7]
	}
      );
      
      $span_node->addExisting($diff_node);
    }
    $chr_map_start = min($chr_map_start,$chr_start);
    $chr_map_end = max($chr_map_end,$chr_end);
    
    $mapping_node->addExisting($span_node);
  }
  # Lastly, add the mapping chromosomal start and end points as attributes to the mapping node
  $mapping_node->addData({'chr_start' => $chr_map_start, 'chr_end' => $chr_map_end});
  
  return $mapping_node;
}

# Convert a pair array structure into an array of FeaturePairs
sub pairs_2_feature_pairs {
  my $pairs = shift;
  my $lrg_id = shift;
  my $chr_id = shift;
  my $chr_slice = shift;
  
  # Separate the DNA chunks from mismatch annotations in the pairs array. We will create a FeaturePair for each DNA chunk and include all mismatch information in this span
  # Gaps can be discarded at this point since they are implicit between DNA chunks
  my @dna_chunks;
  my @mismatches;
  while (my $pair = shift(@{$pairs})) {
    push(@dna_chunks,$pair) if ($pair->[0] eq 'DNA');
    push(@mismatches,$pair) if ($pair->[0] eq 'M');
  }
  
  # Create a FeaturePair for each DNA chunk
  my @feature_pairs;
  foreach my $chunk (@dna_chunks) {
    my $lrg_start = $chunk->[1];
    my $lrg_end = $chunk->[2];
    my $chr_start = $chunk->[3];
    my $chr_end = $chunk->[4];
    my $strand = $chunk->[5];
    my @type;
    
    #ÊGet all the mismatch information that is spanned by this chunk
    foreach my $mismatch (@mismatches) {
      next unless($mismatch->[1] >= $lrg_start && $mismatch->[2] <= $lrg_end);
      push(@type,$mismatch);
    }

    my $fp = Bio::EnsEMBL::FeaturePair->new(
      -start    => $lrg_start,
      -end      => $lrg_end,
      -strand   => 1,
      -display_id => $chr_id,
      -hstart   => $chr_start,
      -hend     => $chr_end,
      -hstrand  => $strand,
      -hseqname => $chr_id,
      -slice    => $chr_slice
    );
    $fp->identical_matches((scalar(@type) == 0));
    $fp->type(\@type);
    
    push(@feature_pairs,$fp);		    
  }
  
  return \@feature_pairs;
}

sub get_overlapping_genes {
  my $lrg_slice = shift;
  
  # Project the LRG slice to the chromosome coordinate system and get the extreme points in that coord system
  my $segments = $lrg_slice->project('chromosome');
  my $chr_name;
  my $chr_start = 1e11;
  my $chr_end = -1;
  my $strand;
  while (my $segment = shift(@{$segments})) {
    $chr_start = min($chr_start,$segment->[2]->{'start'});
    $chr_end = max($chr_end,$segment->[2]->{'end'});
    # Check that one or more segments don't map to a different chromosome
    die('The mapping spans are on different chromosomes!') if (defined($chr_name) && $chr_name ne $segment->[2]->{'seq_region_name'});
    $chr_name = $segment->[2]->{'seq_region_name'};
    # Do an extra check to make sure that one or more of the segments are not flipped (this could be valid but then we must modify the scripts to deal with it)
    die('One or more mapping spans are flipped w.r.t. the others!') if (defined($strand) && $strand != $segment->[2]->{'strand'});
    $strand = $segment->[2]->{'strand'};
  }
  
  # Fetch a chromosome slice spanning the LRG from the reference (read-only db)
  my $sa_ro = $dbCore_ro->get_SliceAdaptor();
  my $ref_slice = $sa_ro->fetch_by_region("chromosome",$chr_name,$chr_start,$chr_end,$strand);
  my $genes = $ref_slice->get_all_Genes();
  
  return $genes;
}


1;
