use strict;
use warnings;

package LRGMapping;

use Bio::EnsEMBL::LRGSlice;
use Bio::EnsEMBL::MappedSliceContainer;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use LRG::LRG;
use LRG::LRGImport;
use List::Util qw (min max);


# global variables

our $MIN_GAP_BETWEEN_SPANS = 1000;

our $SSAHA_BIN = 'ssaha2';

our $input_dir;
our $output_dir;

our $target_dir;
our %rec_seq;

#ÊKeep different db connections for a core db with write access to and one with read-only access
our $dbCore_rw;
our $dbCore_ro;
our $dbFuncGen;

our $lrg_name;
our $registry_file;
our $current_assembly;

our $mapping_num = 1;

our $SSAHA_RUNTIME_LIMIT = 5;
our $FARM_MEMORY = 4000;
our %SSAHA_PARAMETERS = (
  '-align' 	=> '1',
  '-kmer' 	=> '12',
  '-seeds' 	=> '25',
  '-cut' 	=> '1000',
  '-output' 	=> 'vulgar',
  '-depth'	=> '5',
  '-best' 	=> '1',
  '-memory'	=> '4000'
);

our $EXONERATE_BIN = 'exonerate';
our %EXONERATE_PARAMETERS = (
  '--model'		=>	'affine:global',
  '--exhaustive'	=>	'1',
  '--showvulgar'	=>	'1',
  '--showcigar'		=>	'0',
  '--showsugar'		=>	'0',
  '--showalignment'	=>	'1',
  '--bestn'		=>	'1'
);

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
  	
  my $queue = "-q normal -R'select[mem>" . $FARM_MEMORY . "] rusage[mem=" . $FARM_MEMORY . "]' -M" . $FARM_MEMORY . "000";
  	
  my $output_file_name = "ssaha2_output\_$input_file_name";
  my $input_file = "$input_dir/$input_file_name";
  my $output_file ="$output_dir/$output_file_name";
  my $wait_output = "$input_dir/$$.waiting.out";
  
  my ($subject, %rec_find, %input_length, %done);
  	
  $subject = "$target_dir/ref";
  $rec_seq{$name} = $sequence;
  my $seqobj = Bio::PrimarySeq->new(-id => $name, -seq => $rec_seq{$name});
  $rec_seq{$name} = $seqobj;
  	
  my $error_file = bsub_ssaha_job($queue,$input_file,$output_file,$subject);
  
  my $call = "bsub -o $wait_output -K -w 'done($input_file\_ssaha_job)' -J waiting_process -c $SSAHA_RUNTIME_LIMIT sleep 1"; #waits until all ssaha jobs have finished to continue (but maximum $SSAHA_RUNTIME_LIMIT minutes)
  system($call);
  
  # Check that the waiting process wasn't killed by the runtime limit
  die("*** ERROR *** SSAHA2 mapping process didn't finish in time. This is probably because ssaha2 crashed.") if (`grep 'TERM_CPULIMIT' $wait_output`);
  
  # Check the error file
  my @error;
  open(IN,"<",$error_file) or die ("Could not read from error file $error_file");
  while (<IN>) {
    chomp;
    push(@error,$_);
  }
  close(IN);
  
  if (scalar(@error)) {
    print STDERR "*** ERROR ***\nssaha2 generated the following error:\n";
    print STDERR "\t" . join("\n\t",@error) . "\n";
    print STDERR "*************\n";
  }
  
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
  
  unlink($output_file);
  unlink($input_file);
  unlink($error_file);
  unlink($wait_output);
  
  return $mapping;
}

# subroutine to submit a job to the farm
sub bsub_ssaha_job {
  my ($queue, $input_file, $output_file, $subject) = @_;
  
  my $host = `hostname`;
  chomp $host;
  
  my $error_file = "$output_dir/error_ssaha";
  unlink($error_file);
  
  my $ssaha_command = "$SSAHA_BIN";
  while (my ($param,$val) = each(%SSAHA_PARAMETERS)) {
    $ssaha_command .= " $param $val";
  }
  $ssaha_command .= " -save $subject $input_file";
#  my $call = "echo '$ssaha_command; scp $output_file $host:$input_dir/' | bsub -E 'scp $host:$input_file $input_dir/' -J $input_file\_ssaha_job -P ensembl-variation $queue -e $output_dir/error_ssaha -o $output_file";
  my $call = "echo '$ssaha_command' | bsub -J $input_file\_ssaha_job $queue -e $error_file -o $output_file";
	
  system ($call);
#  print $call, "\n";
  return $error_file;
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

}

# Align two input sequences using exonerate. Does not use LSF as this method is only expected to be used for cDNA (=> quick)
sub exonerate_align {
  my $q_seq = shift;
  my $t_seq = shift;
  
  # Temp files in the temp directory for input and output files to exonerate
  my $q_file = $output_dir . "/query_" . time() . ".fa"; 
  my $t_file = $output_dir . "/target_" . time() . ".fa";
  my $o_file = $output_dir . "/exonerate_" . time() . ".out";
  
  # The command that will be executed
  my $cmd = $EXONERATE_BIN . " --query $q_file --target $t_file";
  while (my ($key,$val) = each(%EXONERATE_PARAMETERS)) {
    $cmd .= " $key $val";
  }
  
  # Write the input sequences to the temp files in fasta format
  open(OUT,'>',$q_file) or die("Could not write input sequence to temporary file $q_file");
  print OUT ">query\n$q_seq\n";
  close(OUT);
  open(OUT,'>',$t_file) or die("Could not write input sequence to temporary file $t_file");
  print OUT ">target\n$t_seq\n";
  close(OUT);
  
  # Open the output file for writing
  open(OUT,'>',$o_file) or die("Could not open temporary exonerate output file $o_file");
  
  # Execute exonerate and write the output to the tempfile
  print OUT `$cmd`;
  close(OUT);
  
  return $o_file;
}

#ÊParse the exonerate output
sub parse_exonerate {
  my $output_file = shift;
  
  # Open the output file for reading
  open(IN,'<',$output_file) or die("Could not open exonerate output file $output_file for reading");
  
  my $q_seq = "";
  my $q_start;
  my $q_end;
  my $t_seq = "";
  my $t_start;
  my $t_end;
  my $vulgar = "";
  my $lines = 0;
  
  # Parse the output and concatenate the aligned sequences. Also store the ranges
  while (<IN>) {
    chomp;
    if ($_ =~ m/^\s*\d+\s+\:\s+(\S+)\s+\:\s+\d+\s*$/) {
      $lines++;
      # Odd alignment lines correspond to the query, even to the target
      if ($lines%2 == 0) {
	$t_seq .= $1;
      }
      else {
	$q_seq .= $1;
      }
    }
    elsif ($_ =~ m/Query range\:\s+(\d+)\s+\-\>\s+(\d+)/) {
      $q_start = $1;  
      $q_end = $2;  
    }
    elsif ($_ =~ m/Target range\:\s+(\d+)\s+\-\>\s+(\d+)/) {
      $t_start = $1;  
      $t_end = $2;  
    }
    elsif ($_ =~ m/^\s*(vulgar\:.+)$/) {
      $vulgar = $1;
      last;
    }
  }
  close(IN);
  
  my $align_length = length($q_seq);
  
  # Remove gap characters from the sequences (i.e. revert to unaligned sequences)
  $q_seq =~ s/\-//g;
  $t_seq =~ s/\-//g;
  
  # The parsing requires 1-based indices. The exonerate output is 0-based, so increment the start offsets
  $vulgar =~ s/^(vulgar\:\s+[^\s]+\s+)(\d+)(\s+\d+\s+\+\s+[^\s]+\s+)(\d+)(\s+\d+\s+\+\s+\d+\s+.+)$/$1@{[$2+1]}$3@{[$4+1]}$5/;
  
  # Use the method for parsing the vulgar string to get a pairs hash
  my $data = parse_vulgar_string($vulgar,$q_seq,$t_seq);
  
  #ÊAdd the alignment length to the data hash
  $data->{'length'} = $align_length;
  
  return $data;
}

# 
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
# The input sequences should be unaligned
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
  my ($t_name) = $t_id =~ m/^([^\-]+)/;
  
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
	if (substr($q_seq,1 + $q_offset + $i - 1,1) ne substr($t_seq,1 + $t_offset + $i - 1,1)) {
	  $pair = [
	    $lbl,
	    $q_start + $q_offset + $i,
	    $q_start + $q_offset + $i,
	    ($t_strand > 0 ? $t_start + $t_offset + $i : $t_end - $t_offset - $i),
	    ($t_strand > 0 ? $t_start + $t_offset + $i : $t_end - $t_offset - $i),
	    $q_strand,
	    substr($q_seq,1 + $q_offset + $i - 1,1),
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
	($q_len == 0 ? '' : substr($q_seq,1 + $q_offset -1,$q_len)),
	($t_len == 0 ? '' : substr($t_seq,1 + $t_offset -1,$t_len))
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

sub clear_mapping {
  my $lrg_name = shift;
  my $lrg_coord_system_name = shift;
  
  $LRGImport::dbCore = $dbCore_rw;
  LRGImport::purge_db($lrg_name,$lrg_coord_system_name);
}

sub get_annotations {
  my $lrg_name = shift;
  my $lrg_coord_system_name = shift;
  
  # If not specified, use default value 'lrg' for coord_system_name
  $lrg_coord_system_name ||= 'lrg';
  
  # These parameters need only be defined if we need to add an entry to core db
  my $chr_name = shift;
  my $lrg_len = shift;
  my $mapping = shift;
  my $pairs = $mapping->{'pairs'};
  
  #ÊTry to fetch a slice for the LRG (from the read-write db), this should be possible if an entry exists in the core db
  my $sa_rw = $dbCore_rw->get_SliceAdaptor();
  # Since this method dies if not successful, wrap it in an eval block
  my $lrg_slice;
  eval {
    $lrg_slice = $sa_rw->fetch_by_region($lrg_coord_system_name,$lrg_name);
  };
  if ($@) {
    print "Could not create LRG slice, will add mapping to db\n";
  }
 
  # If it failed, insert mapping data to the core db. 
  if (!defined($lrg_slice)) {
  
    # Return undef if the parameters needed for mapping were not supplied
    return undef unless (defined($chr_name) && defined($lrg_len) && defined($mapping));
    
    # Add a mapping between the LRG and chromosome to the core db
    # For consistency in the annotations, do transfer between LRG and chromosome coord systemsÊeven if there is a perfect match
    $LRGImport::dbCore = $dbCore_rw;
    LRGImport::add_mapping($lrg_name,$lrg_coord_system_name,$lrg_len,$mapping);
    
    # Reload the db connections to get rid of any caching (is this necessary?)
    Bio::EnsEMBL::Registry->clear();
    Bio::EnsEMBL::Registry->load_all( $registry_file );
    $dbCore_rw = Bio::EnsEMBL::Registry->get_DBAdaptor('human','core_rw');
    
    # Now, try again to fetch the LRG slice
    $sa_rw = $dbCore_rw->get_SliceAdaptor();
    $lrg_slice = $sa_rw->fetch_by_region($lrg_coord_system_name,$lrg_name);
    
    die('Could not get LRG slice despite adding it to core db!') if (!defined($lrg_slice));
  }  
  
  my $genes = get_overlapping_features($lrg_slice,'gene');
  my @gene_features;
  
  foreach my $gene (@{$genes}) {
    my $gene_feature = gene_2_feature($gene,$lrg_slice);
    push(@gene_features,$gene_feature);
  }
  
=head Will not go live with this yet
  my $reg_features = get_overlapping_features($lrg_slice,'regulatory');
  foreach my $reg_feature (@{$reg_features}) {
    my $node = regulatory_feature_2_feature($reg_feature,$lrg_slice);
    push(@gene_features,$node);
  }
=cut

  return \@gene_features;
}

# Get the outer limits of a feature when projected onto another slice or undef if the projection failed
sub get_feature_limits {
  my $feature = shift;
  my $slice = shift;
  
  my %limits;
  
  my $feat_start = $feature->start();
  my $feat_end = $feature->end();
  my $feat_slice = $feature->slice();
  my $feat_strand = $feature->strand();
  
  #ÊGet the start position of the feature on the slice
  my $lrg_position = lift_coordinate($feat_start,$feat_slice,$slice);
  if (defined($lrg_position->{'position'})) {
    $limits{'start'} = $lrg_position->{'position'};
    $limits{'strand'} = $lrg_position->{'strand'};
  }
  else {
    #ÊIf start position lies upstream of LRG slice, indicate that the feature is partial in its 3'-end or 5'-end
    if ($lrg_position->{'reason'} eq 'upstream') {
      $limits{'start'} = 1;
      if ($feat_strand > 0) {
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
  $lrg_position = lift_coordinate($feat_end,$feat_slice,$slice);
  if (defined($lrg_position->{'position'})) {
    $limits{'end'} = $lrg_position->{'position'};
    $limits{'strand'} = $lrg_position->{'strand'};
  }
  else {
    #ÊIf end position lies downstream of LRG slice, indicate that the feature is partial in its 5'-end or 3'-end
    if ($lrg_position->{'reason'} eq 'downstream') {
      if ($feat_strand > 0) {
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
  
  return {'reason' => 'failed'} if (!defined($position));
  
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
    warn("Failed to do a project_to_slice from " . $source_slice->coord_system_name() . ":" . $source_slice->seq_region_name() . ":" . $source_slice->start() . "-" . $source_slice->end() . ":" . $source_slice->strand() . " to " . $target_slice->coord_system_name() . ":" . $target_slice->seq_region_name() . ":" . $target_slice->start() . "-" . $target_slice->end() . ":" . $target_slice->strand() . ", will try again...");
    $projections = $source_slice->project_to_slice($target_slice);
    if (!defined($projections) || scalar(@{$projections}) == 0) {
      warn("Still failed!");
      $mapping{'reason'} = 'failed';
      return \%mapping;
    }
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

sub attach_protein {
  my $transcript = shift;
  my $lrg_slice = shift;
  my $transcript_node = shift;
  
  # Get the gene node that is the parent of the transcript node
  my $gene_node = $transcript_node->parent();
  
  # Get the translation of this transcript
  my $protein = $transcript->translation();
  
  if (!defined($protein) || length($protein->seq()) == 0) {
    return;
  }
  
  my $cds_start;
  my $cds_end;
  my $partial_5;
  my $partial_3;
  my $slice_end = $lrg_slice->end();
  # Flag to indicate if protein coding region is entirely outside of the LRG
  my $outside = 0;
  
  # Get the coding region start coordinates in LRG coordinates.
  #ÊNote that this is always the lowest coordinate, so the coding_region_start on a transcript is the stop codon position if the transcript is on the negative strand
  my $lift = lift_coordinate($transcript->coding_region_start(),$transcript->slice(),$lrg_slice);
  
  # In case the position couldn't be found on the LRG slice, investigate why
  if (!defined($lift->{'position'})) {
    
    #ÊIf coding start lies downstream of LRG slice, just set the start to the first position outside of the slice (no coding part is within the LRG)
    if ($lift->{'reason'} eq 'downstream') {
      $cds_start = $slice_end + 1;
      $outside = 1;
    }
    
    #ÊIf coding start lies upstream of LRG slice, set start of protein to the closest exon start
    elsif ($lift->{'reason'} eq 'upstream') {
      $partial_5 = 1;
    }
    
    #ÊIf it is because of a deletion, don't know what to do... (set to closest, non-deleted position?)
    elsif ($lift->{'reason'} eq 'deleted') {}
    
    #ÊElse if it failed to project, this is likely to be because the transcript is entirely outside of LRG. Shouldn't have gotten this far in that case but still, sset start outside of LRG slice
    elsif ($lift->{'reason'} eq 'failed') {
      $cds_start = $slice_end + 1;
    }
  }
  else {
    $cds_start = $lift->{'position'};
  }
  
  # Get the coding region end coordinates in LRG coordinates  
  $lift = lift_coordinate($transcript->coding_region_end(),$transcript->slice(),$lrg_slice);
  
  # In case the position couldn't be found on the LRG slice, investigate why
  if (!defined($lift->{'position'})) {
    #ÊIf coding end lies upstream of LRG slice, set the end to be negative
    if ($lift->{'reason'} eq 'upstream') {
      $cds_end = -1;
      $outside = 1;
    }
    
    #ÊIf coding end lies downstream of LRG slice, set end of protein to the closest exon end
    elsif ($lift->{'reason'} eq 'downstream') {
      $partial_3 = 1;
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
  
  # Create a node for the protein 
  my $protein_node = LRG::Node->new('protein_product');

  my $name_content = "";
  
  # Get all DBEntries
  my $entries = $protein->get_all_DBEntries();
  while(my $entry = shift(@{$entries})) {
    
    # We are only interested in a subset of the possible xrefs	  
    next unless $entry->dbname =~ /RefSeq|Uniprot|CCDS|MIM_GENE|Entrez/;
    	  
    if($entry->dbname eq 'RefSeq_peptide' && defined($entry->description())) {
      $name_content = $entry->description();
    }
    
    # Create an xref node
    my $xref = xref($entry);
    
    #ÊAttach the xref node to the protein if applicable 
    if ($entry->dbname() !~ /MIM_GENE|Entrez/) {
      $protein_node->addExisting($xref);
    }
    # For xrefs that should be attached to the gene node rather than the protein node, do that unless the xref already exists
    else {
      $gene_node->addExisting($xref) unless ($gene_node->nodeExists($xref));
    }
  }
  
  # Create a node for the long name to use if one was found
  if (length($name_content) > 0) {
    my $name_node = LRG::Node->new('long_name');
    $name_node->content($name_content);
    $protein_node->addExisting($name_node);
  }
  
  # Only add the protein node to the transcript if there is a coding exon within the LRG slice
  if ($outside) {
    return;
  }
    
  #ÊIf part of the protein falls outside the LRG slice, adjust the cds_start and cds_end to the closest coding exons
  if ($partial_5 || $partial_3) {
    my $feat_start = $transcript_node->data()->{'start'};
    my $feat_end = $transcript_node->data()->{'end'};
    my $exons = $transcript->get_all_Exons();
    my $exon_overlap = 0;
    if (scalar(@{$exons}) > 0) {
      $feat_start = $slice_end;
      $feat_end = 1;
      foreach my $exon (@{$exons}) {
        my $limits = get_feature_limits($exon,$lrg_slice);
        next if (!defined($limits) || $limits->{'upstream'} || $limits->{'downstream'});
	
	#ÊGet the position of the coding sequence start and end within the exon.
	# On the exon (as opposed to the transcript), coding_region_start ALWAYS referes to the start codon, so we need to pass coding_region_end if the exon is on the negative strand
	$lift = lift_coordinate(($exon->strand > 0 ? $exon->coding_region_start($transcript) : $exon->coding_region_end($transcript)),$exon->slice(),$lrg_slice);
	if ($lift->{'position'}) {
	  $feat_start = min($feat_start,$lift->{'position'});
	}
	$lift = lift_coordinate(($exon->strand > 0 ? $exon->coding_region_end($transcript) : $exon->coding_region_start($transcript)),$exon->slice(),$lrg_slice);
	if ($lift->{'position'}) {
	  $feat_end = max($feat_end,$lift->{'position'});
	}
      }
      if ($feat_start > $feat_end) {
        ($feat_start,$feat_end) = ($feat_end,$feat_start);
      }
    }
    $cds_start = $feat_start if ($partial_5);
    $cds_end = $feat_end if ($partial_3);
  }
  
  #ÊIn case the transcript is on the opposite strand relative to the LRG, flip the partial flags
  ($partial_5,$partial_3) = ($partial_3,$partial_5) if ($transcript->strand() < 0);
    
  # Add partial nodes if necessary
  $protein_node->addNode('partial')->content('5-prime') if ($partial_5);
  $protein_node->addNode('partial')->content('3-prime') if ($partial_3);
  
  # Add the cds data to the protein node
  $protein_node->addData(
    {
      'source' => 'Ensembl',
      'accession' => $protein->stable_id(),
      'cds_start' => $cds_start,
      'cds_end' => $cds_end
    }
  );
  
  $transcript_node->addExisting($protein_node);
}

sub attach_transcripts {
  my $gene = shift;
  my $lrg_slice = shift;
  my $gene_node = shift;
  
  #ÊGet the transcripts from the gene object
  my $transcripts = $gene->get_all_Transcripts();
  
  # Loop over the transcripts
  foreach my $transcript (@{$transcripts}) {

    # Check for partial overlaps
    my $limits = get_feature_limits($transcript,$lrg_slice);
    
    # If the transcript falls entirely outside of the LRG it should be skipped
    next if (!defined($limits) || $limits->{'upstream'} || $limits->{'downstream'});
    
    my $feat_start = $limits->{'start'};
    my $feat_end = $limits->{'end'};
    my $feat_strand = $limits->{'strand'};
    my $partial_5 = $limits->{'partial_5'};
    my $partial_3 = $limits->{'partial_3'}; 
		
    my $transcript_node = LRG::Node->new("transcript", undef, {'source' => 'Ensembl', 'start' => $feat_start, 'end' => $feat_end, 'transcript_id' => $transcript->stable_id()});
		
    # If transcript partially overlaps, it should be indicated
    $transcript_node->addNode('partial')->content('5-prime') if ($partial_5);
    $transcript_node->addNode('partial')->content('3-prime') if ($partial_3);
	  
    my $name_content = "";
	  
    # Extract xref information
    my $entries = $transcript->get_all_DBEntries();
    while (my $entry = shift(@{$entries})) {
  
      #ÊSkip the xref if the source is not one of the allowed ones	  	
      next unless $entry->dbname =~ /GI|RefSeq|MIM_GENE|Entrez|CCDS|RFAM|miRBase|pseudogene.org/;
      next if $entry->dbname =~ /RefSeq_peptide/;
	
      # Get the long name from RefSeq if it's available    
      if ($entry->dbname() eq 'RefSeq_dna' && defined($entry->description())) {
	$name_content = $entry->description();
      }
      
      # Create an xref node
      my $xref = xref($entry);
      # Watch out for the OMIM xrefs, they should go to the gene node      
      if ($entry->dbname !~ /MIM_GENE/) {
	$transcript_node->addExisting($xref);
      }
      else {
	$gene_node->addExisting($xref);
      }
    }
	  
    # Make sure to fill in something in the long_name element
    if (length($name_content) == 0 && defined($transcript->description())) {
      $name_content = $transcript->description();
    }
    if (length($name_content) == 0 && defined($transcript->external_name())) {
      $name_content = $transcript->external_name();
    }
    $name_content .= " (" . $transcript->biotype() . ")";
    
    # Create the long name node
    my $name_node = LRG::Node->new('long_name');
    $name_node->content($name_content);
    $transcript_node->addExisting($name_node);
	
    # Attach the transcript node to the gene node before attaching proteins  
    $gene_node->addExisting($transcript_node);
    attach_protein($transcript,$lrg_slice,$transcript_node);
  }
}

sub gene_2_feature {

  my $gene = shift;
  my $lrg_slice = shift;
  
  # Check if the boundaries of the gene extends beyond the mapped region, in which case the end points should be set to those of the mapped region and a flag indicating partial overlap should be set
  my $limits = get_feature_limits($gene,$lrg_slice);
  
  my $partial_5 = $limits->{'partial_5'};
  my $partial_3 = $limits->{'partial_3'};
  
  my $feat_slice = $gene->feature_Slice();
  
  my $feat_start = $limits->{'start'};
  my $feat_end = $limits->{'end'};
  my $feat_strand = $limits->{'strand'};
  
  # Create the gene node
  my $gene_node = LRG::Node->new("gene", undef, {'symbol' => $gene->external_name(), 'start' => $feat_start, 'end' => $feat_end, 'strand' => $gene->strand()});
  
  # If the gene partially overlaps, this should be indicated
  $gene_node->addNode('partial')->content('5-prime') if ($partial_5);
  $gene_node->addNode('partial')->content('3-prime') if ($partial_3);
  
  my $name_content = "";
  
  # Get xrefs for the gene
  my $entries = $gene->get_all_DBEntries();
  while(my $entry = shift(@{$entries})) {
    	
    # Skip if the xref source is not one of the allowed ones
    next unless $entry->dbname =~ /GI|RefSeq|HGNC$/;
    next if $entry->dbname =~ /RefSeq_peptide/;
    
    # Get synonyms from HGNC entry
    if($entry->dbname() eq 'HGNC') {
      
      # Add a lrg_gene_name node. This should be moved in a later processing step to the LRG branded section
      $gene_node->addNode('lrg_gene_name',{'source' => $entry->dbname()})->content($entry->display_id());
      
      foreach my $synonym (@{$entry->get_all_synonyms()}) {
        $gene_node->addNode('synonym')->content($synonym);
      }
      if (defined($entry->description())) {
        $name_content = $entry->description();
      }
    }
    
    # Add a xref node
    $gene_node->addExisting(xref($entry));
  }
  # Add a xref to Ensembl as well
  $gene_node->addEmptyNode('db_xref', {'source' => 'Ensembl', 'accession' => $gene->stable_id()});
  
  # Make sure to fill in something in the long_name element
  if (length($name_content) == 0 && defined($gene->description())) {
    $name_content = $gene->description();
  }
  if (length($name_content) == 0 && defined($gene->external_name())) {
    $name_content = $gene->external_name();
  }
  $name_content .= " (" . $gene->biotype() . ")";
  my $name_node = LRG::Node->new('long_name');
  $name_node->content($name_content);
  $gene_node->addExisting($name_node);
  
  # Attach the transcripts
  attach_transcripts($gene,$lrg_slice,$gene_node);
  
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
  my $last_strand = 1;
  my $chr_name = $mapping_node->data->{'chr_name'};
  my $chr_offset_start = $mapping_node->data->{'chr_start'};
  my $chr_offset_end = $mapping_node->data->{'chr_end'};
  my $lrg_offset_start = 1e11;
  my $lrg_offset_end = -1;
  
  foreach my $span (@{$mapping_spans}) {
    my $lrg_start = $span->data->{'lrg_start'};
    my $lrg_end = $span->data->{'lrg_end'};
    my $chr_start = $span->data->{'start'};
    my $chr_end = $span->data->{'end'};
    my $strand = $span->data->{'strand'};
    
    my $dna_pair = [
      'DNA',
      $lrg_start,
      $lrg_end,
      $chr_start,
      $chr_end,
      $strand
    ];
    
    # Store the extreme points of the LRG
    $lrg_offset_start = min($lrg_offset_start,$lrg_start);
    $lrg_offset_end = max($lrg_offset_end,$lrg_end);
    
    # If we have changed orientation since last span (or since start when sequence was forward), flip the sequence
    if ($strand != $last_strand) {
      reverse_comp(\$chr_seq);
      $last_strand = $strand;
    }
    
    #ÊAt each indel in the alignment, a new DNA pair is added to the array and a G(ap) pair as well
    my $diffs = $span->findNodeArray('diff');
    foreach my $diff (@{$diffs}) {
      my $pair;
      my $diff_lrg_start = $diff->data()->{'lrg_start'};
      my $diff_lrg_end = $diff->data()->{'lrg_end'};
      my $diff_chr_start = $diff->data()->{'start'};
      my $diff_chr_end = $diff->data()->{'end'};
      my $diff_lrg_seq = $diff->data()->{'lrg_sequence'};
      my $diff_chr_seq = $diff->data()->{'genomic_sequence'};
      if ($diff->data()->{'type'} eq 'mismatch') {
	$pair = [
	  'M',
	  $diff_lrg_start,
	  $diff_lrg_end,
	  $diff_chr_start,
	  $diff_chr_end,
	  $strand,
	  $diff_lrg_seq,
	  $diff_chr_seq
	];
	push(@pairs,$pair);
      }
      elsif ($diff->data()->{'type'} =~ m/[lrg|genomic]_ins/) {
	$pair = [
	  'G',
	  $diff_lrg_start,
	  $diff_lrg_end,
	  $diff_chr_start,
	  $diff_chr_end,
	  $strand,
	  $diff_lrg_seq,
	  $diff_chr_seq
	];
	push(@pairs,$pair);
	
	# End the current DNA pair, add it to the array and start a new one
	$dna_pair->[2] = ($diff_lrg_start - 1);
	$dna_pair->[4] = ($diff_chr_start - 1) if ($strand > 0);
	$dna_pair->[3] = ($diff_chr_end + 1) if ($strand < 0);
	push(@pairs,$dna_pair);
	$dna_pair = [
	  'DNA',
	  ($diff_lrg_end + 1),
	  $lrg_end,
	  ($strand > 0 ? ($diff_chr_end + 1) : $chr_start),
	  ($strand > 0 ? $chr_end : ($diff_chr_start - 1)),
	  $strand
	];
      }
    } 
    push(@pairs,$dna_pair);
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

#ÊConvert a pair array structure to an XML alignment structure
sub pairs_2_alignment {
  my $data = shift;
  
  #ÊA counter for the number of identical nucleotides in the alignment
  my $matches = 0;
  
  my $alignment_node = LRG::Node->new('alignment');
  # Should add the LRG and genomic offset to the coordinates?
  $alignment_node->addData({'query_start' => $data->{'lrg_start'}, 'query_end' => $data->{'lrg_end'}, 'target_start' => $data->{'chr_start'}, 'target_end' => $data->{'chr_end'}});
  my $span_node = LRG::Node->new('alignment_span',undef,$alignment_node->data());
  $span_node->addData({'strand' => $data->{'strand'}});
  $alignment_node->addExisting($span_node);
  
  my $pairs = $data->{'pairs'};
  foreach my $pair (@{$pairs}) {
    if ($pair->[0] eq 'DNA') {
      # Add the length of this alignment block to the count of identical nucleotides
      $matches += ($pair->[2] - $pair->[1] + 1);
    }
    else {
      my $diff_node = LRG::Node->newEmpty('diff');
      if ($pair->[0] eq 'M') {
	$diff_node->addData({'type' => 'mismatch'});
      }
      elsif ($pair->[0] eq 'G') {
	if ($pair->[2] >= $pair->[1]) {
	  $diff_node->addData({'type' => 'query_ins','query_sequence' => $pair->[6]});
	}
	else {
	  $diff_node->addData({'type' => 'target_ins','target_sequence' => $pair->[7]});
	}
      }
      $diff_node->addData({'query_start' => $pair->[1], 'query_end' => $pair->[2], 'target_start' => $pair->[3], 'target_end' => $pair->[4]});
      $span_node->addExisting($diff_node);
    }
  }
  
  # Calculate the similarity (#matches / total alignment length)
  my $similarity = sprintf("%.1f",(100*$matches)/$data->{'length'}) . "%";
  
  $alignment_node->addData({'query_name' => $data->{'lrg_id'}, 'target_name' => $data->{'chr_name'}, 'similarity' => $similarity, 'method' => 'exonerate-' . $EXONERATE_PARAMETERS{'--model'}});
  
  return $alignment_node;
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
  
  #ÊCreate a local deep copy of the array so we won't mess up the elements
  my @pairs_array;
  foreach my $pair (@{$pairs}) {
    my @copy;
    foreach my $p (@{$pair}) {
      push(@copy,$p);
    }
    push(@pairs_array,\@copy);
  }
  
  #ÊOrder the pairs array according to ascending LRG coordinates
  my @sorted_pairs = sort {$a->[1] <=> $b->[1]} @pairs_array;
  
  # Separate the DNA chunks from mismatches in the pairs array. We will create a mapping span for each DNA chunk and include all mismatch information in this span
  # Gaps can be discarded at this point since they are implicit between DNA chunks
  # FIXME: Gaps should be annotated as diffs, only "large" gaps (on the kb-scale) should be represented by different spans
  my @dna_chunks;
  my @mismatches;
  my @gaps;
  while (my $pair = shift(@sorted_pairs)) {
    push(@dna_chunks,$pair) if ($pair->[0] eq 'DNA');
    push(@mismatches,$pair) if ($pair->[0] eq 'M');
    push(@gaps,$pair) if ($pair->[0] eq 'G');
  }
  
  # Loop over the DNA chunk array and create mapping spans and diffs as necessary
  my $chr_map_start = 1e11;
  my $chr_map_end = -1;
  my $last_chr_start;
  my $last_chr_end;
  my $last_lrg_start;
  my $last_lrg_end;
  my $span_node;
  
  # Loop over the dna chunks and join those that are not sufficiently far apart
  for (my $i=1; $i<scalar(@dna_chunks); $i++) {
    # Only join if the chunks are in the same direction (i.e. no inversion)
    if ($dna_chunks[$i]->[5] == $dna_chunks[$i-1]->[5] && ($dna_chunks[$i]->[1] - $dna_chunks[$i-1]->[2]) < $MIN_GAP_BETWEEN_SPANS) {
      # If the orientation is the same, compare start and end coordinates
      if ($dna_chunks[$i]->[5] > 0 && ($dna_chunks[$i]->[3] - $dna_chunks[$i-1]->[4]) < $MIN_GAP_BETWEEN_SPANS) {
	$dna_chunks[$i-1]->[2] = $dna_chunks[$i]->[2];
	$dna_chunks[$i-1]->[4] = $dna_chunks[$i]->[4];
	splice(@dna_chunks,$i,1);
	$i--;
      }
      # If the orientation is reversed, compare the end and start coordinates
      if ($dna_chunks[$i]->[5] < 0 && ($dna_chunks[$i-1]->[3] - $dna_chunks[$i]->[4]) < $MIN_GAP_BETWEEN_SPANS) {
	$dna_chunks[$i-1]->[2] = $dna_chunks[$i]->[2];
	$dna_chunks[$i-1]->[3] = $dna_chunks[$i]->[3];
	splice(@dna_chunks,$i,1);
	$i--;
      }
      
    }
  }
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
    # Loop over the gaps and create diffs as necessary
    foreach my $gap (@gaps) {
      next unless ($gap->[1] >= $lrg_start && $gap->[2] <= $lrg_end);
      
      my $diff_node = LRG::Node->newEmpty('diff');
      $diff_node->addData(
	{
	  'type' => ((($gap->[2] - $gap->[1]) >= 0) ? 'lrg_ins' : 'genomic_ins'),
	  'lrg_start' => $gap->[1],
	  'lrg_end' => $gap->[2],
	  'start' => $gap->[3],
	  'end' => $gap->[4]
	}
      );
      $diff_node->addData({'lrg_sequence' => $gap->[6]}) if (defined($gap->[6]) && length($gap->[6]));
      $diff_node->addData({'genomic_sequence' => $gap->[7]}) if (defined($gap->[7]) && length($gap->[7]));
      
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
  return get_overlapping_features(@_,'gene');
}

sub get_overlapping_features {
  my $lrg_slice = shift;
  my $feature_type = shift;
  
  $feature_type ||= 'gene';
  
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
  my $features;
  if ($feature_type eq 'gene') {
    $features = $ref_slice->get_all_Genes();
  }
  elsif ($feature_type eq 'regulatory') {
  
    # Get a FeatureSet adaptor to get the feature sets
    my $feature_set_adaptor = $dbFuncGen->get_FeatureSetAdaptor();
    #ÊGet all Feature Sets (might want to limit these?)
    my $fsets = $feature_set_adaptor->fetch_all();
    my @reg_features;
    # Loop over all Feature Sets
    while (my $fset = shift(@{$fsets})) {
      # Get all Regulatory Features in the set on the LRG slice
      push(@reg_features,@{$fset->get_Features_by_Slice($ref_slice)});
    }
    $features = \@reg_features;
  }
	
  
  return $features;
}

sub regulatory_feature_2_feature {
  my $reg_feature = shift;
  my $lrg_slice = shift;
  
  # Create a new node for this regulatory feature
  my $reg_feature_node = LRG::Node->new('regulatory_element');
  
  # FIXME: Need to translate coordinates to LRG, for now, use chromosome coordinates
  $reg_feature_node->addData({'name' => $reg_feature->feature_type()->name(),'start' => $reg_feature->start(), 'end' => $reg_feature->end(), 'strand' => $reg_feature->strand()});
  
  # FIXME: Check if there is partial overlap
  
  #ÊAdd an element with the Regulatory Feature class
  if (length($reg_feature->feature_type()->class()) > 0) {
    my $class_node = LRG::Node->new('class');
    $class_node->content($reg_feature->feature_type()->class());
    $reg_feature_node->addExisting($class_node);
  }
  #ÊAdd an element with the Regulatory Feature display label
  if (length($reg_feature->display_label()) > 0) {
    my $label_node = LRG::Node->new('label');
    $label_node->content($reg_feature->display_label());
    $reg_feature_node->addExisting($label_node);
  }
  #ÊAdd an element with the Regulatory Feature description
  if (length($reg_feature->feature_type()->description()) > 0) {
    my $desc_node = LRG::Node->new('description');
    $desc_node->content($reg_feature->feature_type()->description());
    $reg_feature_node->addExisting($desc_node);
  }
  
  # Add cell type data if available
  if (defined($reg_feature->cell_type())) {
    my $ct = $reg_feature->cell_type();
    my $ct_node = LRG::Node->new('cell_type');
    $ct_node->addData({'name' => $ct->name()});
    
    # Add gender information if available
    if (defined($ct->gender())) {
      my $node = LRG::Node->new('gender');
      $node->content($ct->gender());
      $ct_node->addExisting($node);
    }
    # Add label if available
    if (defined($ct->display_label())) {
      my $node = LRG::Node->new('label');
      $node->content($ct->display_label());
      $ct_node->addExisting($node);
    }
    # Add description if available
    if (defined($ct->description())) {
      my $node = LRG::Node->new('description');
      $node->content($ct->description());
      $ct_node->addExisting($node);
    }
    
    $reg_feature_node->addExisting($ct_node);
  }
  
  # If this is a RegulatoryFeature object which have an Ensembl stable id, store this as a xref
  if ($reg_feature->isa('Bio::EnsEMBL::Funcgen::RegulatoryFeature')) {
    my $xref_node = LRG::Node->new('db_xref');
    $xref_node->addData({'source' => 'Ensembl', 'accession' => $reg_feature->stable_id()});
    $reg_feature_node->addExisting($xref_node);
  }
  
  return $reg_feature_node;
}

1;
