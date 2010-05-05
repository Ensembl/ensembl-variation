use strict;
use warnings;

package LRGMapping;

use Bio::EnsEMBL::MappedSliceContainer;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use LRG;
use List::Util qw (min max);


# global variables
our $input_dir = '/tmp';
our $output_dir = '/tmp';
our $target_dir;
our %rec_seq;

our $lrg_name;
our $dbCore;
our $registry_file;
our $current_assembly;

our $mapping_num = 1;


sub mapping {
  my $sequence = shift;
  my $seq_length = length($sequence);
	
  #print "LENGTH $seq_length\n";
	
  # if we need to split it
  if($seq_length > 25000) {
    my $split_size = 10000;
		
    my $current_start = 0;
		
    my @maps = ();
	
    while($current_start < $seq_length) {
			
      # split_size is conservative; if this is the penultimate split, just do a bigger
      # one to make sure we don't try and map too small a sequence
      if($seq_length - $current_start < (2 * $split_size)) {
		$split_size = 2 * $split_size;
      }
		
      #print "\nMapping from $current_start to ", $current_start + $split_size, "\n";
      push @maps, mapping(substr($sequence, $current_start, $split_size));
			
      $current_start += $split_size;
    }
	
	#print "\n>>> Done running ssaha\n";

	my $prev_end = 0;
	my @dna_pairs = ();
    my @pairs = ();
	my @mis_pairs = ();
    
	my $full_match = 1;

	# get all pairs from all maps into arrays
	foreach my $map(@maps) {
		
		$full_match = 0 if ($map->identical_matches == 0);
		
		# sort by query start
		foreach my $pair(sort {$a->[2] <=> $b->[2]} @{$map->type}) {
			#print "full_match is $full_match and ",$map->identical_matches,'-',$pair->[0],'-',$pair->[1],'-',$pair->[2],'-',$pair->[3],'-',$pair->[4],'-',$pair->[5],'-',$pair->[6],'-',$pair->[7],"\n";
			
			#print "@$pair\n";
			
			# add the previous DNA end to each query coordinate
			# since we split up the query sequence
			$pair->[1] += $prev_end;
			$pair->[2] += $prev_end;
			
			# put the DNA pairs in a separate array
			if($pair->[0] eq 'DNA') {
				push @dna_pairs, $pair;
			}
			
			# put mismatch pairs in the @mis_pairs array
			elsif(defined($pair->[6]) && defined($pair->[7]) && $pair->[6].$pair->[7] =~ /.+/) {
				push @mis_pairs, $pair;
			}
			
			# put the match pairs in the @pairs array
			else {
				push @pairs, $pair;
			}
		}
		
		# get prev_end from the last added DNA pair
		$prev_end = $dna_pairs[-1]->[2];
	}
	
	# now join the pairs using join_pairs subroutine
	my @joined_pairs = (@{join_pairs(\@pairs)}, @{join_pairs(\@mis_pairs)}, @{join_pairs(\@dna_pairs)});
	
	# join the maps
	my $main_map = shift @maps;
	
	my $prev_q_end = $main_map->end;
	
	foreach my $map(@maps) {
		if($map->hend < $main_map->hend && $map->hstart < $main_map->hstart) {
		  $main_map->hstart($map->hstart);
		}
		
		else {
		  $main_map->hend($map->hend);
		}
	  
		$main_map->end( $map->end + $prev_q_end);
		  
		$prev_q_end = $main_map->end;
	}
	
    $main_map->identical_matches($full_match);
    my $seqname = $main_map->seqname;
    chop($seqname);
    $main_map->seqname($seqname);
    
	#print "main_map is ",$main_map->seqname,'-',$main_map->hseqname,'-',$main_map->start,'-',$main_map->end,'-',$main_map->hstart,'-',$main_map->hend,'-',$main_map->identical_matches,"\n";
    
	# add the pairs onto the joined map
    $main_map->type(\@joined_pairs);

    #foreach my $pair (@{$main_map->type}) {
    #  print "checking DNA bit : ",$pair->[0],'-',$pair->[1],'-',$pair->[2],'-',$pair->[3],'-',$pair->[4],'-',$pair->[5],'-',$pair->[6],'-',$pair->[7],"\n";
    #}
	
    return $main_map;
  }
	
  else {
	
    my $name = "temp$$".$mapping_num++;
	#my $name = "temp19098".$mapping_num++;
    #my $name = "temp214871".$mapping_num++;
    #my $name = "temp_will";
    my $input_file_name = $name.'.fa';
		
    open OUT, ">$input_dir/$name\.fa" or die $!;
    #open OUT, ">$name\.fa" or die $!;
    print OUT '>', $name, "\n", $sequence;
    close OUT;
		
    my $queue = "-q normal -R'select[mem>4000] rusage[mem=4000]' -M4000000";
		
    my $output_file_name = "ssaha2_output\_$input_file_name";
    my $input_file = "$input_dir/$input_file_name";
    my $output_file ="$output_dir/$output_file_name";
		
    my ($subject, %rec_find, %input_length, %done);
		
    $subject = "$target_dir/ref";
    $rec_seq{$name} = $sequence;
    #print "q_seq length is ",length($sequence),"\n";			
    my $seqobj = Bio::PrimarySeq->new(-id => $name, -seq => $rec_seq{$name});
    $rec_seq{$name} = $seqobj;
		
    bsub_ssaha_job($queue,$input_file,$output_file,$subject);

    my $call = "bsub -o $input_dir/$$.waiting.out -P ensembl-variation -K -w 'done($input_file\_ssaha_job)' -J waiting_process sleep 1"; #waits until all ssaha jobs have finished to continue
    system($call);
	
	#unlink $input_file;
	
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
	
    return $mapping->[0];
  }
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

# subroutine to join mapping elements
sub join_pairs {
	my $pairs = shift;
	my $prev_pair;
	my @joined_pairs;
	
	foreach my $pair(@$pairs) {
		if(defined $prev_pair) {
			
			# if the start of this pair is 1 more than the
			# end of the previous one we want to join them
			if($pair->[1] - $prev_pair->[2] == 1) {
				
				# copy this pair's query end to prev_pair
				$prev_pair->[2] = $pair->[2];
				
				# add any bases that were in a mismatch pair
				$prev_pair->[6] .= $pair->[6] unless (!defined($pair->[6]));
				$prev_pair->[7] .= $pair->[7] unless (!defined($pair->[7]));
				
				# do the same with the target end
				# but we have to work out which it is since target
				# coords are always reported start < end
				if($pair->[3] < $prev_pair->[3] && $pair->[4] < $prev_pair->[4]) {
					$prev_pair->[3] = $pair->[3];
				}
				
				else {
					$prev_pair->[4] = $pair->[4];
				}
			}
			
			# if the two pairs shouldn't be joined
			else {
				push @joined_pairs, $prev_pair;
				
				# clear prev_pair so that code below works
				$prev_pair = undef;
			}
		}
		
		# copy this pair to prev_pair if we're not joining to previous
		$prev_pair = $pair unless defined $prev_pair;
	}
	
	# clean up remaining prev_pair by adding to joined_pairs
	# also allows for situation where there's only 1 pair
	push @joined_pairs, $prev_pair if defined $prev_pair;
	
	return \@joined_pairs;
}

# subroutine to parse the output from ssaha2
sub parse_ssaha2_out {

#   open OUT, ">$output_dir/mapping_file_$input_file_name" or die "can't open output_file : $!\n";
  my $output_file = shift;

  open SSAHA, "$output_file" or die "can't open $output_file : $!\n";
  my @rec_find;
	
  while (<SSAHA>) {
    #print ;
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
}

# processes and creates mapping pairs
sub make_feature_pair {

  my (@pairs,@fps);
  my ($h) = @_;

  my $q_id = $h->{'q_id'};
  my $q_start = $h->{'q_start'};
  my $q_end = $h->{'q_end'};
  my $q_strand = $h->{'q_strand'};
  my $t_id = $h->{'t_id'};
  my $t_start = $h->{'t_start'};
  my $t_end = $h->{'t_end'};
  my $t_strand = $h->{'t_strand'};
  my $score = $h->{'score'};
  my $match = $h->{'match'};

  #print "$q_id,$q_start,$q_end,$q_strand,$t_id,$t_start,$t_end,$t_strand,$score and match $match\n";
  my ($f_q_start,$f_q_end) ;#feature start-end used to make feature_pairs
  if ($q_strand ==1) {
    ($f_q_start,$f_q_end) = ($q_start,$q_end);
  }
  else {
    ($f_q_start,$f_q_end) = ($q_end,$q_start);
  }

  my ($seq_region_name) = split /\-/, $t_id;
  
  my $sa = $dbCore->get_SliceAdaptor();
  my $slice = $sa->fetch_by_region('chromosome',$seq_region_name);

  my $q_seqobj = $rec_seq{$q_id};

  #warning that query_seq length != matched query_length
  if ($q_seqobj and length($q_seqobj->seq) != abs($q_end - $q_start) + 1){
    #print "q_start is $q_start and q_end is $q_end and length is ",length($q_seqobj->seq),"\n";
    die("query sequence length not equal matched query sequence length");
  }

  my @match_components = split /\s+/, $match;

  my ($full_match,$new_q_start, $new_q_end, $new_t_start, $new_t_end);

  $full_match=1; #initial full_match

  while (@match_components) {

    my $type = shift @match_components;
    my $query_match_length = shift @match_components;
    my $target_match_length = shift @match_components;
    my ($tmp_q_start,$tmp_q_end);
      
    if ($type eq 'M') {#go through each base to check SNP,target strand always 1
      my ($q_seq);
      if ($q_strand == 1) {
	$new_q_start = $q_start + $query_match_length - 1 ;
	$q_seq = substr($q_seqobj->seq,$q_start-1,$query_match_length);
      }
      else {
	$new_q_start = $q_start - $query_match_length + 1;
	$q_seq = substr($q_seqobj->seq,$new_q_start-1,$query_match_length);
	reverse_comp(\$q_seq);
      }
      $new_t_start = $t_start + $target_match_length -1;
      my $tmp_qs = ($q_strand ==1) ? $q_start : $new_q_start;
      my $tmp_qe = ($q_strand ==1) ? $new_q_start : $q_start;
      push @pairs, ['DNA',$tmp_qs,$tmp_qe,$t_start,$new_t_start,$q_strand];

      #print "$q_start,$new_q_start,$t_start,$new_t_start\n";
      #my $t_seq = substr($slice->seq,$t_start-1,$new_t_start-$t_start+1);
	  
	  my $t_end = $t_start + $target_match_length -1;
	  my $t_seq = $slice->sub_Slice($t_start, $t_end)->seq;
	  
      #print "$q_start,$new_q_start,$t_start,$new_t_start and length of q_seq ",length($q_seq),'-',length($t_seq),"\n";
      #print "q_seq_is $q_seq\n";
      #print "t_seq is $t_seq\n";
      my $q_count = 1;
      my $t_count = 1;
	  
      my %q_seqs = map {$q_count++,$_} split '', $q_seq;
      my %t_seqs = map {$t_count++,$_} split '', $t_seq;
      my ($sub_q_end,$new_count);
      my $sub_q_start = $q_start;#initial sub_q_start
      my $sub_t_start = $t_start;
      foreach my $count (sort {$a<=>$b} keys %q_seqs) {
	if (defined $q_seqs{$count+1} and defined $t_seqs{$count+1} and $q_seqs{$count+1} !~ /$t_seqs{$count+1}/i and $q_seqs{$count} =~ /$t_seqs{$count}/i and $count < length($q_seq)) {
	  push @pairs, [$type,$q_start,$q_start+$count-1,$t_start,$t_start+$count-1,$q_strand] if $q_strand==1;
	  push @pairs, [$type,$q_start,$q_start-$count+1,$t_start,$t_start+$count-1,$q_strand] if $q_strand==-1;
	}
	if ($q_seqs{$count} !~ /$t_seqs{$count}/i) {
	  #if there is a mismatch, we need to record the base based on query sequence
          if ($q_strand==-1) {
            reverse_comp(\$q_seqs{$count});
            reverse_comp(\$t_seqs{$count});
          }

	  #print "count is $count\n";
	  $full_match =0;
	  my $sub_t_end = $t_start + ($count) - 1;
	  $sub_t_start = $sub_t_end ;

	  if ($q_strand ==1) {
	    $sub_q_end = $q_start + ($count -1) - 1;
	  }
	  else {
	    $sub_q_end = $q_start - ($count -1) + 1;
	  }
	  if ($q_strand==1) {
	    $tmp_q_start = $sub_q_end+1;#change start to end
	    $tmp_q_end = $sub_q_end+1; #change start to end
	    $sub_q_start = $sub_q_end+2;#+1 cover SNP base, +1 cover next base
	  }
	  else {
	    $tmp_q_start = $sub_q_end-1;
	    $tmp_q_end = $sub_q_end-1;
	    $sub_q_start = $sub_q_end-2;
	  }
	
	  ($tmp_q_start,$tmp_q_end) = ($tmp_q_end,$tmp_q_start) if($tmp_q_end<$tmp_q_start);
	  push @pairs, [$type,$tmp_q_start,$tmp_q_end,$sub_t_start,$sub_t_end,$q_strand,uc($q_seqs{$count}),uc($t_seqs{$count})];

	  #print "in Mismatch, q_start is $tmp_q_start and q_end is $tmp_q_end and t_start is $sub_t_start and t_end is $sub_t_end and q_base is,uc($q_seqs{$count}) and target_base is ,uc($t_seqs{$count}) \n";
	  $sub_t_start = $sub_t_end+1;
	}
      }

      $q_start = $sub_q_start;
      $t_start = $sub_t_start;
	
      #print "q_start is $q_start and t_start is $t_start\n";
      #print "in M, q_start is $q_start and q_end is $new_q_start and t_start is $t_start and t_end is $new_t_start\n";

      if  ($q_strand==1) {
	$tmp_q_start = $q_start;
	$tmp_q_end = $new_q_start;
      }
      else {
	$tmp_q_start = $new_q_start;
	$tmp_q_end = $q_start;
      }
      ($tmp_q_start,$tmp_q_end) = ($tmp_q_end,$tmp_q_start) if($tmp_q_end<$tmp_q_start);
      push @pairs, [$type,$tmp_q_start,$tmp_q_end,$t_start,$new_t_start,$q_strand];
    }
    elsif ($type eq 'G') {
    my ($q_seq);
      if ($q_strand ==1) {
	$new_q_start = $q_start + $query_match_length + 1;
        $q_seq = substr($q_seqobj->seq,$q_start,$query_match_length);
        $q_seq = ($q_seq) ? $q_seq : '-';
      }
      else {
	$new_q_start = $q_start - $query_match_length - 1;
        $q_seq = substr($q_seqobj->seq,$new_q_start,$query_match_length);
        $q_seq = ($q_seq) ? $q_seq : '-';
      }
      $new_t_start = $t_start + $target_match_length +1;
      my $t_seq = substr($slice->seq, $t_start, $target_match_length);
      if ($q_strand ==-1) {
        reverse_comp(\$t_seq) if $t_seq;
      }
      $t_seq = ($t_seq) ? $t_seq : '-';

      #the following part just to make tmp start-end to put into pairs
      if ($q_strand ==-1) {
	$tmp_q_start = $new_q_start+1;
	$tmp_q_end = $q_start-1;
      }
      else {
	$tmp_q_start = $q_start+1;
	$tmp_q_end = $new_q_start-1;
      }
      ($tmp_q_start,$tmp_q_end) = ($tmp_q_end,$tmp_q_start) if($tmp_q_end<$tmp_q_start);
      push @pairs, [$type,$tmp_q_start,$tmp_q_end,$t_start+1,$new_t_start-1,$q_strand,uc($q_seq),uc($t_seq)];

      #print "in G, q_start is $q_start and q_end is $new_q_start and t_start is $t_start and t_end is $new_t_start\n";
    }

    $q_start = $new_q_start;
    $t_start = $new_t_start;
  }

  foreach  my $pair (sort {$a->[3]<=>$b->[3]} @pairs) {
    print "pairs in ssaha_mapping are @$pair\n";
  }
  my $fp = Bio::EnsEMBL::FeaturePair->new(-start    => $f_q_start,
					  -end      => $f_q_end,
					  -strand   => $q_strand,
					  -display_id => $h->{'q_id'},
					  -hstart   => $h->{'t_start'},
					  -hend     => $h->{'t_end'},
					  -hstrand  => $h->{'t_strand'},
					  -hseqname => $h->{'t_id'},
					  -slice    => $slice,
					 );
  $fp->seqname($h->{'q_id'});
  $fp->type(\@pairs);
  $fp->identical_matches($full_match);

  push @fps,$fp;

  return \@fps;
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
  my $fp = shift;
  my $q_seq = shift;
  
  die("LRG name not defined\n") unless defined $lrg_name;

# For now, just assume one mapping span so get the first element from the feature pair array argument
  $fp = $fp->[0];

  my ($q_start,$q_end,$t_start,$t_end,$q_strand);
  my $slice = $fp->slice;
  $q_start = $fp->start;
  $q_end = $fp->end;
  $t_start = $fp->hstart;
  $t_end = $fp->hend;
  $q_strand = $fp->strand;
  my $pairs = $fp->type;
  
 # my $sub_slice = $slice->sub_Slice($t_start,$t_end,$q_strand);
  my $full_match = $fp->identical_matches;
  
  my @nodes;

  # For consistency in the annotations, do transfer between LRG and chromosome coord systems
  #Êeven if there is a perfect match
  my $lrg_coord_system = $lrg_name;
  my $additions = add_mapping($lrg_name,$lrg_coord_system,length($q_seq),$current_assembly,$slice->seq_region_name(),$pairs,$dbCore);
  
  Bio::EnsEMBL::Registry->clear();
  Bio::EnsEMBL::Registry->load_all( $registry_file );
  $dbCore = Bio::EnsEMBL::Registry->get_DBAdaptor('human','core');
  
  my $sa = $dbCore->get_SliceAdaptor();
  my $lrg_slice = $sa->fetch_by_region($lrg_coord_system,$lrg_name);

  my $min = 9999999999;
  my $max = -9999999999;
  my $chrom;
  my $strand;

  foreach my $segment (@{$lrg_slice->project('chromosome',$current_assembly)}) {
    my $from_start = $segment->from_start();
    my $from_end    = $segment->from_end();
    my $to_name    = $segment->to_Slice->seq_region_name();
    $chrom = $to_name;

    my $to_start    = $segment->to_Slice->start();
    my $to_end    = $segment->to_Slice->end();
    if($to_start > $max){
	$max = $to_start;
    }
    if($to_start < $min){
	$min = $to_start;
    }
    if($to_end > $max){
	$max = $to_end;
    }
    if($to_end <  $min){
	$min = $to_end;
    }
    my $ori        = $segment->to_Slice->strand();
    $strand = $ori;   
  
    #print "$from_start-$from_end  => $to_name $to_start-$to_end ($ori) \n";
  }

   # foreach my $gene (@{$lrg_slice->get_all_Genes()}){
      #print "gene_name is ",$gene->stable_id,"\n";
    #}
    #print "q_seq from lrg_slice is ",$lrg_slice->start,'-'.$lrg_slice->end,"\n";
    #print "length of q_seq is ",length($q_seq), " and length of hseq is ", length($lrg_slice->seq),"\n";
    
    #if ($lrg_slice->seq eq $q_seq) {
    #  print "hseq is same as q_seq\n";
    #}
    #else {
    #  print "hseq is different from q_seq\n";
    #  #print "q_seq is ",$q_seq,"\n";
    #  #print "hseq is ",$lrg_slice->seq,"\n";
    #}

   my $ref_slice = $sa->fetch_by_region("chromosome",$chrom, $min, $max, $strand);
   my @genes = @{$ref_slice->get_all_Genes()};
	
   foreach my $gene (@genes) {
#     next unless $gene->biotype() eq 'protein_coding';
		
     #print "before transfprm g start-end ",$gene->start,'-',$gene->end,"\n";
     #my $new_gene = $gene->transform('LRG');
#     my $new_gene = $gene->transfer($lrg_slice);
     #print "after transform g start-end ",$new_gene->start,'-',$new_gene->end,"\n" if $new_gene;
	  
    # Pass the gene as it is on the reference slice instead and do a transformation for each transcript separately in the subroutine	  
     my @new_nodes = @{get_gene_annotation($gene,$lrg_slice)};
	  
     push @nodes, @new_nodes;
   }
 # }
  clean_mapping($additions,$dbCore);
  
  return \@nodes;
}

sub get_gene_annotation {

	my $gene = shift;
	my $lrg_slice = shift;
	my $slice_end = $lrg_slice->end;
	my ($current, $entry);
	my @nodes;
	
	my $lrg_gene = $gene->transfer($lrg_slice);
	
# Check if the boundaries of the gene extends beyond the mapped region, in which case the end points should be set to those of the mapped region and a flag indicating partial overlap should be set

	my $feat_start = max($gene->start,1);
	my $feat_end = min($gene->end,$slice_end);
	my $feat_strand = $lrg_gene->strand;
	
	my $partial_5 = ($feat_start != $gene->start);
	my $partial_3 = ($feat_end != $gene->end);
	if ($feat_strand < 0) {
	  ($partial_5,$partial_3) = ($partial_3,$partial_5);
	}
	
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
	
	foreach my $trans(@{$gene->get_all_Transcripts}) {
		
		#print "Gene/Trans ", $gene->stable_id, " ", $trans->stable_id, "\n";
		
    ### There is still a problem with this in case sub-features (exons) are crossing the coord system border ###
    ### Use Ian's suggestion instead? LRG_12 is an example of this (ENSG00000141698). 
	  my $lrg_trans = $trans->transfer($lrg_slice);
	  next if (!defined($lrg_trans));
	  
# Check for partial overlaps
		$feat_start = max($trans->start,1);
		$feat_end = min($trans->end,$slice_end);
		$feat_strand = $trans->strand;
		$partial_5 = ($feat_start != $trans->start);
		$partial_3 = ($feat_end != $trans->end);
		if ($feat_strand < 0) {
		  ($partial_5,$partial_3) = ($partial_3,$partial_5);
		}
	      	
# If the transcript falls entirely outside of the LRG it should be skipped
		next if ($feat_end < 1 || $feat_start > $slice_end);
		
#ÊIf the transcript only partially overlaps the LRG, we need to set the end points of the transcript to the end point of the last exon that falls within the LRG record
		my $exon_overlap = 1;
		
		if ($partial_5 || $partial_3) {
		  my $exons = $trans->get_all_Exons();
		  $exon_overlap = 0;
		  if (scalar(@{$exons}) > 0) {
		    $feat_start = $slice_end;
		    $feat_end = 1;
		    foreach my $exon (@{$exons}) {
		      my $lrg_exon = $exon->transfer($lrg_slice);
		      next if (!defined($lrg_exon));
# Set a flag to indicate that there is in fact an exon in this transcript within the scope of the LRG
		      $exon_overlap = 1;
		      if ($exon->end > 1) {
			$feat_start = min($feat_start,max(1,$exon->start));
		      }
		      if ($exon->start < $slice_end) {
			$feat_end = max($feat_end,min($slice_end,$exon->end));
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
			
			next unless $entry->dbname =~ /GI|RefSeq|MIM_GENE|Entrez|CCDS|RFAM|miRBase/;
			next if $entry->dbname =~ /RefSeq_peptide/;
			
			if($entry->dbname eq 'RefSeq_dna') {
			  my $note_node = LRG::Node->new('long_name');
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
		  		  
# Check for partial overlaps
		  $feat_start = max($trans->coding_region_start,$feat_start);
		  $feat_end = min($trans->coding_region_end,$feat_end);
		  $partial_5 = ($feat_start != $trans->coding_region_start);
		  $partial_3 = ($feat_end != $trans->coding_region_end);
		  if ($feat_strand < 0) {
		    ($partial_5,$partial_3) = ($partial_3,$partial_5);
		  }
# Some gene-level data (e.g. NCBI GeneID) is buried in the translation object. So before determining if we skip the transcript, parse the xrefs
		  my $prot_node = LRG::Node->new('protein_product');
		  $prot_node->data(
		    {
		      'source' => 'Ensembl',
		      'accession' => $protein->stable_id,
		      'cds_start' => $feat_start,
		      'cds_end' => $feat_end
		    }
		  );
		  
		  $entries = $protein->get_all_DBEntries();
			  
		  while($entry = shift @$entries) {
				  
		    next unless $entry->dbname =~ /RefSeq|Uniprot|CCDS|MIM_GENE|Entrez/;
				  
		    if($entry->dbname eq 'RefSeq_peptide') {
		      my $note_node = LRG::Node->new('long_name');
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
		  unless ($feat_end < 1 || $feat_start > $slice_end) {
	      
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
	
	unshift @nodes, $gene_node;
	
	return \@nodes;
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

# Add all the necessary data to the db in order to project and transfer between the LRG and chromosome
# Returns a hashref holding the tables affected and conditions for retrieving the changes
sub add_mapping {
  my $lrg_name = shift;
  my $lrg_coord_system = shift;
  my $q_seq_length = shift;
  my $assembly = shift;
  my $target_chr = shift;
  my $pairs = shift;
  my $dbCore = shift;

  my %additions;
  
  # Get the coord system id for this LRG or add a new entry if not present
  my $stmt = qq{
    SELECT
      coord_system_id
    FROM
      coord_system
    WHERE
      name = '$lrg_coord_system'
    LIMIT 1
  };
  my $cs_id = $dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
  if (!defined($cs_id)) {
    $stmt = qq{
      SELECT
	MAX(rank)
      FROM
	coord_system
    };
    my $rank = ($dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0]+1);
    $stmt = qq{
      INSERT INTO
	coord_system (
	  species_id,
	  name,
	  rank,
	  attrib
	)
      VALUES (
	1,
	'$lrg_coord_system',
	$rank,
	'default_version'
      )
    };
    $dbCore->dbc->do($stmt);
    $cs_id = $dbCore->dbc->db_handle->{'mysql_insertid'};
    $additions{'coord_system'} = 'coord_system_id = ' . $cs_id;
  }

  # In order to project the slice, we need to add an entry to the meta table (if not present)
  my $meta_value = $lrg_coord_system . '#chromosome:' . $assembly;
  $stmt = qq{
    SELECT
      meta_id
    FROM
      meta
    WHERE
      meta_key = 'assembly.mapping' AND
      meta_value LIKE '$meta_value%'
  };
  my $ref = $dbCore->dbc->db_handle->selectall_arrayref($stmt);
  my @meta_id;
  foreach my $id (@{$ref}) {
    push(@meta_id,$id->[0]);
  }
  if (scalar(@meta_id) == 0) {
    foreach my $val (($meta_value, $meta_value . '#contig', $meta_value . '#supercontig')) {
      $stmt = qq{
        INSERT INTO
          meta (
            species_id,
            meta_key,
            meta_value
          )
        VALUES (
          1,
          'assembly.mapping',
          '$val'
        )
      };
      $dbCore->dbc->do($stmt);
      push(@meta_id,$dbCore->dbc->db_handle->{'mysql_insertid'});
    }
    $additions{'meta'} = 'meta_id = ' . join(' OR meta_id = ',@meta_id);
  }

  $stmt = qq{
    SELECT
      seq_region_id
    FROM
      seq_region
    WHERE
      name = '$lrg_name' AND
      coord_system_id = '$cs_id'
  };
  my $q_seq_region_id = $dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
  if (!defined($q_seq_region_id)) {
    $stmt = qq{
      INSERT INTO
	seq_region (
	  name,
	  coord_system_id,
	  length
	)
      VALUES (
	'$lrg_name',
	$cs_id,
	$q_seq_length
      )
    };
    $dbCore->dbc->do($stmt);
    $q_seq_region_id = $dbCore->dbc->db_handle->{'mysql_insertid'};
    $additions{'seq_region'} = 'seq_region_id = ' . $q_seq_region_id;
  }
  
  # Get the seq_region_id for the target chromsome
  $stmt = qq{
    SELECT
      sr.seq_region_id
    FROM
      seq_region sr,
      coord_system cs
    WHERE
      sr.name = '$target_chr' AND
      sr.coord_system_id = cs.coord_system_id AND
      cs.name = 'chromosome' AND
      cs.version = '$assembly'
    LIMIT 1
  };
  my $t_seq_region_id = $dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
  
  # Add mappings between the LRG and chromosome coordinate systems
  my @conditions;
  foreach  my $pair (sort {$a->[3]<=>$b->[3]} @{$pairs}) {
    if ($pair->[0] eq 'DNA') {
      if ($pair->[2]-$pair->[1] == $pair->[4]-$pair->[3]) {
	my $condition = qq{
	  cmp_seq_region_id = $t_seq_region_id AND 
	  asm_seq_region_id = $q_seq_region_id AND 
	  cmp_start = $pair->[3] AND 
	  cmp_end = $pair->[4] AND 
	  asm_start = $pair->[1] AND 
	  asm_end = $pair->[2] AND 
	  ori = $pair->[5]
	};
	$stmt = qq{
	  SELECT EXISTS (
	    SELECT
	      *
	    FROM
	      assembly
	    WHERE
	      $condition
	  )
	};
	if (!$dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0]) {
	  $stmt = qq{
	    INSERT IGNORE INTO
	      assembly (
	        cmp_seq_region_id,
	        asm_seq_region_id,
	        cmp_start,
	        cmp_end,
	        asm_start,
	        asm_end,
	        ori
	      )
	    VALUES (
	      $t_seq_region_id,
	      $q_seq_region_id,
	      $pair->[3],
	      $pair->[4],
	      $pair->[1],
	      $pair->[2],
	      $pair->[5]
	    )  
	  };
	  $dbCore->dbc->do($stmt);
	  push(@conditions,$condition);
	}
      }
      else {
	die("distance between query and target is not the same, there is a indel");
      }
    }
  }
  $additions{'assembly'} = '(' . join(') OR (',@conditions) . ')';
  return \%additions;
}

sub clean_mapping {
  my $additions = shift;
  my $dbCore = shift;
  
  foreach my $table (keys %{$additions}) {
    next if ($additions->{$table} eq '()'); 
    my $stmt = qq{
      DELETE FROM
	$table
      WHERE
	$additions->{$table}
    };
    $dbCore->dbc->do($stmt);
  }
}
