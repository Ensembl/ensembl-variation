#!/software/bin/perl

use strict;
use lib '/nfs/acari/dr2/projects/src/ensembl/ensembl/modules';
use Getopt::Long;
use LRG;
use Bio::EnsEMBL::Registry;
use Bio::Seq;
use Bio::EnsEMBL::MappedSliceContainer;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::FeaturePair;
use ImportUtils qw(dumpSQL debug create_and_load load );
use FindBin qw( $Bin );
use Bio::EnsEMBL::Variation::Utils::Sequence qw(unambiguity_code);
use Bio::EnsEMBL::Utils::Sequence qw(expand reverse_comp);
use Data::Dumper;

# set default options
my $out_file_stem;
my $in_file_name;
my $template_file_name = 'template.xml';
my $registry_file = "ensembl.registry";
my $input_dir;
my $output_dir;
my $target_dir;

# get options from command line
GetOptions(
	   'in_file_name=s' => \$in_file_name,#file in format : LRG_5   LEPRE1  NG_008123   NM_022356   NP_071751
	   'output_file_stem=s' => \$out_file_stem,
	   'template_file_name=s' => \$template_file_name,
	   'registry_file=s' => \$registry_file,
	   'input_dir=s' => \$input_dir,
	   'output_dir=s' => \$output_dir,
	   'target_dir=s' => \$target_dir,
);

#$input_dir ||= "tempin";
#$output_dir ||= "tempout";
$input_dir ||= "/lustre/work1/ensembl/yuan/SARA/LRG/input_dir";
$output_dir ||= "/lustre/work1/ensembl/yuan/SARA/LRG/output_dir";
$target_dir ||= "/lustre/work1/ensembl/yuan/SARA/human/ref_seq_hash";
our $template_file = $template_file_name;
our $in_file = $in_file_name;
our $mapping_num = 1;
my %rec_seq;

# get registry and a gene adaptor
#Bio::EnsEMBL::Registry->load_all( $registry_file );
Bio::EnsEMBL::Registry->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');

my $species = 'human';
my $cdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $dbCore = $cdb;
my $sa = $dbCore->get_SliceAdaptor();

my $geneAd = Bio::EnsEMBL::Registry->get_adaptor($species ,'core', 'gene');
my $transAd = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'transcript');
my $protAd = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'protein');
print "genead is ",ref($geneAd),"\n";

open IN, $in_file or die "Could not read from input file ", $in_file, "\n";

my @cols = qw/lrg gene sequence cdna protein/;

# just a variable to hold the current node
my (%data,$current);

while(<IN>) {
    chomp;

    my @split = split /\s+/, $_;

    # check number of columns
    if(scalar @cols != scalar @split) {
		die "Incorrect number of columns in file ", $in_file, "\n";
    }


    # add data to a hash using @cols array for column names
    foreach my $col(@cols) {
		$data{$col} = shift @split;
    }

    # define the output XML file name
    my $out_file = (defined $out_file_stem ? $out_file_stem.'_' : '').$data{'lrg'}.'.xml';

    # create a root LRG object from the template file
    my $root = LRG::LRG::newFromFile($template_file, $out_file);
    
    
    
    
    ############################
    # FIXED ANNOTATION SECTION #
    ############################

    # get a ref to the fixed_annotation node
    my $fixed = $root->findOrAdd("fixed_annotation");

    # add ID and creation date
    $fixed->findOrAdd("id")->content($data{'lrg'});
    $fixed->findOrAdd("creation_date")->content($root->date);

    # change the URL
    $current = $fixed->findOrAdd("source/url");
    my $url = $current->content();
    $url =~ s/COL1A1/$data{'gene'}/;
    $current->content($url);

    # find and add genomic DNA sequence
    my $genomic_sequence = $root->pfetch($data{'sequence'});

    $fixed->findOrAdd("sequence")->content($genomic_sequence);

    my $trans_num = 1;

    my $trans_node = $fixed->addNode("transcript", {'name' => 't'.$trans_num++});

	my $trans;

	if($data{'cdna'} =~ /^ENS/) {
		$trans = $transAd->fetch_by_stable_id($data{'cdna'});
	}
	
	else {
		# get the transcript from EnsEMBL
		$trans = $transAd->fetch_all_by_external_name($data{'cdna'})->[0];
		
		# now refetch using the stable ID - hack for stupid bug
		$trans = $transAd->fetch_by_stable_id($trans->display_id);
	}
	
	print "Fetched transcript ", $trans->stable_id(), "\n";

    my $total_exon_seq = '';
    my $first = 1;
    my $first_exon_start;
    my $last_exon_end;

    # get all the exons in this transcript
    my @exons = @{$trans->get_all_Exons()};
	
	# get the cDNA and protein sequence from pfetch
	my $cdna_seq = $root->pfetch($data{'cdna'});
	my $prot_seq = $root->pfetch($data{'protein'});
	
    # add coding region node
    $current = $trans_node->addNode('coding_region');

    # add the translation node and the protein sequence
    $current->addNode('cdna')->addNode('sequence')->content($cdna_seq);
	
	# add the cdna node and cdna sequence
	$current->addNode('translation')->addNode('sequence')->content($prot_seq);
	
	# now work out where and what frame the AA sequence starts in the cDNA
	my $seq_obj = Bio::Seq->new(-name => 'cdna', -seq => $cdna_seq);
	
	my ($frame, $prot_start, $prot_end);
	
	foreach my $i(0..2) {
		$frame = $i;
		
		my $trans_prot_seq = $seq_obj->translate(undef, undef, $frame)->seq;
		
		if($trans_prot_seq =~ m/$prot_seq/g) {
			$prot_end = pos($trans_prot_seq);
			$prot_start = $prot_end - length($&);
			
			pos($trans_prot_seq) = 0;
			
			last;
		}
	}
	
	my $cdna_start = ($prot_start * 3) + $frame;
	my $cdna_end = $cdna_start + (3 * length($prot_seq));
	
	#print "CDNA $prot_start $prot_end $cdna_start $cdna_end $frame LENGTH ", length($prot_seq), "\n";
	
	my $prev_cds_end = 0;
	my $exon_number = 0;

    while(my $exon = shift @exons) {
		my ($start, $end);
	
		# get the sequence of this exon
		my $seq = $exon->seq->seq;
		$total_exon_seq .= $seq;
	
		# try and map the sequence to the genomic sequence
		if($genomic_sequence =~ m/$seq/g) {
			$end = pos($genomic_sequence);
			$start = $end - length($&) + 1;
		}
	
		# make sure to reset pos so the next m// searches from the beginning
		pos($genomic_sequence) = 0;
	
		# this bit of code allows for a 1bp mismatch between the exon sequence
		# and the genomic sequence if the whole didn't match
		if(!defined $start) {
			my $pos = 0;
	
			# take a copy of the sequence
			my $backup = $seq;
	
			# iterate through each position in the exon's sequence
			for my $pos(0..(length($seq)-1)) {
				# change that pos to '.' - this will match anything in the regexp
				substr($seq, $pos, 1) = '.';
		
				# now try the match again
				if($genomic_sequence =~ m/$seq/g) {
					$end = pos($genomic_sequence);
					$start = $end - length($&) + 1;
				}
				pos($genomic_sequence) = 0;
		
				# end the loop if matched
				last if defined $start;
		
				# otherwise recopy the sequence and continue the loop
				$seq = $backup;
			}
		}
		
		# store the coords of the first exon - needed later
		if($first) {
			$first = 0;
			$first_exon_start = $start;
		}
		
		# create exon node
		my $exon = $current->addNode('exon', {'lrg_number' => ++$exon_number});
		
		
		### EXON IN LRG COORDS
		$exon->addNode('lrg_start')->content($start);
		$exon->addNode('lrg_end')->content($end);
		
		
		### EXON IN CDNA COORDS
		my ($cds_start, $cds_end);
		
		if(defined $prev_cds_end) {
			$cds_start = $prev_cds_end + 1;
			$cds_end = $cds_start + ($end - $start);
		}
		
		else {
			$cds_start = ($start - $first_exon_start) + 1;
			$cds_end = $end - ($first_exon_start - 1);
		}
		
		$exon->addNode('cds_start')->content($cds_start);
		$exon->addNode('cds_end')->content($cds_end);
		
		
		### EXON IN AA COORDS
		my $aa_start = ((($cds_start - 1) - $cdna_start) / 3) + 1;
		my $aa_end = ((($cds_end - 1) - $cdna_start) / 3) + 1;
		
		# only record if both ends are not in either UTR
		if(($aa_end >= 1) && ($aa_start <= length($prot_seq))) {
			
			$aa_start = 1 if $aa_start < 1;
			$aa_end = length($prot_seq) if $aa_end > length($prot_seq);
			
			my $phase_start = sprintf("%.0f", 3 * ($aa_start - int($aa_start)));
			my $phase_end = sprintf("%.0f", 3 * ($aa_end - int($aa_end)));
			
			$aa_start = int($aa_start);
			$aa_end = int($aa_end);
			
			$exon->addNode('translation_start', {'phase' => $phase_start})->content($aa_start);
			$exon->addNode('translation_end', {'phase' => $phase_end})->content($aa_end);
		}
		
		# store the coords of the last exon - needed later
		$last_exon_end = $end;
		$prev_cds_end = $cds_end; 
    }
	
    # get coords for the coding region relative to the genomic sequence - relies on Ensembl
    my $temp = $current->addNode('cds_end');
	$temp->content($last_exon_end - (length($total_exon_seq) - $trans->cdna_coding_end));
	$temp->position(0);
	
    $temp = $current->addNode('cds_start');
	$temp->content($trans->cdna_coding_start + $first_exon_start - 1);
	$temp->position(0);
	
    
    
    ################################
    # UPDATABLE ANNOTATION SECTION #
    ################################
    
	$current = $root->findOrAdd('updatable_annotation/mapping');

	print "Starting mapping\n";
	
	substr($genomic_sequence, 50, 3) = "AAA";

	my $mapping = mapping($genomic_sequence);
	
	#print "DUMPY: ", Dumper($mapping), "\n";
	
	foreach my $pair(@{$mapping->type}) {
		print "@{$pair}\n";
		
		next unless $pair->[0] eq 'M';
		
		# switch chromosome start/end if on negative strand
		if($pair->[5] < 0 && $pair->[3] < $pair->[4]) {
			($pair->[3], $pair->[4]) = ($pair->[4], $pair->[3]);
		}
		
		$current->addEmptyNode(
			'align',
			{
				'chromosome' => $mapping->hseqname,
				'start' => $pair->[3],
				'end' => $pair->[4],
				'lrg_start' => $pair->[1],
				'lrg_end' => $pair->[2],
				'strand' => $pair->[5],
			}
		);
	}
	
	
	my @feature_nodes = @{get_annotations($mapping,$genomic_sequence)};
	$current = $root->findOrAdd('updatable_annotation/features');
	
    foreach my $feature(@feature_nodes) {
		$current->addExisting($feature);
	}
	
    # export
    $root->printAll();
}

close IN;


sub addExon {
    my ($current, $start, $end, $number) = @_;
    my $exon;

    if(defined $number) {
		$exon = $current->addNode('exon', {'lrg_number' => $number});
    }

    else {
		$exon = $current->addNode('exon');
		$exon->data({'lrg_number' => $exon->position + 1});
    }

    $exon->addNode('start')->content($start);
    $exon->addNode('end')->content($end);

    return $exon;
}

sub mapping {
  my $sequence = shift;

  my $seq_length = length($sequence);
	
  print "LENGTH $seq_length\n";
	
  # if we need to split it
  if($seq_length > 25000) {
    my $split_size = 10000;
		
    my $current_start = 0;
		
    my @maps = ();
    my @pairs = ();
	
    while($current_start < $seq_length) {
			
      # split_size is conservative; if this is the penultimate split, just do a bigger
      # one to make sure we don't try and map too small a sequence
      if($seq_length - $current_start < (2 * $split_size)) {
		$split_size = 2 * $split_size;
      }
			
      print "Mapping from $current_start to ", $current_start + $split_size, "\n";
      push @maps, mapping(substr($sequence, $current_start, $split_size));
			
      $current_start += $split_size;
    }

	my $prev_end = 0;
	my @dna_pairs;
        my $full_match =1;

	# get all pairs from all maps into arrays
	foreach my $map(@maps) {
		$full_match=0 if ($map->identical_matches==0);
		# sort by query start
		foreach my $pair(sort {$a->[2] <=> $b->[2]} @{$map->type}) {
			#print "full_match is $full_match and ",$map->identical_matches,'-',$pair->[0],'-',$pair->[1],'-',$pair->[2],'-',$pair-[3],'-',$pair->[4],'-',$pair->[5],"\n";
			# add the previous DNA end to each query coordinate
			# since we split up the query sequence
			$pair->[1] += $prev_end;
			$pair->[2] += $prev_end;
			
			# put the DNA pairs in a separate array
			if($pair->[0] eq 'DNA') {
				push @dna_pairs, $pair;
			}
			
			# put the match pairs in the @pairs array
			else {
				push @pairs, $pair;
			}
		}
		
		# get prev_end from the last added DNA pair
		$prev_end = $dna_pairs[-1]->[2];
	}
	
	# now join the pairs
	my @joined_pairs = (@{join_pairs(\@pairs)}, @{join_pairs(\@dna_pairs)});
	
	
	
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
    print "main_map is ",$main_map->seqname,'-',$main_map->hseqname,'-',$main_map->start,'-',$main_map->end,'-',$main_map->hstart,'-',$main_map->hend,'-',$main_map->identical_matches,"\n";
    # add the pairs onto the joined map
    $main_map->type(\@joined_pairs);

    #foreach my $pair (@{$main_map->type}) {
    #  print "checking DNA bit : ",$pair->[0],'-',$pair->[1],'-',$pair->[2],'-',$pair->[3],'-',$pair->[4],"\n";
    #}
    return $main_map;
  }
	
  else {
	
    my $name = "temp$$".$mapping_num++;
    #my $name = "temp23614".$mapping_num++;
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
    print "q_seq length is ",length($sequence),"\n";			
    my $seqobj = Bio::PrimarySeq->new(-id => $name, -seq => $rec_seq{$name});
    $rec_seq{$name} = $seqobj;
		
    bsub_ssaha_job($queue,$input_file,$output_file,$subject);

    my $call = "bsub -P ensembl-variation -K -w 'done($input_file\_ssaha_job)' -J waiting_process sleep 1"; #waits until all ssaha jobs have finished to continue
    system($call);
		
    my $mapping = parse_ssaha2_out($output_file);
		
    return $mapping->[0];
  }
}

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

sub bsub_ssaha_job {
  my ($queue, $input_file, $output_file, $subject) = @_;
	
  my $ssaha_command = "/nfs/acari/yuan/ensembl/src/ensembl-variation/scripts/ssahaSNP/ssaha2/ssaha2_v1.0.9_x86_64/ssaha2";
  $ssaha_command .= " -align 1 -kmer 12 -seeds 4 -cut 1000 -output vulgar -depth 10 -best 1 -save $subject $input_file";
  my $call = "bsub -J $input_file\_ssaha_job -P ensembl-variation $queue -e $output_dir/error_ssaha -o $output_file $ssaha_command";
	
  system ($call);
  print $call, "\n";
}

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
  # if we have more than one hit
  if(scalar @rec_find > 1 && $rec_find[0]->{'score'} == $rec_find[1]->{'score'}) {
    die "Mapped ",scalar @rec_find, "times\n";
  }
  my $mapping = make_feature_pair($rec_find[0]);
  return $mapping;
}

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

  print "$q_id,$q_start,$q_end,$q_strand,$t_id,$t_start,$t_end,$t_strand,$score and match $match\n";
  my ($f_q_start,$f_q_end) ;#feature start-end used to make feature_pairs
  if ($q_strand ==1) {
    ($f_q_start,$f_q_end) = ($q_start,$q_end);
  }
  else {
    ($f_q_start,$f_q_end) = ($q_end,$q_start);
  }

  my ($seq_region_name) = split /\-/, $t_id;

  my $slice = $sa->fetch_by_region('chromosome',$seq_region_name);

  my $q_seqobj = $rec_seq{$q_id};

  #warning that query_seq length != matched query_length
  if ($q_seqobj and length($q_seqobj->seq) != abs($q_end - $q_start) + 1){
    print "q_start is $q_start and q_end is $q_end and length is ",length($q_seqobj->seq),"\n";
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
      my $t_seq = substr($slice->seq,$t_start-1,$new_t_start-$t_start+1);
      #print "$q_start,$new_q_start,$t_start,$new_t_start and length of q_seq ",length($q_seq),'-',length($t_seq),"\n";
      #print "q_seq_is $q_seq\n";
      #print "t_seq is $t_seq\n";
      my $q_count = 1;
      my $t_count = 1;my ($q_seq);
            

      my %q_seqs = map {$q_count++,$_} split '', $q_seq;
      my %t_seqs = map {$t_count++,$_} split '', $t_seq;
      my ($sub_q_end);
      my $sub_q_start = $q_start;#initial sub_q_start
      my $sub_t_start = $t_start;
      foreach my $count (sort {$a<=>$b} keys %q_seqs) {
	if ($q_seqs{$count} !~ /$t_seqs{$count}/i) {
	  #if there is a mismatch, we need to record the base based in query sequence
          if ($q_strand==-1) {
            reverse_comp(\$q_seqs{$count});
            reverse_comp(\$t_seqs{$count});
          }

	  print "count is $count\n";
	  $full_match =0;
	  my $sub_t_end = $t_start + ($count -1) - 1;
	  if ($q_strand ==1) {
	    $sub_q_end = $q_start + ($count -1) - 1;
	  }
	  else {
	    $sub_q_end = $q_start - ($count -1) + 1;
	  }
	  if ($q_strand==1) {
	    $tmp_q_start = $sub_q_start;
	    $tmp_q_end = $sub_q_end;
	    $sub_q_start = $sub_q_end+2;#+1 cover SNP base, +1 cover next base
	  }
	  else {
	    $tmp_q_start = $sub_q_end;
	    $tmp_q_end = $sub_q_start;
	    $sub_q_start = $sub_q_end-2;
	  }
	
	  ($tmp_q_start,$tmp_q_end) = ($tmp_q_end,$tmp_q_start) if($tmp_q_end<$tmp_q_start);
	  push @pairs, [$type,$tmp_q_start,$tmp_q_end,$sub_t_start,$sub_t_end,$q_strand,uc($q_seqs{$count}),uc($t_seqs{$count})];

	  $sub_t_start = $sub_t_end+2;
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

  #foreach  my $pair (sort {$a->[3]<=>$b->[3]} @pairs) {
  #  print "pairs in ssaha_mapping are ",$pair->[0],'-',$pair->[1],'-',$pair->[2],'-',$pair->[3],'-',$pair->[4],"\n";
  #}
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
  my $lrg_name = $data{'lrg'};

  my ($q_start,$q_end,$t_start,$t_end,$q_strand);
  my $slice = $fp->slice;
  $q_start = $fp->start;
  $q_end = $fp->end;
  $t_start = $fp->hstart;
  $t_end = $fp->hend;
  $q_strand = $fp->strand;
  my $pairs = $fp->type;
  
  my $sub_slice = $slice->sub_Slice($t_start,$t_end);
  my $full_match = $fp->identical_matches;
  
  my @nodes;

  if ($full_match) {
    $sub_slice->seq_region_name($lrg_name);
    foreach my $gene (@{$sub_slice->get_all_Genes()}) {
      my @new_nodes = @{get_gene_annotation($gene)};
	  
      push @nodes, @new_nodes;
    }
  }
  else {
    my $csa = $dbCore->get_CoordSystemAdaptor();
    my $cs = $csa->fetch_all_by_name('LRG');
    my $cs_id = $cs->[0]->dbID();
    my $q_seq_region_id;
    my $q_seq_length = length($q_seq);
    my $t_seq_region_id = $sa->get_seq_region_id($slice);
    my $lrg_name_ref =$dbCore->dbc->db_handle->selectall_arrayref(qq{SELECT name from seq_region WHERE name="$lrg_name"});
    if (! $lrg_name_ref->[0][0]) {
      $dbCore->dbc->do(qq{INSERT INTO seq_region(name,coord_system_id,length)values("$lrg_name",$cs_id,$q_seq_length)});
      $q_seq_region_id = $dbCore->dbc->db_handle->{'mysql_insertid'};
    }
    else {
      my $q_seqid_ref = $dbCore->dbc->db_handle->selectall_arrayref(qq{SELECT seq_region_id from seq_region WHERE name="$lrg_name"});
      $q_seq_region_id = $q_seqid_ref->[0][0];
    }
    #we don't need dna sequence
    #my $dna_ref = $dbCore->dbc->db_handle->selectall_arrayref(qq{SELECT seq_region_id from dna WHERE seq_region_id = $q_seq_region_id});
    #if (! $dna_ref->[0][0]) {
    #  $dbCore->dbc->do(qq{INSERT INTO dna(seq_region_id,sequence)values($q_seq_region_id,"$q_seq")});
    #}
    
    foreach  my $pair (sort {$a->[3]<=>$b->[3]} @$pairs) {
      print "pairs are ",$pair->[0],'-',$pair->[1],'-',$pair->[2],'-',$pair->[3],'-',$pair->[4],'-', $pair->[5],'-',$pair->[6],'-',$pair->[7],"\n";
      
      if ($pair->[0] eq 'DNA') {
	if ($pair->[2]-$pair->[1] == $pair->[4]-$pair->[3]) {
	  $dbCore->dbc->do(qq{INSERT IGNORE INTO assembly(asm_seq_region_id,cmp_seq_region_id,asm_start,asm_end,cmp_start,cmp_end,ori)values($t_seq_region_id,$q_seq_region_id,$pair->[3],$pair->[4],$pair->[1],$pair->[2],$q_strand)});
	}
	else {
	  die("distance between query and target is not the same, there is a indel");
	}
      }
    }

    my $lrg_slice = $sa->fetch_by_region('LRG',"$lrg_name");

    my $min = 99999999;
    my $max = -9999999;
    my $chrom;
    my $strand;

    foreach my $segment (@{$lrg_slice->project('chromosome')}) {
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
    
      print "$from_start-$from_end  => $to_name $to_start-$to_end ($ori) \n";
    }

    foreach my $gene (@{$lrg_slice->get_all_Genes()}){
      print "gene_name is ",$gene->stable_id,"\n";
    }
    print "q_seq from lrg_slice is ",$lrg_slice->start,'-'.$lrg_slice->end,"\n";
    print "length of q_seq is ",length($q_seq), " and length of hseq is ", length($lrg_slice->seq),"\n";
    
    if ($lrg_slice->seq eq $q_seq) {
      print "hseq is same as q_seq\n";
    }
    else {
      print "hseq is different from q_seq\n";
      #print "q_seq is ",$q_seq,"\n";
      #print "hseq is ",$lrg_slice->seq,"\n";
    }

    my $ref_slice = $sa->fetch_by_region("chromosome",$chrom, $min, $max, $strand);
    my @genes = @{$ref_slice->get_all_Genes()};
	
    foreach my $gene (@genes) {
      next unless $gene->biotype eq 'protein_coding';
		
      print "before transfprm g start-end ",$gene->start,'-',$gene->end,"\n";
      #my $new_gene = $gene->transform('LRG');
      my $new_gene = $gene->transfer($lrg_slice);
      print "after transform g start-end ",$new_gene->start,'-',$new_gene->end,"\n" if $new_gene;
	  
      my @new_nodes = @{get_gene_annotation($new_gene)};
	  
      push @nodes, @new_nodes;
    }
  }
  
  return \@nodes;
}

sub get_gene_annotation {

	my $gene = shift;
	my ($current, $entry);
	my @nodes;
	
	# create the gene node
	my $gene_node = LRG::Node->new("gene", undef, {'name' => $gene->external_name, 'start' => $gene->start, 'end' => $gene->end});
	
	# get xrefs for the gene
	my $entries = $gene->get_all_DBEntries();

	print "\n\n>>> ", $gene->stable_id, " ", $gene->biotype, " ", $gene->status, "\n";

	while($entry = shift @$entries) {
		
		# get synonyms from HGNC entry
		if($entry->dbname eq 'HGNC') {
			foreach my $synonym(@{$entry->get_all_synonyms}) {
				$gene_node->addNode('synonym')->content($synonym);
			}
			
			$gene_node->addNode('note')->content($entry->description) if length($entry->description) > 1;
		}
		
		next unless $entry->dbname =~ /GI|RefSeq|HGNC$/;
		next if $entry->dbname =~ /RefSeq_peptide/;
		
		$gene_node->addExisting(xref($entry));
		
		print $entry->dbname, " ", $entry->primary_id, ".", $entry->version, " ", $entry->description, "\n";
	}

	# get the xrefs for the transcript
	my %ext = ();
	my %extdesc = ();
	
	foreach my $trans(@{$gene->get_all_Transcripts}) {
		
		#print "Gene/Trans ", $gene->stable_id, " ", $trans->stable_id, "\n";
		
		my $cds_node = LRG::Node->new("cds", undef, {'source' => 'Ensembl', 'codon_start' => $trans->coding_region_start, 'transcript_id' => $trans->stable_id});
		
		$entries = $trans->get_all_DBLinks();
		
		print "\n\n>>> ", $trans->stable_id, "\n";
		
		while($entry = shift @$entries) {
			
			next unless $entry->dbname =~ /GI|RefSeq|MIM_GENE|Entrez/;
			next if $entry->dbname =~ /RefSeq_peptide/;
			
			$cds_node->addExisting(xref($entry));
			
			print $entry->dbname, " ", $entry->primary_id, ".", $entry->version, " ", $entry->description, "\n";
			
			$ext{$entry->dbname} = $entry->primary_id;
			$extdesc{$entry->dbname} = $entry->description;
		}
		
		my $protein = $trans->translation;
		
		if($protein ne '') {
			my $prot_node = $cds_node->addNode('protein_product');
			
			$prot_node->addEmptyNode('protein_id', {'source' => 'Ensembl', 'accession' => $protein->stable_id});
			
			$entries = $protein->get_all_DBLinks();
			
			print "\n\n>>> ", $protein->stable_id, "\n";
		
			while($entry = shift @$entries) {
				
				next unless $entry->dbname =~ /RefSeq|Uniprot|CCDS/;
				
				if($entry->dbname eq 'RefSeq_peptide') {
					my $note_node = $prot_node->addNode('note');
					$note_node->content($entry->description);
					$note_node->moveTo(1);
				}
				
				$prot_node->addExisting(xref($entry));
				
				print $entry->dbname, " ", $entry->primary_id, ".", $entry->version, " ", $entry->description, "\n";
			}
		}
		
		push @nodes, $cds_node;
	}
	
	# finish the gene with xrefs
	$gene_node->addEmptyNode('db_xref', {'source' => 'Ensembl', 'accession' => $gene->stable_id});
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
