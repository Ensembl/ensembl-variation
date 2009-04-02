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
    'in_file_name=s' => \$in_file_name,
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
our $template_file = "$input_dir/$template_file_name";
our $in_file = "$input_dir/$in_file_name";
our $mapping_num = 1;
my %rec_seq;

# get registry and a gene adaptor
Bio::EnsEMBL::Registry->load_all( $registry_file );
#Bio::EnsEMBL::Registry->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');

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
my $current;

while(<IN>) {
    chomp;

    my @split = split /\s+/, $_;

    # check number of columns
    if(scalar @cols != scalar @split) {
		die "Incorrect number of columns in file ", $in_file, "\n";
    }

    my %data;

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

    # used to store exon coords so we don't need to fetch them again later
    my @starts;
    my @ends;

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
	
		# add an exon node
		addExon($trans_node, $start, $end);
	
		# store the coords of the first exon - needed later
		if($first) {
			$first = 0;
			$first_exon_start = $start;
		}
	
		# store the coords of the last exon - needed later
		$last_exon_end = $end;
	
		push @starts, $start;
		push @ends, $end;
    }

    # add coding region node
    $current = $trans_node->addNode('coding_region');

    # get coords for the coding region relative to the genomic sequence - relies on Ensembl
    $current->addNode('cds_start')->content($trans->cdna_coding_start + $first_exon_start - 1);
    $current->addNode('cds_end')->content($last_exon_end - (length($total_exon_seq) - $trans->cdna_coding_end));

    # add the translation node and the protein sequence
    $current->addNode('translation')->addNode('sequence')->content($root->pfetch($data{'protein'}));

    # add cDNA node
    $current = $trans_node->addNode('cdna');

    # go through the exons again
    # this time we want coords relative to the coding region
    while(my $start = shift @starts) {
		my $end = shift @ends;
	
		$start -= $first_exon_start - 1;
		$end -= $first_exon_start - 1;
	
		addExon($current, $start, $end);
    }

    # add the cdna sequence
    $current = $current->addNode('sequence');
    $current->content($root->pfetch($data{'cdna'}));
    $current->position(0);
    
    
    
    
    ################################
    # UPDATABLE ANNOTATION SECTION #
    ################################
    
    $current = $root->findOrAdd('updatable_annotation/features');
    
    # gene
#     print "Fetching gene using ", $trans->stable_id, "\n";
    my $gene = $geneAd->fetch_by_transcript_stable_id($trans->stable_id);
    
    $current = $current->addNode('gene', {'start' => $first_exon_start, 'end' => $last_exon_end, 'name' => $gene->external_name()});
    
    my $entries = $gene->get_all_DBEntries();
    
    my $hgnc_entry;
    while($hgnc_entry = shift @$entries) {
    	last if $hgnc_entry->dbname eq 'HGNC';
    }
    
    foreach my $synonym(@{$hgnc_entry->get_all_synonyms}) {
    	$current->addNode('synonym')->content($synonym);
    }
    
    $current->addNode('note')->content($hgnc_entry->description) if length($hgnc_entry->description) > 1;

	# get the xrefs for the transcript
    my %ext = ();
    my %extdesc = ();
    
    foreach my $temptrans(@{$gene->get_all_Transcripts}) {
		$entries = $temptrans->get_all_DBLinks();
    
		while(my $entry = shift @$entries) {
			$ext{$entry->dbname} = $entry->primary_id;
			$extdesc{$entry->dbname} = $entry->description;
# 			print $gene->external_name, " ", $entry->dbname, " ", $entry->description, " ", $entry->primary_id, "\n";
		}
    }
    
    # finish the gene with MIM and HGNC xrefs
	$current->addEmptyNode('db_xref', {'source' => 'GeneID', 'accession' => $ext{'EntrezGene'}}) if defined $ext{'EntrezGene'};
    $current->addEmptyNode('db_xref', {'source' => 'HGNC', 'accession' => $ext{'HGNC'}}) if defined $ext{'HGNC'};
    $current->addEmptyNode('db_xref', {'source' => 'MIM', 'accession' => $ext{'MIM_GENE'}}) if defined $ext{'MIM_GENE'};
    
    $current = $current->parent;
    
    # now add the cds node
	$current = $current->addNode('cds', {'source' => 'RefSeq', 'transcript_id' => $ext{'RefSeq_dna'}, 'codon_start' => ($trans->cdna_coding_start + $first_exon_start - 1)});
	
	$current = $current->addNode('protein_product');
	
	$current->addEmptyNode('protein_id', {'source' => 'RefSeq', 'accession' =>  $data{'protein'}});
	$current->addNode('note')->content($extdesc{'RefSeq_peptide'}) if defined $extdesc{'RefSeq_peptide'};
	$current->addEmptyNode('db_xref', {'source' => 'CCDS', 'accession' => $ext{'CCDS'}}) if defined $ext{'CCDS'};
	$current->addEmptyNode('db_xref', {'source' => 'GeneID', 'accession' => $ext{'EntrezGene'}}) if defined $ext{'EntrezGene'};
	
	$current = $root->findOrAdd('updatable_annotation/mapping');

	my $mapping = mapping($genomic_sequence);
	
	foreach my $pair(@{$mapping->type}) {
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
	
    # export
    $root->printAll();
    get_annotations ($mapping);
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

	# get all pairs from all maps into arrays
	foreach my $map(@maps) {
		
		# sort by query start
		foreach my $pair(sort {$a->[2] <=> $b->[2]} @{$map->type}) {
			
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
	my $prev_pair;
	my @joined_pairs;
	my $pair;
	
	foreach $pair(@pairs) {
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
	
	# add the DNA pairs again
	push @joined_pairs, @dna_pairs;
		
	# join the maps
    my $main_map = shift @maps;
		
    my $prev_q_end = $main_map->end;
		
    foreach my $map(@maps) {
		$main_map->hend($map->hend);
		
		$main_map->end( $map->end + $prev_q_end);
		
		$prev_q_end = $main_map->end;
    }
	
	# add the pairs onto the joined map
	$main_map->type(\@joined_pairs);

    return $main_map;
  }
	
  else {
	
    my $name = "temp5863".$mapping_num++;
    #my $name = "temp_will";
    my $input_file_name = $name.'.fa';
		
    #open OUT, ">$input_dir/$name\.fa" or die $!;
	open OUT, ">$name\.fa" or die $!;
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
		
    #bsub_ssaha_job($queue,$input_file,$output_file,$subject);
		
    my $call = "bsub -P ensembl-variation -K -w 'done($input_file\_ssaha_job)' -J waiting_process sleep 1"; #waits until all ssaha jobs have finished to continue
    #system($call);
		
    my $mapping = parse_ssaha2_out($output_file);
		
    return $mapping->[0];
  }
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
      my $t_count = 1;

      my %q_seqs = map {$q_count++,$_} split '', $q_seq;
      my %t_seqs = map {$t_count++,$_} split '', $t_seq;
      my ($sub_q_end);
      my $sub_q_start = $q_start;#initial sub_q_start
      my $sub_t_start = $t_start;
      foreach my $count (sort {$a<=>$b} keys %q_seqs) {
	if ($q_seqs{$count} !~ /$t_seqs{$count}/i) {
	  #next;
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
	  push @pairs, [$type,$tmp_q_start,$tmp_q_end,$sub_t_start,$sub_t_end,$q_strand];

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
      if ($q_strand ==1) {
	$new_q_start = $q_start + $query_match_length + 1;
      }
      else {
	$new_q_start = $q_start - $query_match_length - 1;
      }
      $new_t_start = $t_start + $target_match_length +1;
      if ($q_strand ==-1) {
	$tmp_q_start = $new_q_start+1;
	$tmp_q_end = $q_start-1;
      }
      else {
	$tmp_q_start = $q_start+1;
	$tmp_q_end = $new_q_start-1;
      }
      ($tmp_q_start,$tmp_q_end) = ($tmp_q_end,$tmp_q_start) if($tmp_q_end<$tmp_q_start);
      #push @pairs, [$type,$tmp_q_start,$tmp_q_end,$t_start+1,$new_t_start-1,$q_strand];

      #print "in G, q_start is $q_start and q_end is $new_q_start and t_start is $t_start and t_end is $new_t_start\n";
    }

    $q_start = $new_q_start;
    $t_start = $new_t_start;
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
  my $lrg_name = 'LRG5';

  my ($q_start,$q_end,$t_start,$t_end,$q_strand);
  my $slice = $fp->slice;
  my $q_seqobj = $rec_seq{$fp->seqname};
  my $q_seq = $q_seqobj->seq;
  $q_start = $fp->start;
  $q_end = $fp->end;
  $t_start = $fp->hstart;
  $t_end = $fp->hend;
  $q_strand = $fp->strand;
  my $pairs = $fp->type;
  
  my $sub_slice = $slice->sub_Slice($t_start,$t_end);
  my $full_match = $fp->identical_matches;

  if ($full_match) {
    $sub_slice->seq_region_name($lrg_name);
    my @genes = @{$sub_slice->get_all_Genes()};
    my @transcripts = @{$genes[0]->get_all_Transcripts()};
    print "gene name is ", $genes[0]->stable_id,'-',$genes[0]->start,'-',$genes[0]->end,"\n" if defined $genes[0];
    print "trans name is ", $transcripts[0]->stable_id,'-',$transcripts[0]->start,'-',$transcripts[0]->end,"\n" if defined $transcripts[0];
    my $hnum_exons;

    foreach my $exon (@{$transcripts[0]->get_all_Exons }) {
      print "  ", $exon->stable_id,"\t",$exon->start,," ", $exon->end, "\n";
      $hnum_exons++;
    }
    print "There are ", $hnum_exons," exons\n";
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
    my $dna_ref = $dbCore->dbc->db_handle->selectall_arrayref(qq{SELECT seq_region_id from dna WHERE seq_region_id = $q_seq_region_id});
    if (! $dna_ref->[0][0]) {
      $dbCore->dbc->do(qq{INSERT INTO dna(seq_region_id,sequence)values($q_seq_region_id,"$q_seq")});
    }
    
    foreach  my $pair (sort {$a->[3]<=>$b->[3]} @$pairs) {
      print "pairs are ",$pair->[0],'-',$pair->[1],'-',$pair->[2],'-',$pair->[3],'-',$pair->[4],"\n";
      
      if ($pair->[0] eq 'DNA') {
	if ($pair->[2]-$pair->[1] == $pair->[4]-$pair->[3]) {
	  $dbCore->dbc->do(qq{INSERT IGNORE INTO assembly(asm_seq_region_id,cmp_seq_region_id,asm_start,asm_end,cmp_start,cmp_end,ori)values($t_seq_region_id,$q_seq_region_id,$pair->[3],$pair->[4],$pair->[1],$pair->[2],$q_strand)});
	}
	else {
	  die("distance between query and target is not the same, there is a indel");
	}
      }
    }

    my $hslice = $sa->fetch_by_region('LRG',"$lrg_name");
    print "q_seq from hslice is ",$hslice->start,'-'.$hslice->end,"\n";
    print "length of q_seq is ",length($q_seq), " and length of hseq is ", length($hslice->seq),"\n";
    
    if ($hslice->seq eq $q_seq) {
      print "hseq is same as q_seq\n";
    }
    else {
      print "hseq is different from q_seq\n";
      print "q_seq is ",$q_seq,"\n";
      print "hseq is ",$hslice->seq,"\n";
    }


    my $msc = Bio::EnsEMBL::MappedSliceContainer->new(
						      -SLICE => $hslice
						     );
    my $asa = $dbCore->get_AssemblySliceAdaptor();
    $msc->set_AssemblySliceAdaptor($asa);
    #$msc->attach_AssemblySlice('NCBI36');
    $msc->attach_LRG('LRG');
    
    foreach my $mapped_slice (@{$msc->get_all_MappedSlices()}){
      if ($mapped_slice->seq eq $q_seq) {
	print "mapped_seq is same as q_seq\n";
      }
      else {
	print "mapped_seq is different from q_seq\n";
	#print "mapped_seq is ",$mapped_slice->seq,"\n";
      }

      print "mapped_slice name is ",$mapped_slice->name,"\t",$mapped_slice->start,"\t",$mapped_slice->end,"\n";
      my $num_exons = 0;
      my $num_trans = 0;

      foreach my $exon (@{ $mapped_slice->get_all_Exons }) {#working here???
	print "  ", $exon->stable_id,"\t",$exon->start,," ", $exon->end, "\n";
	$num_exons++;
      }
      print "There are ", $num_exons," exons\n"; 

=head
#need mapped_slice->get_all_Exons to work for multi seq_region_name defined as LRG
	foreach my $gene (@{ $mapped_slice->get_all_Genes }) {#not working here
	  print "gene is ",ref($gene),"\n";
	  my @gene_dbentries = @{$gene->get_all_DBEntries()};
	  foreach my $dbe (@gene_dbentries) {
	    my @gene_syns = @{$dbe->get_all_synonyms()};
	    print "all_synonyms are @gene_syns\n";
	    print $dbe->db_display_name(),'-',$dbe->description(),'-',$gene->external_db(),"\n";
	  }
	  foreach my $transcript (@{$mapped_slice->get_all_Transcripts()}) {
	    print "Tanscript cdna_start_end  ", $transcript->stable_id,"\t",$transcript->cdna_coding_start,," ", $transcript->cdna_coding_end, "\n";
	    $num_trans++;
	    print "There are ", $num_trans," transcripts\n";
	    my $translation = $transcript->translation();
	    print "The protein sequence is ",$translation->seq,"\n";
	    my @transl_dbentries = @{$translation->get_all_DBEntries()};
	    foreach my $dbe (@transl_dbentries) {
	      print $dbe->db_display_name(),'-',$dbe->get_all_synonyms(),'-',$dbe->description(),'-',$translation->external_db,"\n";
	    }
	  }
	}
=cut
    }




    my @genes = @{$sub_slice->get_all_Genes()};
    foreach my $gene (@genes) {
      my $num_exons = 0;
      my $num_trans = 0;
      print "before transfprm g start-end ",$gene->start,'-',$gene->end,"\n";
      my $new_gene = $gene->transform('LRG');
      print "after transform g start-end ",$new_gene->start,'-',$new_gene->end,"\n" if $new_gene;
      my @transcripts = @{$new_gene->get_all_Transcripts()};
      print "trs start-end ",$transcripts[0]->start,'-',$transcripts[0]->end,"\n" if $transcripts[0];
      print "new_gene name is ", $new_gene->stable_id,'-',$new_gene->start,'-',$new_gene->end,"\n" if defined $new_gene;
      #my @genes = @{$hslice->get_all_Genes()};
      my @transcripts = @{$new_gene->get_all_Transcripts()} if $new_gene;
      #print "gene name is ", $genes[0]->stable_id,'-',$genes[0]->start,'-',$genes[0]->end,"\n" if defined $genes[0];
      print "trans name is ", $transcripts[0]->stable_id,'-',$transcripts[0]->start,'-',$transcripts[0]->end,"\n" if defined $transcripts[0];
      my $hnum_exons;

      foreach my $exon (@{$transcripts[0]->get_all_Exons }) {
	print "  ", $exon->stable_id,"\t",$exon->start,," ", $exon->end, "\n";
	$hnum_exons++;
      }
      print "There are ", $hnum_exons," exons\n";
      
      foreach my $transcript (@{$new_gene->get_all_Transcripts()}) {
	print "Tanscript cdna_start_end  ", $transcript->stable_id,"\t",$transcript->cdna_coding_start,," ", $transcript->cdna_coding_end, "\n";
	$num_trans++;
	print "There are ", $num_trans," transcripts\n";
	my $translation = $transcript->translation();
	print "The protein sequence is ",$translation->seq,"\n";
	my @transl_dbentries = @{$translation->get_all_DBEntries()};
	foreach my $dbe (@transl_dbentries) {
	  print $dbe->db_display_name(),'-',$dbe->get_all_synonyms(),'-',$dbe->description(),"\n";
	}
      }
    }
  }
}



