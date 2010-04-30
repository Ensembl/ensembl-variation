#!perl

use strict;
use Getopt::Long;
use LRG;
use LRGMapping qw(mapping get_annotations);
use Bio::Seq;
use Bio::EnsEMBL::Registry;
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

$LRGMapping::input_dir = $input_dir if defined $input_dir;
$LRGMapping::output_dir = $output_dir if defined $output_dir;
$LRGMapping::target_dir = $target_dir if defined $target_dir;

our $template_file = $template_file_name;
our $in_file = $in_file_name;
our $mapping_num = 1;
my %rec_seq;

# get registry and a gene adaptor
Bio::EnsEMBL::Registry->load_all( $registry_file );
#Bio::EnsEMBL::Registry->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');

my $species = 'human';
my $cdb = Bio::EnsEMBL::Registry->get_DBAdaptor('human','core');
my $dbCore = $cdb;

$LRGMapping::dbCore = $dbCore;

my $sa = $dbCore->get_SliceAdaptor();

my $geneAd = Bio::EnsEMBL::Registry->get_adaptor($species ,'core', 'gene');
my $transAd = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'transcript');
my $protAd = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'protein');

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
	
	$LRGMapping::lrg_name = $root->findNode('fixed_annotation/id')->content;

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
	
	# add the cdna node and cdna sequence
    $trans_node->addNode('cdna')->addNode('sequence')->content($cdna_seq);
	
    # add coding region node
    $current = $trans_node->addNode('coding_region');

    # add the translation node and the protein sequence
	$current->addNode('translation')->addNode('sequence')->content($prot_seq);
	
	
	$current = $trans_node;
	
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
		my $exon = $current->addNode('exon');#, {'lrg_number' => ++$exon_number});
		
		
		### EXON IN LRG COORDS
		$exon->addEmptyNode('lrg_coords', {'start' => $start, 'end' => $end});
		
		
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
		
		
		$exon->addEmptyNode('cdna_coords', {'start' => $cds_start, 'end' => $cds_end});
		
		
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
			
			$exon->addEmptyNode(
				'peptide_coords',
				{
					'start' => $aa_start,
					'end' => $aa_end,
					'start_phase' => $phase_start,
					'end_phase' => ($aa_end == length($prot_seq) ? 2 : $phase_end)
				}
			);
		}
		
		# store the coords of the last exon - needed later
		$last_exon_end = $end;
		$prev_cds_end = $cds_end; 
    }
	
	# add transcript start and end
	$trans_node->addData({'start' => $first_exon_start, 'end' => $last_exon_end});
	
	# add coding region start and end
	$trans_node->findNode('coding_region')->addData(
		{
			'start' => $trans->cdna_coding_start + $first_exon_start - 1,
			'end' => ($last_exon_end - (length($total_exon_seq) - $trans->cdna_coding_end)),
		}
	);
	
    
    
    ################################
    # UPDATABLE ANNOTATION SECTION #
    ################################
    
	# add a mapping node
	my $mapping_node = $root->findOrAdd('updatable_annotation/annotation_set/mapping');
	
	# update the modification date
	$root->findNode('updatable_annotation/annotation_set/modification_date')->content($root->date);
	
	# run the mapping sub-routine
	my $mapping = LRGMapping::mapping($genomic_sequence);
	my $current = $mapping_node->addNode('mapping_span');
	
	# create three arrays to hold the different pair types
	my (@dna_pairs, @match_pairs, @mis_pairs);
	
	# sort the pairs into the arrays
	foreach my $pair(sort {$a->[1] <=> $b->[1]} @{$mapping->type}) {
		
		# DNA pair (should be only 1)
		if($pair->[0] eq 'DNA') {
			push @dna_pairs, $pair;
		}
		
		# mismatch pairs - have bases in [6] and [7]
		elsif($pair->[6].$pair->[7]) {
			push @mis_pairs, $pair;
		}
		
		# normal align pairs
		else {
			push @match_pairs, $pair
		}
	}
	
	# go through DNA, then align, then mismatch pairs
	foreach my $pair(@dna_pairs, @match_pairs, @mis_pairs) {
		
		my ($l_s, $l_e, $c_s, $c_e);
		
		($l_s, $l_e, $c_s, $c_e) = ($pair->[1], $pair->[2], $pair->[3], $pair->[4]);
		
		# switch chromosome start/end if on negative strand
		if($pair->[5] < 0 && $c_s < $c_e) {
			($c_s, $c_e) = ($c_e, $c_s);
		}
		
		# switch LRG start/end if needed - start should _always_ be less than end since LRG-centric
		($l_s, $l_e) = ($l_e, $l_s) if $l_s > $l_e;
		
		# DNA pair - whole alignment
		if($pair->[0] eq 'DNA') {
			
			# this gets added as attributes to the mapping node
			$mapping_node->addData(
				{
					'chr_name' => $mapping->hseqname,
					'chr_start' => ($c_s < $c_e ? $c_s : $c_e),
					'chr_end' => ($c_s < $c_e ? $c_e : $c_s)
				}
			);
			
			$current->addData(
				{
					'start' => $c_s,
					'end' => $c_e,
					'lrg_start' => $l_s,
					'lrg_end' => $l_e,
					'strand' => $pair->[5],
				}
			)
		}
		
		# mismatch pair
		elsif($pair->[6].$pair->[7] =~ /.+/) {
			$current->addEmptyNode(
				'diff',
				{
					'type' => 'mismatch',
					'start' => $c_s,
					'end' => $c_e,
					'lrg_start' => $l_s,
					'lrg_end' => $l_e,
					'lrg_sequence' => $pair->[6],
					'genomic_sequence' => $pair->[7],
				}
			);
		}
		
		# match pair
		#else {
		#	$current->addEmptyNode(
		#		'align',
		#		{
		#			'chromosome' => $mapping->hseqname,
		#			'start' => $c_s,
		#			'end' => $c_e,
		#			'lrg_start' => $l_s,
		#			'lrg_end' => $l_e,
		#			'strand' => $pair->[5],
		#		}
		#	);
		#}
	}
	
	# get annotations
	my @feature_nodes = @{LRGMapping::get_annotations($mapping,$genomic_sequence)};
	
	# add features node
	$current = $root->findOrAdd('updatable_annotation/annotation_set/features');
	
	# add the features returned from get_annotations
    foreach my $feature(@feature_nodes) {
		$current->addExisting($feature);
	}
	
    # export
    $root->printAll();
}

close IN;
