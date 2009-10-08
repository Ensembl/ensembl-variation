#!/software/bin/perl

use strict;
use lib '/nfs/users/nfs_w/wm2/Ensembl/CheckedOut/ensembl-variation/scripts/import/LRG/';
use lib '/nfs/users/nfs_i/ianl/ensembl-live/ensembl/modules/';
use Getopt::Long;
use LRG;
use LRGMapping qw(mapping get_annotations);
use Bio::Seq;
use Bio::EnsEMBL::Registry;

# set default options
my $output_file = 'make.lrg.xml';
my $template_file = 'template.xml';
my $registry_file = "ensembl.registry";
my ($config_file, $target_dir, $exon_file, $genomic_file, $cdna_file, $peptide_file, $lrg_id, $help);

usage() unless scalar @ARGV;

# get options from command line
GetOptions(
	   'xml_template=s' => \$template_file,
	   'registry_file=s' => \$registry_file,
	   'target_dir=s' => \$target_dir,
	   'exons=s' => \$exon_file,
	   'lrg=s' => \$genomic_file,
	   'cdna=s' => \$cdna_file,
	   'peptide=s' => \$peptide_file,
	   'out=s' => \$output_file,
	   'id=s' => \$lrg_id,
	   'help' => \$help,
	   'file=s' => \$config_file,
);

usage() if $help;

if(defined $config_file) {
	open IN, $config_file or warn "Could not read from config file $config_file\...attempting to continue\n";
	
	my (%config, $option, $value);
	
	while(<IN>) {
		chomp;
		
		my @data = split /\t/, $_;
		
		$option = shift @data;
		$value = join "\t", @data;
		
		$config{$option} = $value;
	}
	close IN;
	
	# now update the options
	$template_file ||= $config{'xml_template'};
	$registry_file ||= $config{'registry_file'};
	$target_dir ||= $config{'target_dir'};
	$exon_file ||= $config{'exons'};
	$cdna_file ||= $config{'cdna'};
	$genomic_file ||= $config{'lrg'};
	$peptide_file ||= $config{'peptide'};
	$lrg_id ||= $config{'id'};
	$output_file ||= $config{'out'};
}

# we need an ID
die "No LRG ID supplied (-id) or ID in wrong format (LRG_[0-9]+)\n" unless $lrg_id =~ /LRG_[0-9]+/;

$LRGMapping::lrg_name = $lrg_id;

# set the target dir
$LRGMapping::target_dir = $target_dir if defined $target_dir;

# get registry and a gene adaptor
Bio::EnsEMBL::Registry->load_all( $registry_file );
#Bio::EnsEMBL::Registry->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
$LRGMapping::dbCore = Bio::EnsEMBL::Registry->get_DBAdaptor('human','core');


my $root = LRG::LRG::newFromFile($template_file, $output_file);
my $fixed = $root->findOrAdd("fixed_annotation");

# add LRG ID
$fixed->findNode("id")->content($lrg_id);

# add creation date
$fixed->findNode("creation_date")->content($root->date);



## GENOMIC SEQUENCE
my $genomic_sequence;
	
if(defined $genomic_file) {

	# try to read from file
	if(open IN, $genomic_file) {
		while(<IN>) {
			next if /^>/;
			chomp;
			s/\s+//g;
			tr/acgt/ACGT/;
			$genomic_sequence .= $_;
		}
		close IN;
	}
	
	# if can't read from file, try to pfetch
	else {
		$genomic_sequence = $root->pfetch($genomic_file);
	}
	
	# check to see we have something
	if($genomic_sequence !~ /a|c|g|t/i) {
		die "Could not get genomic sequence, or error in genomic sequence\n";
	}
	
	$fixed->findOrAdd("sequence")->content($genomic_sequence);
}


## CDNA SEQUENCE
my $cdna_sequence;
my $trans_node;

if(defined $cdna_file) {

	# try to read from file
	if(open IN, $cdna_file) {
		while(<IN>) {
			next if /^>/;
			chomp;
			s/\s+//g;
			tr/acgt/ACGT/;
			$cdna_sequence .= $_;
		}
		close IN;
	}
	
	# if can't read from file, try to pfetch
	else {
		$cdna_sequence = $root->pfetch($cdna_file);
	}
	
	# check to see we have something
	if($cdna_sequence !~ /a|c|g|t/i) {
		die "Could not get cdna sequence, or error in cdna sequence\n";
	}
	
	$trans_node = $fixed->findOrAdd("transcript", {'name' => 't1'});
	$trans_node->addNode("cdna")->addNode("sequence")->content($cdna_sequence);
}


## PEPTIDE SEQUENCE
my $peptide_sequence;

if(defined $peptide_file) {

	# try to read from file
	if(open IN, $peptide_file) {
		while(<IN>) {
			next if /^>/;
			chomp;
			s/\s+//g;
			$peptide_sequence .= $_;
		}
		close IN;
	}
	
	# if can't read from file, try to pfetch
	else {
		$peptide_sequence = $root->pfetch($peptide_file);
	}
	
	# check to see we have something
	if($peptide_sequence !~ /[a-z]+/i) {
		die "Could not get peptide sequence, or error in peptide sequence\n";
	}
}


## EXONS
my @exons;

if(defined $exon_file) {
	open IN, $exon_file or die "Could not read from exon file $exon_file\n";
	
	while(<IN>) {
		chomp;
		
		my @data = split /\s+/, $_;
		
		my %exon = (
			'start' => $data[0],
			'end' => $data[1],
		);
		
		push @exons, \%exon;
	}
}

## FRAME CALCULATION
my $seq_obj = Bio::Seq->new(-name => 'cdna', -seq => $cdna_sequence);

my ($frame, $prot_start, $prot_end);

foreach my $i(0..2) {
	$frame = $i;
	
	my $trans_peptide_sequence = $seq_obj->translate(undef, undef, $frame)->seq;
	
	print "Frame $frame Comparing:\n", substr($trans_peptide_sequence, 0, 20), "\nTo:\n", substr($peptide_sequence, 0, 20), "\n\n";
	
	if($trans_peptide_sequence =~ m/$peptide_sequence/g) {
		$prot_end = pos($trans_peptide_sequence);
		$prot_start = $prot_end - length($&);
		
		pos($trans_peptide_sequence) = 0;
		
		last;
	}
}

die "Peptide sequence does not match translation of cDNA sequence\n" unless defined $prot_start;

# calculate coding start
my $coding_start = ($prot_start * 3) + $frame;

# calculate coding end - add 1 to include the STOP codon if necessary
my $coding_end = $coding_start + (3 * (length($peptide_sequence) + ($peptide_sequence =~ /\-$/ ? 0 : 1)));

print "CDS: $coding_start $coding_end\n";

my $coding_region_node = $trans_node->addNode("coding_region");
$coding_region_node->addNode("translation")->addNode("sequence")->content($peptide_sequence);

my $prev_cdna_end = 0;
my $exon_number = 0;
my $first = 1;
my $coding_end_lrg = $coding_end;
my ($first_exon_start, $last_exon_end, $start, $end, $cdna_start, $cdna_end);

while(my $exon = shift @exons) {
	($start, $end) = ($exon->{'start'}, $exon->{'end'});
	
	# store the coords of the first exon - needed later
	if($first) {
		$first = 0;
		$first_exon_start = $start;
	}
	
	# create exon node
	my $exon = $trans_node->addNode('exon');#, {'lrg_number' => ++$exon_number});
	
	
	### EXON IN LRG COORDS
	$exon->addEmptyNode('lrg_coords', {'start' => $start, 'end' => $end});
	
	
	### EXON IN CDNA COORDS
	if(defined $prev_cdna_end) {
		$cdna_start = $prev_cdna_end + 1;
		$cdna_end = $cdna_start + ($end - $start);
	}
	
	else {
		$cdna_start = ($start - $first_exon_start) + 1;
		$cdna_end = $end - ($first_exon_start - 1);
	}
	
	
	$exon->addEmptyNode('cdna_coords', {'start' => $cdna_start, 'end' => $cdna_end});
	
	
	### EXON IN AA COORDS
	my $aa_start = ((($cdna_start - 1) - $coding_start) / 3) + 1;
	my $aa_end = ((($cdna_end - 1) - $coding_start) / 3) + 1;
	
	# only record if both ends are not in either UTR
	if(($aa_end >= 1) && ($aa_start <= length($peptide_sequence))) {
		
		$aa_start = 1 if $aa_start < 1;
		$aa_end = length($peptide_sequence) if $aa_end > length($peptide_sequence);
		
		my $phase_start = sprintf("%.0f", 3 * ($aa_start - int($aa_start)));
		my $phase_end = sprintf("%.0f", 3 * ($aa_end - int($aa_end)));
		
		$aa_start = int($aa_start);
		$aa_end = int($aa_end);
		
		$exon->addEmptyNode(
			'peptide_coords',
			{
				'start' => $aa_start,
				'end' => $aa_end,
				#'start_phase' => $phase_start,
				#'end_phase' => ($aa_end == length($peptide_sequence) ? 2 : $phase_end)
			}
		);
		
		#$exon->addNode('intron_phase')->content($phase_end < 2 ? $phase_end + 1 : 0) if scalar @exons;
		$trans_node->addEmptyNode('intron', {'phase' => ($phase_end < 2 ? $phase_end + 1 : 0)}) if scalar @exons;
	}
	
	# update the coding_end coord with the intron length
	$coding_end_lrg += ($start - $last_exon_end) - 1 if $coding_end >= $cdna_start;
	
	#print "$start $end\t$cdna_start $cdna_end\t$coding_end $coding_end_lrg\n";
	
	# store the coords of the last exon - needed later
	$last_exon_end = $end;
	$prev_cdna_end = $cdna_end;
}

# add transcript start and end
$trans_node->addData({'start' => $first_exon_start, 'end' => $last_exon_end});

# add CDS start and end
$coding_region_node->addData({'start' => $first_exon_start + $coding_start, 'end' => $coding_end_lrg});



################################
# UPDATABLE ANNOTATION SECTION #
################################

my $current = $root->findOrAdd('updatable_annotation');
my $annotation_set_node = $current->findOrAdd('annotation_set');

# source
if(!$annotation_set_node->findNode('source')) {
	$current = $annotation_set_node->addNode('source');
	
	$current->addNode('name')->content('Ensembl');
	
	$current = $current->addNode('contact');
	$current->addNode('name')->content('Ensembl Variation');
	$current->addNode('address')->content('European Bioinformatics Institute');
	$current->addNode('email')->content('ensembl-variation@ebi.ac.uk');
}

$current = $annotation_set_node;

# update the modification date
$current->findOrAdd('modification_date')->content($root->date);

# other exon naming
$current->findOrAdd('other_exon_naming');

# add a mapping node
my $mapping_node = $current->addNode('mapping');

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
				'assembly' => 36,
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
}


# other exon naming
$annotation_set_node->findOrAdd('alternate_amino_acid_numbering');

# get annotations
my @feature_nodes = @{LRGMapping::get_annotations($mapping,$genomic_sequence)};

# add features node
$current = $annotation_set_node->findOrAdd('features');

# add the features returned from get_annotations
foreach my $feature(@feature_nodes) {
	$current->addExisting($feature);
}


$root->printAll();


sub usage() {
	
	print
	"Usage:
	
	-x || --xml_template	XML template file
	-r || --registry_file	registry file
	-t || --target_dir	target directory
	-e || --exons		exon coordinate file
	-l || --lrg		genomic sequence file
	-c || --cdna		cDNA sequence file
	-p || --peptide		peptide sequence file
	-o || --out		output file
	-i || --id		LRG identifier
	-h || --help		print this message\n";
	
	die;
}