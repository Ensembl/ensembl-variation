#!perl -w

use strict;

use Getopt::Long;
use List::Util qw(min max);
use LRG::LRG;
use LRG::LRGMapping qw(mapping get_annotations);
use Bio::Seq;
use Bio::EnsEMBL::Registry;

# This variable should indicate the name of the most recent assembly available
my $CURRENT_ASSEMBLY = 'GRCh37';
my $CURRENT_SCHEMA_VERSION = '1.6';

my ($template_file, $output_file, $registry_file, $config_file, $target_dir, $exon_file, $genomic_file, $cdna_file, $peptide_file, $lrg_id, $help, $skip_fixed, $skip_updatable, $skip_unbranded, $skip_transcript_matching, $use_existing_mapping, $replace_annotations, $skip_host_check);
my $tmpdir;
my $ssaha2_bin;
my $exonerate_bin;
&usage() if (!scalar(@ARGV));

# get options from command line
GetOptions(
	   'tmpdir=s' => \$tmpdir, 
	   'ssaha2=s' => \$ssaha2_bin, 
	   'exonerate=s' => \$exonerate_bin, 
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
	   'skip_fixed!' => \$skip_fixed,
	   'skip_updatable!' => \$skip_updatable,
	   'skip_unbranded!' => \$skip_unbranded,
	   'skip_transcript_matching!' => \$skip_transcript_matching,
	   'use_existing_mapping!' => \$use_existing_mapping,
	   'replace_annotations!' => \$replace_annotations,
	   'skip_host_check!' => \$skip_host_check,
);

&usage() if $help;

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

# Revert unspecified options to default values
$template_file ||= 'template.xml';
$registry_file ||= 'ensembl.registry';
$output_file ||= 'make.lrg.xml';
	

$LRGMapping::current_assembly = $CURRENT_ASSEMBLY;
$LRGMapping::input_dir = $tmpdir;
$LRGMapping::output_dir = $tmpdir;
$LRGMapping::SSAHA_BIN = $ssaha2_bin;
$LRGMapping::EXONERATE_BIN = $exonerate_bin;

# we need an ID
die "No LRG ID supplied (-id) or ID in wrong format (LRG_[0-9]+)\n" unless $lrg_id =~ /LRG_[0-9]+|tmp_.+/;
$LRGMapping::lrg_name = $lrg_id;

# This must correspond to the assembly used to map the LRG and where annotations are fetched
my $assembly = $CURRENT_ASSEMBLY;

if (!$skip_updatable) {
# set the target dir
	$LRGMapping::target_dir = $target_dir if defined $target_dir;

# get registry and a gene adaptor
	$LRGMapping::registry_file = $registry_file;
	Bio::EnsEMBL::Registry->load_all( $registry_file );
	$LRGMapping::dbCore_ro = Bio::EnsEMBL::Registry->get_DBAdaptor('human','core_ro');
	$LRGMapping::dbCore_rw = Bio::EnsEMBL::Registry->get_DBAdaptor('human','core_rw');
	$LRGMapping::dbFuncGen = Bio::EnsEMBL::Registry->get_DBAdaptor('human','funcgen');
	my $host = $LRGMapping::dbCore_rw->dbc->host();
	if ($host !~ m/variation/ && !$skip_host_check) {
		die('Host is ' . $host . '! Changes will be written to the database, make sure you want to use the database on this host. If so, skip this check by using -skip_host_check on command line');
	}

# Get the assembly used to fetch annotations
#	my $csa = Bio::EnsEMBL::Registry->get_adaptor('Human','Core','CoordSystem');
#	$assembly = $csa->fetch_top_level()->version();
}

my $root = LRG::LRG::newFromFile($template_file, $output_file);

# Set the schema version
$root->findNode('lrg')->data({'schema_version' => $CURRENT_SCHEMA_VERSION});

################################
#   FIXED ANNOTATION SECTION   #
################################	
	
create_fixed_annotation($root,$lrg_id,$genomic_file,$exon_file,$cdna_file,$peptide_file) unless $skip_fixed;

################################
# UPDATABLE ANNOTATION SECTION #
################################

if (!$skip_updatable) {

	# This will always be called to update the modification date and Ensembl version info. If $use_existing_mapping is true and $replace_annotations is false, immediately returns
	create_updatable_annotation($root,$use_existing_mapping,$replace_annotations,$LRGMapping::dbCore_rw->get_SliceAdaptor());

	# Move annotations that should be "unbranded" to the LRG section and do consistency checking between these NCBI and Ensembl annotations
	move_to_unbranded($root) unless $skip_unbranded;
	
	align_updatable_to_fixed_transcripts($root,$LRGMapping::dbCore_ro->get_TranscriptAdaptor()) unless $skip_transcript_matching;
	
	# Find transcripts in the updatable section that correspond to transcripts in the fixed section (only for Ensembl annotations for now)
	match_fixed_annotation_transcripts($root,$LRGMapping::dbCore_ro->get_TranscriptAdaptor()) unless $skip_transcript_matching;

}

# Dump XML to output_file
$root->printAll();



sub create_fixed_annotation {
	my $root = shift;
	my $lrg_id = shift;
	my $genomic_file = shift;
	my $exon_file = shift;
	my $cdna_file = shift;
	my $peptide_file = shift;
	
	my $fixed = $root->findOrAdd("fixed_annotation");

# add LRG ID
	$fixed->findOrAdd("id")->content($lrg_id);

# add creation date
	$fixed->findOrAdd("creation_date")->content($root->date);

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
		if($genomic_sequence !~ /^[acgt]+$/i) {
			die "Could not get genomic sequence, or error in genomic sequence\n";
		}
	
		$fixed->findOrAdd("sequence")->content($genomic_sequence);
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
	}
	else {
## EXTRACT cDNA SEQUENCE FROM GENOME SEQUENCE (If not specified)
	  die ("Exon coordinates must be specified in order to infer cDNA seuqnece!") unless (scalar(@exons));
	  $cdna_sequence = make_cdna($genomic_sequence,\@exons);
	}
# check to see we have something
	if($cdna_sequence !~ /^[acgt]+$/i) {
		die "Could not get cdna sequence, or error in cdna sequence\n";
	}
	
	$trans_node = $fixed->findOrAdd("transcript", {'name' => 't1'});
	$trans_node->addNode("cdna")->addNode("sequence")->content($cdna_sequence);

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

	my $coding_start;
	my $coding_end;

## FRAME CALCULATION
	{
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
		$coding_start = ($prot_start * 3) + $frame;

# calculate coding end - add 1 to include the STOP codon if necessary
		$coding_end = $coding_start + (3 * (length($peptide_sequence) + ($peptide_sequence =~ /\-$/ ? 0 : 1)));
	}	

	print "CDS: $coding_start $coding_end\n";

	my $coding_region_node = $trans_node->addNode("coding_region");
	$coding_region_node->addNode("translation")->addNode("sequence")->content($peptide_sequence);

	my $prev_cdna_end = 0;
	my $exon_number = 0;
	my $first = 1;
	my $coding_end_lrg = $coding_end;
	my $coding_start_lrg = $coding_start;
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
	
# update the coding_start/end coord with the intron length
		$coding_start_lrg += ($start - $last_exon_end) - 1 if $coding_start > $cdna_start;
		$coding_end_lrg += ($start - $last_exon_end) - 1 if $coding_end >= $cdna_start;
	
		print "$start $end\t$cdna_start $cdna_end\t$coding_start $coding_start_lrg\t$coding_end $coding_end_lrg\n";
	
	# store the coords of the last exon - needed later
		$last_exon_end = $end;
		$prev_cdna_end = $cdna_end;
	}

# add transcript start and end
	$trans_node->addData({'start' => $first_exon_start, 'end' => $last_exon_end});

# add CDS start and end
	$coding_region_node->addData({'start' => $coding_start_lrg, 'end' => $coding_end_lrg});
}

#ÊThis routine will create the updatable annotation section if necessary and the Ensembl annotation_set.
# If second parameter is true, no mapping will be performed but the pre-existing one will be used. If it is false, a new mapping to the db assembly will be made
#ÊIf the third parameter is true, new annotations will be fetched from the database and replace the existing ones
sub create_updatable_annotation {
	my $root = shift;
	my $use_existing_mapping = shift;
	my $replace_annotations = shift;
	my $slice_adaptor = shift;
	
	# If no annotations will be added and no new mapping will be performed, we should return here
	return if ($use_existing_mapping && !$replace_annotations);
	
	my $genomic_sequence = $root->findNode('fixed_annotation/sequence')->content();
	
	my $current;
	
	#ÊAdd the updatable_annotation element if not present
	my $annotation_node = $root->findOrAdd('updatable_annotation');
	
	#ÊGet the Ensembl annotation_set if it's present
	my $annotation_set_nodes = $annotation_node->findNodeArray('annotation_set');
	my $annotation_set_ensembl;
	foreach my $asn (@{$annotation_set_nodes}) {
		my $src = $asn->findNode('source');
		if ($src && $src->findNode('name')->content() eq 'Ensembl') {
			$annotation_set_ensembl = $asn;
		}
	}
	# If not, create it
	if(!$annotation_set_ensembl) {
		$annotation_set_ensembl = $annotation_node->addNode('annotation_set');
		$current = $annotation_set_ensembl->addNode('source');	
		$current->addNode('name')->content('Ensembl');
		$current->addNode('url')->content('http://www.ensembl.org/');
		$current = $current->addNode('contact');
		$current->addNode('name')->content('Ensembl Variation');
		$current->addNode('address')->content('European Bioinformatics Institute');
		$current->addNode('email')->content('helpdesk@ensembl.org');
	}
	$current = $annotation_set_ensembl;
	
# Update the comment field to indicate which Ensembl release was used
	$annotation_set_ensembl->findOrAdd('comment')->content('Annotation is based on Ensembl release ' . $LRGMapping::dbCore_ro->get_MetaContainerAdaptor()->get_schema_version());
	
# Update the modification date
	$annotation_set_ensembl->findOrAdd('modification_date')->content($root->date);

# Add other exon naming element
	$annotation_set_ensembl->findOrAdd('other_exon_naming');

# Get the mapping node corresponding to the db assembly if it exists
	my $mapping_node = $annotation_node->findNode('mapping',{'assembly' => $assembly});
	my $chr_id;
	
# Run the mapping sub-routine. unless we want to use the mapping already made, this should allow for multiple mappings. Only the most_recent mapping should be used for getting annotations though.
	my $mapping;
	#ÊIf a new mapping should be made
	if (!$use_existing_mapping) {
		
		# Clear the old mapping from db
		LRGMapping::clear_mapping($LRGMapping::lrg_name,'lrg');
		
		# Run ssaha2 to map and the result is a reference to a hash
		$mapping = LRGMapping::mapping($genomic_sequence);
		my $pairs = $mapping->{'pairs'};
		$chr_id = $mapping->{'chr_name'};
		
		# Warn if the mapping could not be found over the entire LRG region
		warn("*** WARNING *** Could not map the entire LRG region to the $CURRENT_ASSEMBLY assembly!") if ($mapping->{'lrg_start'} > 1 || $mapping->{'lrg_end'} < length($genomic_sequence));
		
		# Create a mapping node from the pairs array
		my $new_node = LRGMapping::pairs_2_mapping($pairs,$assembly,$chr_id,($assembly eq $CURRENT_ASSEMBLY));
		
		# Remove the old one and replce it with the newly created one
		$mapping_node->remove() if (defined($mapping_node));
		$annotation_set_ensembl->addExisting($new_node);
	}
	#ÊElse, get the already existing one, corresponding to the db assembly and parse it
	else {
		# Die if no existing mapping could be found
		die("Could not find any pre-existing genomic mapping for " . $LRGMapping::lrg_name . " to $assembly assembly") unless (defined($mapping_node));
		
		# Else, parse the mapping
		$chr_id = $mapping_node->data->{'chr_name'};
		my $chr_start = $mapping_node->data->{'chr_start'};
		my $chr_end = $mapping_node->data->{'chr_end'};
		my $chr_slice = $slice_adaptor->fetch_by_region('chromosome',$chr_id);
		my $chr_seq = $chr_slice->subseq($chr_start,$chr_end);
		$mapping = LRGMapping::mapping_2_pairs($mapping_node,$genomic_sequence,$chr_seq);
	}
	
	# Add alternative amino acid numbering element
	$annotation_set_ensembl->findOrAdd('alternate_amino_acid_numbering');
	
	#ÊIf we will not fetch annotations, we're done
	return if (!$replace_annotations);
	
	# Get annotations
	my @feature_nodes = @{LRGMapping::get_annotations($LRGMapping::lrg_name,'lrg',$chr_id,length($genomic_sequence),$mapping)};

	# add features node
	my $feat_node = $annotation_set_ensembl->findNode('features');
	$feat_node->remove() if(defined($feat_node));
	$feat_node = LRG::Node->new('features');
	
	my $fixed_transcripts = $root->findNode('lrg/fixed_annotation')->findNodeArray('transcript');

# add the features returned from get_annotations
	foreach my $feature (@feature_nodes) {
# Need to do an extra check to remove any lrg_gene_name tags that might have been set within genes that are not the ones that this LRG was created for. Determine this by looking for a transcript that has the same cds_start and cds_end as any of the transcripts in the fixed section
		my $lrgs = $feature->findNodeArray('lrg_gene_name');
		foreach my $tag (@{$lrgs}) {
			my $lrg_gene = 0;
			my $up_transcripts = $feature->findNodeArray('transcript');
			if (defined($up_transcripts)) {
				foreach my $transcript (@{$up_transcripts}) {
					my $protein = $transcript->findNode('protein_product');
					next unless (defined($protein));
					my $up_cds_start = $protein->{'data'}{'cds_start'};
					my $up_cds_end = $protein->{'data'}{'cds_end'};
					foreach my $fix_transcript (@{$fixed_transcripts}) {
						my $cds = $fix_transcript->findNode('coding_region');
						my $fix_cds_start = $cds->{'data'}{'start'};
						my $fix_cds_end = $cds->{'data'}{'end'};
						if ($up_cds_start == $fix_cds_start && $up_cds_end == $fix_cds_end) {
							$lrg_gene = 1;
							last;
						}
					}
					last if ($lrg_gene);
				}
			}
			# Move the lrg_gene_name tag to the annotation_set node
			$tag->remove();
			if ($lrg_gene) {
				#ÊRemove any existing lrg_gene_name tag
				my $old_tag = $annotation_set_ensembl->findNode('lrg_gene_name');
				$old_tag->remove() if (defined($old_tag));
				$annotation_set_ensembl->addExisting($tag);
			}
		}
		$feat_node->addExisting($feature);
	}
	$annotation_set_ensembl->addExisting($feat_node);
}

sub get_annotation_set {
	my $root = shift;
	my $source = shift;
	
	# Get the annotation sets in the updatable section and loop through them to get the desired section
	my $annotation_sets = $root->findNodeArray('lrg/updatable_annotation/annotation_set');
	my $desired_set = undef;
	foreach my $annotation_set (@{$annotation_sets}) {
		next if ($annotation_set->findNode('source/name')->content() ne $source);
		$desired_set = $annotation_set;
		last;
	}
	
	return $desired_set;
}

sub get_gene {
	my $annotation_set = shift;
	my $gene_symbol = shift;
	
	# Get the gene corresponding to the gene symbol
	my $gene = $annotation_set->findNode('features/gene',{'symbol' => $gene_symbol});
	
	return $gene;
}

sub align_updatable_to_fixed_transcripts {
	my $root = shift;
	my $transcript_adaptor = shift;
	
	# Get the HGNC gene symbol for this LRG record
	my $hgnc_symbol_node = get_annotation_set($root,'LRG')->findNode('lrg_gene_name',{'source' => 'HGNC'});
	my $hgnc_symbol;
	if (!defined($hgnc_symbol_node)) {
		warn("Could not find lrg_gene_name node with source 'HGNC' in LRG-branded annotation set");
		return;
	}
	my $hgnc_symbol = $hgnc_symbol_node->content();
	
	# Get the transcripts in the fixed annotation section		
	my $fixed_transcripts = $root->findNode('lrg/fixed_annotation')->findNodeArray('transcript');
	
	# Get the Ensembl annotation set
	my $annotation_set = get_annotation_set($root,'Ensembl') or die("Could not get the Ensembl annotation_set");
	
	# Get the gene corresponding to the LRG record
	my $gene = get_gene($annotation_set,$hgnc_symbol);
	my $q_transcripts = $gene->findNodeArray('transcript',{'source' => 'Ensembl'});
	
	# Loop over the transcripts and align each one to each of the fixed transcripts
	foreach my $q_transcript (@{$q_transcripts}) {
		# Get the transcript object from Ensembl API
		my $stable_id = $q_transcript->data()->{'transcript_id'};
		my $ens_transcript = $transcript_adaptor->fetch_by_stable_id($stable_id) or warn("Could not fetch Ensembl transcript with stable id $stable_id");
		my $q_seq = $ens_transcript->spliced_seq();
		
		# Loop over the fixed transcripts
		foreach my $t_transcript (@{$fixed_transcripts}) {
			my $t_seq = $t_transcript->findNode('cdna/sequence')->content();
			
			my $o_file = LRGMapping::exonerate_align($q_seq,$t_seq);
			my $data = LRGMapping::parse_exonerate($o_file);
			#unlink($o_file);
			$data->{'chr_name'} = $t_transcript->data()->{'name'};
			$data->{'lrg_id'} = $stable_id;
			my $alignment = LRGMapping::pairs_2_alignment($data);
			$q_transcript->addExisting($alignment);
		}
	}
}
	
sub match_fixed_annotation_transcripts {
	my $root = shift;
	my $transcript_adaptor = shift;
	
# Get the transcripts in the fixed annotation section		
	my $fixed_transcripts = $root->findNode('lrg/fixed_annotation')->findNodeArray('transcript');
	
# Get the annotation sets in the updatable section and loop through them to get the Ensembl and LRG sections	
	my $annotation_sets = $root->findNode('lrg/updatable_annotation')->findNodeArray('annotation_set');	
	my $ensembl_annotation;
	my $ncbi_annotation;
	my $lrg_annotation;
	foreach my $as (@{$annotation_sets}) {
		if ($as->findNode('source/name')->content() eq 'Ensembl') {
			$ensembl_annotation = $as;
		}
		elsif ($as->findNode('source/name')->content() eq 'LRG') {
			$lrg_annotation = $as;
		}
		elsif ($as->findNode('source/name')->content() =~ m/NCBI/i) {
			$ncbi_annotation = $as;
		}
	}
	die('Could not find Ensembl annotation set!') unless (defined($ensembl_annotation));
#	die('Could not find LRG annotation set!') unless (defined($lrg_annotation));

# Get the HGNC symbol so we can look at the transcripts for the correct gene	
#	my $hgnc_symbol = $lrg_annotation->findNode('lrg_gene_name')->content();

# Get the annotated NCBI genes
	my $feature_node = $ncbi_annotation->findNode('features') unless (!defined($ncbi_annotation));
	my $ncbi_genes = $feature_node->findNodeArray('gene') unless (!defined($feature_node));

# Get the annotated Ensembl genes	
	$feature_node = undef;
	$feature_node = $ensembl_annotation->findNode('features') unless (!defined($ensembl_annotation));
	my $ensembl_genes = $feature_node->findNodeArray('gene') unless (!defined($feature_node));

# Keep track of which transcript has been matched to which fixed_id
	my %matched;
	
# Loop over the transcripts in the fixed annotation and see if the exon coordinates match for any of the updatabale transcripts
	foreach my $fixed_transcript (@{$fixed_transcripts}) {
		my $fixed_id = $fixed_transcript->data()->{'name'};
		my $fixed_exons = $fixed_transcript->findNodeArray('exon');
	
# Loop over the updatable genes
		foreach my $ensembl_gene (@{$ensembl_genes}) {
		
# Get the annotated transcripts		
			my $ensembl_transcripts = $ensembl_gene->findNodeArray('transcript');
			foreach my $ensembl_transcript (@{$ensembl_transcripts}) {
# Fetch the transcript from the Core database using the annotated stable id
				my $core_transcript = $transcript_adaptor->fetch_by_stable_id($ensembl_transcript->{'data'}{'transcript_id'});
# Get the exons for the transcript from the Core database			
				my $core_exons = $core_transcript->get_all_Exons();
# Check if the number of exons is equal
                                next unless (scalar(@{$core_exons}) == scalar(@{$fixed_exons}));
# Loop until no more exons or the coordinates don't match
				my $i = 0;
				while ($i < scalar(@{$core_exons}) && $core_exons->[$i]->cdna_start($core_transcript) == $fixed_exons->[$i]->findNode('cdna_coords')->{'data'}{'start'} && $core_exons->[$i]->cdna_end($core_transcript) == $fixed_exons->[$i]->findNode('cdna_coords')->{'data'}{'end'}) {
					$i++;
				}
				if ($i == scalar(@{$core_exons})) {
# If all exon coordinates match, store the transcript as the match for this fixed id
					$matched{$fixed_id} = $ensembl_transcript;
					last;
				}
			}
			last if (exists($matched{$fixed_id}));
		}
		
# If no match could be made, check if a cross-referenced NCBI transcript has a fixed id assigned and in that case use that.
		if (!exists($matched{$fixed_id})) {
			
#ÊGet the NCBI transcript for this fixed_id
			my $ncbi_transcript;
			foreach my $ncbi_gene (@{$ncbi_genes}) {
				$ncbi_transcript = $ncbi_gene->findNode('transcript',{'fixed_id' => $fixed_id}) or next;
				last;
			}
			
# Loop over the updatable genes
			if (defined($ncbi_transcript)) {
				foreach my $ensembl_gene (@{$ensembl_genes}) {
		
# Get the annotated transcripts		
					my $ensembl_transcripts = $ensembl_gene->findNodeArray('transcript');
					foreach my $ensembl_transcript (@{$ensembl_transcripts}) {
						# Skip to next if the transcript is not xref:d to the NCBI transcript for this fixed id
						my $refseq_xref = $ensembl_transcript->findNode('db_xref',{'source' => 'RefSeq', 'accession' => $ncbi_transcript->data()->{'transcript_id'}}) or next;
# It is possible for several Ensembl transcripts to be xref:d to the same NCBI transcript. So in that case, pick the longest one
						if (!exists($matched{$fixed_id}) || ($ncbi_transcript->{'data'}{'end'} - $ncbi_transcript->{'data'}{'start'}) > ($matched{$fixed_id}->{'data'}{'end'} - $matched{$fixed_id}->{'data'}{'start'})) {
							$matched{$fixed_id} = $ensembl_transcript;
						}
					}
				}
			}
		}
		$matched{$fixed_id}->addData({'fixed_id' => $fixed_id}) if (exists($matched{$fixed_id}));
	}
}

sub move_to_unbranded {
	my $root = shift;
	
	my $annotation_node = $root->findNode('updatable_annotation');
	my $annotation_set_nodes = $annotation_node->findNodeArray('annotation_set');
	my $annotation_set_ensembl;
	my $annotation_set_ncbi;
	my $annotation_set_lrg;
	my $current;
	
	foreach my $asn (@{$annotation_set_nodes}) {
		my $src = $asn->findNode('source');
		if ($src && $src->findNode('name')->content() eq 'Ensembl') {
			$annotation_set_ensembl = $asn;
		}
		elsif ($src && $src->findNode('name')->content() =~ m/NCBI/i) {
			$annotation_set_ncbi = $asn;
		}
		elsif ($src && $src->findNode('name')->content() eq 'LRG') {
			$annotation_set_lrg = $asn;
		}
	}

	my $old_annotation_set = LRG::Node::newFromNode($annotation_set_lrg);
	
# Add an LRG branded annotation_set if none is available
	if (!$annotation_set_lrg) {
		$annotation_set_lrg = $annotation_node->addNode('annotation_set');
		$current = $annotation_set_lrg->addNode('source');
		$current->addNode('name')->content('LRG');
		$current->addNode('url')->content('http://www.lrg-sequence.org/');
		$current = $current->addNode('contact');
		$current->addNode('name')->content('Locus Reference Genomic');
		$current->addNode('url')->content('http://www.lrg-sequence.org/page.php?page=contact');
		$current->addNode('email')->content('feedback@lrg-sequence.org');
		$annotation_set_lrg->addNode('modification_date')->content($root->date);
	}

# Attempt to get mapping data from annotation sets
	
	my $mapping_ncbi = $annotation_set_ncbi->findNodeArray('mapping') unless(!$annotation_set_ncbi);
	my $mapping_ensembl = $annotation_set_ensembl->findNodeArray('mapping') unless(!$annotation_set_ensembl);
	my $mapping_lrg = $annotation_set_lrg->findNodeArray('mapping') unless(!$annotation_set_lrg);

	my @mapping_all;
	foreach my $mapping (($mapping_ncbi,$mapping_ensembl,$mapping_lrg)) {
		if (defined($mapping)) {
			@mapping_all = (@mapping_all,@{$mapping});
		}
	}
# Remove all old mappings
	foreach my $mapping (@mapping_all) {
# Replace NCBI's 'NCBI37' assembly label with 'GRCh37'
		$mapping->{'data'}{'assembly'} =~ s/^NCBI37$/GRCh37/i;
		$mapping->remove();
	}
	
# Check the mappings for the different annotations.	
	if (scalar(@mapping_all) > 0) {
# Loop over the mappings and compare the matching assemblies
		for (my $i=0; $i<(scalar(@mapping_all)-1); $i++) {
			for (my $j=($i+1); $j<scalar(@mapping_all); $j++) {
				if ($mapping_all[$i]->{'data'}{'assembly'} eq $mapping_all[$j]->{'data'}{'assembly'}) {
					$mapping_all[$i] = compare_mapping($mapping_all[$i],$mapping_all[$j]);
					splice(@mapping_all,$j,1);
					$j--;
				}
			}
		}
# Insert the mappings in the unbranded section
		foreach my $mapping (@mapping_all) {
			$annotation_set_lrg->addExisting(LRG::Node::newFromNode($mapping));
		}
	}

# Attempt to get HGNC identifier from annotation sets
	my $hgnc_ncbi = $annotation_set_ncbi->findNode('lrg_gene_name') unless(!$annotation_set_ncbi);
	my $hgnc_ensembl = $annotation_set_ensembl->findNode('lrg_gene_name') unless(!$annotation_set_ensembl);
	my $hgnc_lrg = $annotation_set_lrg->findNode('lrg_gene_name') unless(!$annotation_set_lrg);
	
	my @hgnc_all = ($hgnc_ncbi,$hgnc_ensembl,$hgnc_lrg);
	
# Remove all old HGNC entries
	for (my $i=0; $i<scalar(@hgnc_all); $i++) {
		if (defined($hgnc_all[$i])) {
			$hgnc_all[$i]->remove();
		}
		else {
			splice(@hgnc_all,$i,1);
			$i--;
		}
	}
	
# Check to see if annotations agree. Die if they do not. Otherwise, add the mapping
	if (scalar(@hgnc_all) > 0) {
		while (scalar(@hgnc_all) > 1) {
			die("HGNC symbols for " . $hgnc_all[0]->{'parent'}->findNode('source/name')->content() . " and " . $hgnc_all[1]->{'parent'}->findNode('source/name')->content() . " annotations do not match!") unless($hgnc_all[0]->content() eq $hgnc_all[1]->content());
			my %data = (%{$hgnc_all[0]->{'data'}},%{$hgnc_all[1]->{'data'}});
			$hgnc_all[0]->{'data'} = \%data;
			splice(@hgnc_all,1,1);
		}
		$annotation_set_lrg->addExisting(LRG::Node::newFromNode($hgnc_all[0]));
	}
	
# Check if the annotations were changed, in that case, update the modification date
	if (!$old_annotation_set || !$annotation_set_lrg->identical($old_annotation_set)) {
# update the modification date
		$annotation_set_lrg->findOrAdd('modification_date')->content($root->date);
	}
# Finally, move the LRG annotation set to be the first among the updatable annotations
	$annotation_set_lrg->moveTo(0);
}

# Compare mapping data between mappings and mapping_spans. If not equal, the script will die with an error message
sub compare_mapping {
	my $mapping_1 = shift;
	my $mapping_2 = shift;
	
# Fields that are required to match in mapping and mapping_span tag, respectively
	my @map_fields = ('chr_name','chr_start','chr_end');
	my @span_fields = ('lrg_start','lrg_end','start','end');
	
# Check that each required field in the mapping tag is identical. Die if not
	foreach my $field (@map_fields) {
		if (!defined($mapping_1->{'data'}{$field}) || !defined($mapping_2->{'data'}{$field}) || $mapping_1->{'data'}{$field} ne $mapping_2->{'data'}{$field}) {
			die("There is a difference in the mappings for " . $mapping_1->{'data'}{'assembly'} . " between " . $mapping_1->{'parent'}->findNode('source/name')->content() . " and " . $mapping_2->{'parent'}->findNode('source/name')->content() . " " . $field . " attribute differs! (" . $mapping_1->{'data'}{$field} . " vs " . $mapping_2->{'data'}{$field} . ")");
		}
	}

# Merge the data from the mappings in order not to loose any information
	my %data_merged = (%{$mapping_1->{'data'}},%{$mapping_2->{'data'}});
	$mapping_1->{'data'} = \%data_merged;
	
# Check that the mapping contains equal number of spans. Die if not
	my @spans_1 = @{$mapping_1->findNodeArray('mapping_span')};
	my @spans_2 = @{$mapping_2->findNodeArray('mapping_span')};
	if (scalar(@spans_1) != scalar(@spans_2)) {
		die("The number of mapping spans is different for mappings in " . $mapping_1->{'parent'}->findNode('source/name')->content() . " and " . $mapping_2->{'parent'}->findNode('source/name')->content() . "! (" . scalar(@spans_1) . " vs " . scalar(@spans_2));
	}
# For each span, find the matching span from the other source and compare the required fields. The matching span is found if any of the start or end coordinates match (lrg or chromosome)
	foreach my $span_1 (@spans_1) {
		my $found = 0;
		foreach my $span_2 (@spans_2) {
			my $sum = 0;
			foreach my $field (@span_fields) {
				$sum += ($span_1->{'data'}{$field} == $span_2->{'data'}{$field});
			}
			if ($sum > 0) {
				$found = 1;
				die("There is a difference in the mapping span coordinates for " . $mapping_1->{'parent'}->findNode('source/name')->content() . " and " . $mapping_2->{'parent'}->findNode('source/name')->content() . "!") if ($sum < scalar(@span_fields) || $span_1->{'data'}{'strand'} != $span_2->{'data'}{'strand'});
				my %span_merged = (%{$span_1->{'data'}},%{$span_2->{'data'}});
				$span_1->{'data'} = \%span_merged;
# For each diff within a span, require that the same diff is present in the other source
				my @diffs_1 = $span_1->findNodeArray('diff');
				my @diffs_2 = $span_2->findNodeArray('diff');
				die("The number of diffs within one of the spans are not the same between " . $mapping_1->{'parent'}->findNode('source/name')->content() . " and " . $mapping_2->{'parent'}->findNode('source/name')->content() . "!") unless (scalar(@diffs_1) == scalar(@diffs_2));
				foreach my $diff_1 (@diffs_1) {
					my $matched = 0;
					foreach my $diff_2 (@diffs_2) {
						my $n = scalar(keys %{$diff_1->{'data'}});
						my $i=0;
						while (my ($key,$val) = each(%{$diff_1->{'data'}})) {
							$i += ($diff_2->{'data'}{$key} eq $val);
						}
						if ($i == $n) {
							$matched = 1;
							last;
						}
						
					}
					die("The diffs within one of the spans don't match between " . $mapping_1->{'parent'}->findNode('source/name')->content() . " and " . $mapping_2->{'parent'}->findNode('source/name')->content() . "!") unless ($matched);
				}
			}
		}
		die("The spans does not match between " . $mapping_1->{'parent'}->findNode('source/name')->content() . " and " . $mapping_2->{'parent'}->findNode('source/name')->content() . "!") unless($found);
	}

# Return the first mapping node where the union of all information has been stored
	return $mapping_1;
}

# Extract cDNA sequence from a reference sequence using specified exon coordinates
sub make_cdna {
	my $reference = shift;
	my $exons = shift;
	
	my $cdna;
	foreach my $exon (@{$exons}) {
		my $start = $exon->{'start'};
		my $end = $exon->{'end'};
		$cdna .= substr($reference,$start,$end-$start+1);
	}
	return $cdna;
}

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