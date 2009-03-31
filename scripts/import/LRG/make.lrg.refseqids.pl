#!/software/bin/perl

use strict;
use lib './';
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
my $in_file;
my $template_file = 'template.xml';
my $registry_file = "ensembl.registry";
my $input_dir;
my $output_dir;
my $target_dir;

# get options from command line
GetOptions(
    'input_file=s' => \$in_file,
    'output_file_stem=s' => \$out_file_stem,
    'template_file=s' => \$template_file,
    'registry_file=s' => \$registry_file,
    'input_dir=s' => \$input_dir,
	'output_dir=s' => \$output_dir,
	'target_dir=s' => \$target_dir,
);

$input_dir ||= "tempin";
$output_dir ||= "tempout";
$target_dir ||= "/lustre/work1/ensembl/yuan/SARA/human/ref_seq_hash";

our $mapping_num = 1;

# get registry and a gene adaptor
# Bio::EnsEMBL::Registry->load_all( $registry_file );
Bio::EnsEMBL::Registry->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');

my $species = 'human';
my $cdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $dbCore = $cdb;
my $sa = $dbCore->get_SliceAdaptor();

my $geneAd = Bio::EnsEMBL::Registry->get_adaptor('Homo sapiens' ,'Core', 'Gene');
my $transAd = Bio::EnsEMBL::Registry->get_adaptor('Homo sapiens', 'Core', 'Transcript');
my $protAd = Bio::EnsEMBL::Registry->get_adaptor('Homo sapiens', 'Core', 'Protein');

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
	$current->addEmptyNode(
		'align',
		{
			'chromosome' => $mapping->{'t_id'},
			'start' => $mapping->{'t_start'},
			'end' => $mapping->{'t_end'},
			'lrg_start' => $mapping->{'q_start'},
			'lrg_end' => $mapping->{'q_end'},
			'strand' => $mapping->{'q_strand'},
		}
	);
	
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
# 	die;
	
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
			
			print "Mapping from $current_start to ", $current_start + $split_size, "\n";
			push @maps, mapping(substr($sequence, $current_start, $split_size));
			
			$current_start += $split_size;
		}
		
		for my $i(0..$#maps) {
			my $map = $maps[$i];
			print $i, ": ", $map->{'t_start'}, " ", $map->{'t_end'}, "\t", $map->{'q_start'}, " ", $map->{'q_end'}, "\n";
		}
		
		# join the maps
		my $main_map = shift @maps;
		
		my $prev_q_end = $main_map->{'q_end'};
		
		foreach my $map(@maps) {
			$main_map->{'t_end'} = $map->{'t_end'};
			
			$main_map->{'q_end'} = $map->{'q_end'} + $prev_q_end;
			
			$prev_q_end = $main_map->{'q_end'};
		}
		
		return $main_map;
	}
	
	else {
	
		my $name = "temp".$$.$mapping_num++;
	# 	my $name = "temp_will";
		my $input_file_name = $name.'.fa';
		
		open OUT, ">$input_dir/$name\.fa" or die ;
		print OUT '>', $name, "\n", $sequence;
		close OUT;
		
		my $queue = "-q normal -R'select[mem>4000] rusage[mem=4000]' -M4000000";
		
		my $output_file_name = "ssaha2_output\_$input_file_name";
		my $input_file = "$input_dir/$input_file_name";
		my $output_file ="$output_dir/$output_file_name";
		
		my ($subject, %rec_seq, %rec_find, %input_length, %done);
		
		$subject = "$target_dir/ref";
		
		my %rec_seq;
		$rec_seq{$name} = $sequence;
			
		foreach my $name (keys %rec_seq) {#get query seq object
			#print "name is $name and seq is $rec_seq{$name}\n";
			my $seqobj = Bio::PrimarySeq->new(-id => $name, -seq => $rec_seq{$name});
			$rec_seq{$name} = $seqobj;
		}
		
		bsub_ssaha_job($queue,$input_file,$output_file,$subject);
		
		my $call = "bsub -P ensembl-variation -K -w 'done($input_file\_ssaha_job)' -J waiting_process sleep 1"; #waits until all ssaha jobs have finished to continue
		system($call);
		
		my $mapping = parse_ssaha2_out($output_file);
		
		# if we have more than one hit
		if(scalar @$mapping > 1 && $mapping->[0]->{'score'} == $mapping->[1]->{'score'}) {
			my $count = 0;
		
			foreach my $hit(@$mapping) {
				$count++ if $mapping->[0]->{'score'} == $hit->{'score'};
			}
			
			warn "Mapped $count times\n";
		}
		
		else {
			return $mapping->[0];
		}
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
			#next unless $1 eq "rs16822290"
# 			print;
			
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
			
			if($h->{'q_start'} > $h->{'q_end'}) {
				($h->{'q_start'}, $h->{'q_end'}) = ($h->{'q_end'}, $h->{'q_start'});
				($h->{'t_start'}, $h->{'t_end'}) = ($h->{'t_end'}, $h->{'t_start'});
			}
			
			push @rec_find, $h if $h->{'q_id'};
		}
	}
	
	close SSAHA;
	
	@rec_find = sort {$b->{'score'} <=> $a->{'score'}} @rec_find;

# 	print @rec_find, "\n";
	
	return \@rec_find;
}
