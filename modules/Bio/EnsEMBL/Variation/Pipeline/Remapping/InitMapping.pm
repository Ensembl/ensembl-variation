package Bio::EnsEMBL::Variation::Pipeline::Remapping::InitMapping;

use strict;
use warnings;

use FileHandle;
use Bio::DB::Fasta;
use Bio::EnsEMBL::Registry;

use base ('Bio::EnsEMBL::Hive::Process');


# use variation database
sub fetch_input {
	my $self = shift;
	my $species = $self->param('species');

	if ($self->param('generate_fasta_files')) {
		my $fasta_files_dir             = $self->param('fasta_files_dir');
		my $old_assembly_fasta_file_dir = $self->param('old_assembly_fasta_file_dir');
		my $registry_file               = $self->param('registry_file');		

		my $db = Bio::DB::Fasta->new($old_assembly_fasta_file_dir);
		$self->param('fasta_db', $db);
		my $registry = 'Bio::EnsEMBL::Registry';
		$registry->load_all($registry_file);
		my $sa  = $registry->get_adaptor($species, 'core', 'slice');
		my $vfa = $registry->get_adaptor($species, 'variation', 'variationfeature');
		$vfa->db->include_failed_variations(1);
	
		$self->param('sa', $sa);
		$self->param('vfa', $vfa);
	}
	1;	
}

sub run {
	my $self = shift;
	# write report, how many variation feature
	# collect variants with multiple mappings;
	my $fasta_files_dir = $self->param('fasta_files_dir');
	my $file_count;	

	if ($self->param('generate_fasta_files')) {
		$file_count       = 1;
		my $count_entries = 0;
		my $fh_query_sequences = FileHandle->new("$fasta_files_dir/$file_count.fa", 'w');
		my $slice_adaptor             = $self->param('sa');
		my $variation_feature_adaptor = $self->param('vfa');

		# get toplevel seq_region_names
		my $slices = $slice_adaptor->fetch_all('toplevel', undef, 0, 1);

		if ($self->param('debug')) {	
			my $debug_slice = $slice_adaptor->fetch_by_region('chromosome', 1);	
			@$slices = ();
			push @$slices, $debug_slice;	
		}
		
		my $vfs_with_multi_map;

		my $count_vfs = 0;
		foreach my $slice (@$slices) {
			my $slice_end       = $slice->end;
			my $seq_region_name = $slice->seq_region_name();	
			my $vfi             = $variation_feature_adaptor->fetch_Iterator_by_Slice($slice);	
			while (my $vf = $vfi->next)	{
				# build only one query_seq for variants with multiple mappings
				if ($vf->map_weight > 1) {
					next if ($vfs_with_multi_map->{$vf->variation_name});
				} else {
					$vfs_with_multi_map->{$vf->variation_name} = 1;
				}
				$count_vfs++;
				my $entry = $self->get_query_sequence($vf, $slice_end);
				if ($count_entries >= 100_000) {
					$fh_query_sequences->close();
					$file_count++;
					$count_entries = 0;
					$fh_query_sequences = FileHandle->new("$fasta_files_dir/$file_count.fa", 'w');
				}
				$count_entries++;
				print $fh_query_sequences $entry;	
			}
		}
		$fh_query_sequences->close();	
	} else {
		$file_count = 0;
		opendir(DIR, $fasta_files_dir) or die $!;
		while (my $file = readdir(DIR)) {
			if ($file =~ m/\.fa$/) {
				$file_count++;
			}
		}
		closedir(DIR);
	}
	$self->param('file_count', $file_count);	
	1;
}


sub write_output {
	my $self = shift;
	# initialise mapping jobs 
	my $file_count      = $self->param('file_count');
	my $fasta_files_dir = $self->param('fasta_files_dir');
	my $bam_files_dir   = $self->param('bam_files_dir');
	my @jobs;
	my $i = 1;
	while ($i <= $file_count) {
		push @jobs, {
            'file_number'   => $i,
			'bam_files_dir' => $bam_files_dir,
			'fasta_file'    => "$fasta_files_dir/$i.fa",
			'sam_file'      => "$bam_files_dir/$i.sam",
			'bam_file'      => "$bam_files_dir/$i.bam",
			'err_file'      => "$bam_files_dir/$i.err",
			'out_file'      => "$bam_files_dir/$i.out",
		};	
		$i++;
	}
	$self->dataflow_output_id(\@jobs, 2);
	1;
}

sub get_query_sequence {
	my $self = shift;
	my $vf   = shift;
	my $slice_end = shift;
	my $fasta_db  = $self->param('fasta_db');
	my $flank_seq_length = $self->param('flank_seq_length');

	
	my $sequence_name   = $vf->seq_region_name;
	my $variation_name  = $vf->variation_name;
	my $vf_id           = $vf->dbID;
	my $variation_start = $vf->seq_region_start;
	my $variation_end   = $vf->seq_region_end;
	# ...-|B|--...--|VS|-...-|VE|--...--|A|-...
	my $before_variation_start = $variation_start - $flank_seq_length; # |B| 
	my $after_variation_end    = $variation_end + $flank_seq_length; # |A|
	
	if ($before_variation_start <= 0) {
		$before_variation_start = 1;
		$after_variation_end = $variation_end + $flank_seq_length + ($flank_seq_length - $variation_start);
	}

	if ($after_variation_end > $slice_end) {
		$after_variation_end = $slice_end;
		my $diff = $slice_end - $variation_end;
		$before_variation_start = $variation_start - (2 * $flank_seq_length) + $diff;
	}
 
	my $query_sequence = $fasta_db->seq("$sequence_name:$before_variation_start,$after_variation_end");

	my $length_flank_seq_before = $variation_start - $before_variation_start;
	my $length_variation = 0; # only for insertions length is 0
	if ($variation_end >= $variation_start) {
		$length_variation = $variation_end - $variation_start + 1;
	}
	my $length_flank_seq_after = $after_variation_end - $variation_end;

	my $id = ">$vf_id-$length_flank_seq_before-$length_variation-$length_flank_seq_after-$variation_name";

	return "$id\n$query_sequence\n";	

}

1;
