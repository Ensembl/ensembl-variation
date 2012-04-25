use strict;
use warnings;

package AncestralAlleleCaller;
use Bio::DB::Fasta;

sub new {
	my $class = shift;
	my $fasta_files_dir = shift;
	my $db = Bio::DB::Fasta->new($fasta_files_dir, -reindex => 1);
	my @sequence_ids = $db->get_all_ids;
	my %sequence_id_2_chr_number;

	foreach my $sequence_id (@sequence_ids) {
		my @split = split(/:/, $sequence_id);
		$sequence_id_2_chr_number{$split[2]} = $sequence_id;
		print $sequence_id, "\t", $split[2], "\n";
	}
	my $self = bless {
		'db'                       => \$db,
		'sequence_id_2_chr_number' => \%sequence_id_2_chr_number,
	}, $class;
	return $self;
}

sub get_ancestral_allele {
	my $self        = shift;
	my $chr         = shift;
	my $sequence_id = $self->_get_sequence_id($chr);
	my $position_from = shift;
    my $position_to   = shift;
	return ${$self->{'db'}}->seq("$sequence_id:$position_from,$position_to"); 
}

sub _get_sequence_id {
	my $self = shift;	
	my $chr  = shift;
	my $hash = $self->{'sequence_id_2_chr_number'};
	return ${$self->{'sequence_id_2_chr_number'}}{$chr};
}
1;

