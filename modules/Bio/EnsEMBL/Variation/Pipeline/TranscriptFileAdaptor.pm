package Bio::EnsEMBL::Variation::Pipeline::TranscriptFileAdaptor;

use strict;

use Digest::MD5 qw(md5_hex);

sub new {
    my $class = shift;
    my $args  = {@_};
    return bless $args,  $class;
}

sub get_protein_seq {
    my ($self, $transcript_stable_id) = @_;
    
    my $fasta = $self->{fasta_file};
    
    my $seq = `samtools faidx $fasta $transcript_stable_id`;

    $seq =~ s/>.*\n//m;
    $seq =~ s/\s//mg;

    return $seq;
}

sub get_protein_md5 {
    my ($self, $transcript_stable_id) = @_;

    return md5_hex($self->get_protein_seq($transcript_stable_id));
}

sub get_all_transcript_stable_ids {
    my $self = shift;

    my $fasta = $self->{fasta_file};

    my @ids = map {/>(.+)\n/; $1} `grep '>' $fasta`;

    return \@ids;
}

1;
