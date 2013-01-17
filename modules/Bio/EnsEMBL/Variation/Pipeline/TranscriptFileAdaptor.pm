=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/legal/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut
package Bio::EnsEMBL::Variation::Pipeline::TranscriptFileAdaptor;

use strict;

use Digest::MD5 qw(md5_hex);

sub new {
    my $class = shift;

    my %args = @_;

    my $self = bless {}, $class;

    if ($args{fasta_file}) {
        $self->{fasta_file} = $args{fasta_file};
    }
    
    if ($args{transcripts}) {
        $self->_dump_translations($args{transcripts});
    }
    
    return $self;
}

sub get_translation_seq {
    my ($self, $translation_md5) = @_;

    my $fasta = $self->get_translation_fasta($translation_md5);

    $fasta =~ s/>.*\n//m;
    $fasta =~ s/\s//mg;

    return $fasta;
}

sub get_translation_fasta {
    my ($self, $translation_md5) = @_;
    
    my $file = $self->{fasta_file};
    
    my $fasta = `samtools faidx $file $translation_md5`;

    return $fasta;
}

sub get_all_translation_md5s {
    my $self = shift;

    my $fasta = $self->{fasta_file};

    my @ids = map {/>(.+)\n/; $1} `grep '>' $fasta`;

    return \@ids;
}

sub _dump_translations {

    my ($self, $transcripts) = @_;

    # dump the translations out to the FASTA file

    my $fasta = $self->{fasta_file};

    open my $FASTA, ">$fasta" or die "Failed to open $fasta for writing";

    # get rid of any existing index file

    if (-e "$fasta.fai") {
        unlink "$fasta.fai" or die "Failed to delete fasta index file";
    }

    my %seen_md5;

    for my $transcript (@$transcripts) {

        my $tl = $transcript->translation;
    
        next unless $tl;

        my $protein = $tl->seq;

        my $md5 = md5_hex($protein);

        next if $seen_md5{$md5}++;

        $protein =~ s/(.{80})/$1\n/g;

        # get rid of any trailing newline
        chomp $protein;

        print $FASTA ">$md5\n$protein\n";
    }

    close $FASTA;

    # index the file

    `samtools faidx $fasta`;
}
 

1;
