=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

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
