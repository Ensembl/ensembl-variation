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

=head1 NAME

Bio::EnsEMBL::Variation::Utils::ComparaUtils

=head1 DESCRIPTION

This module exports two subroutines, dump_alignment_for_polyphen and
dump_alignment_for_sift that write Compara alignments to files in the 
formats expected by PolyPhen and SIFT. This allows you to use
Compara alignments in place of both these tools' alignment pipelines

=cut

package Bio::EnsEMBL::Variation::Utils::ComparaUtils;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Registry;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::LocatableSeq;

use base qw(Exporter);

our @EXPORT_OK = qw(dump_alignment_for_polyphen dump_alignment_for_sift);

my $MAX_PSIC_SEQS   = 8190;
my $MAX_PSIC_SEQLEN = 409650;

sub _ungap_alignment {
    
    # turn a gapped alignment into an ungapped alignment
    # with respect to the given query_id,
    #
    # e.g. if our query sequence is no. 1 below, we want to 
    # turn this:
    #
    # 1: -ABC--DEFG-H---
    # 2: ---IJKLMN--OPQR
    # 3: ST------UVWXYZ-
    #
    # into:
    #
    # 1: ABCDEFGH
    # 2: --ILMN-O
    # 3: T----UVX

    my ($gapped_alignment, $query_id, $include_query) = @_;

    my @seqs = $gapped_alignment->each_seq;

    # find the Seq object corresponding to the query id 

    my $query_seq;

    for my $seq (@seqs) {
        $query_seq = $seq;
        last if $seq->display_id eq $query_id;
    }

    throw("Could not find query sequence '$query_id' in the alignment")
        unless $query_seq->display_id eq $query_id;

    my $qseq = $query_seq->seq;

    #print $qseq, " (".length($qseq).")\n";

    # first we split the query sequence into its component parts and gaps,
    # so -ABC--DEFG-H--- will become ('ABC','DEFG','H') and ('-','--','-','---')
    
    # we grep because we don't want any empty strings that otherwise seem to 
    # make it through

    my @parts = grep {$_} split /-+/, $qseq;
    my @gaps = grep {$_} split /[^-]+/, $qseq;

    my @chunks;

    # then we compute the coordinates of each sequence chunk,
    # taking into account that we might start with a gap

    my $offset = substr($qseq,0,1) eq '-' ? length(shift @gaps) : 0;

    # we build up a list of [start, length] pairs for each chunk
    # of sequence, incrementing our offset by the length of each 
    # chunk and the gap that follows it, for our example above this 
    # will result in ([1,3],[6,4],[11,1])

    for my $part (@parts) {

        my $l = length($part);

        push @chunks, [$offset, $l];

        #print "$part: $offset - $l\n";

        $offset += $l + length(shift @gaps || '');
    }

    # we then use this list of chunks to obtain the aligned portions
    # of all other members of the alignment and create a new alignment,
   
    my @seq_objs;

    for my $seq (@seqs) {
        
        my $old_seq = $seq->seq;
        my $new_seq;

        for my $chunk (@chunks) {
            $new_seq .= substr($old_seq, $chunk->[0], $chunk->[1]);
        }

        next if $new_seq =~ /^-*$/;

        #my $gap_count = ($new_seq =~ tr/-//);

        #next if $gap_count / length($new_seq) > 0.1;

        if ($seq->display_id eq $query_id) {

            $query_seq =  Bio::LocatableSeq->new(
                -SEQ    => $new_seq,
                -START  => 1,
                -END    => length($new_seq),
                -ID     => 'QUERY',
                -STRAND => 0
            );
            
            unshift @seq_objs, $query_seq if $include_query;
        }
        else {
            push @seq_objs,  Bio::LocatableSeq->new(
                -SEQ    => $new_seq,
                -START  => 1,
                -END    => length($new_seq),
                -ID     => $seq->display_id,
                -STRAND => 0
            )
        }
    }

    my $new_align = Bio::SimpleAlign->new;

    $new_align->add_seq($_) for @seq_objs;

    # sometimes we want the query sequence back as well as the new alignment

    return wantarray ? ($query_seq, $new_align) : $new_align;
}

sub _get_ungapped_alignment {

    # get an ungapped Bio::SimpleAlign for the given query translation

    my ($translation_stable_id, $include_query) = @_;

    my $compara_dba = Bio::EnsEMBL::Registry->get_DBAdaptor('multi', 'compara')
        or throw("Failed to get compara DBAdaptor");

    my $sma = $compara_dba->get_SeqMemberAdaptor
        or throw("Failed to get seq_member adaptor");

    my $fa = $compara_dba->get_FamilyAdaptor
        or throw("Failed to get family adaptor");

    my $seq_member = $sma->fetch_by_stable_id($translation_stable_id)
        or throw("Didn't find family seq_member for $translation_stable_id");

    my $fam = $fa->fetch_by_SeqMember($seq_member)
        or throw("Didn't find a family for $translation_stable_id");

    my $orig_align = $fam->get_SimpleAlign;

    $compara_dba->dbc->disconnect_if_idle;

    return _ungap_alignment(
        $orig_align, 
        $translation_stable_id, 
        $include_query
    );
}

sub _percent_id {
    my ($q, $a) = @_;

    my @q = split //, $q->seq;
    my @a = split //, $a->seq;

    my $num = scalar(@q);
        
    my $tot = 0.0;

    for (my $i = 0; $i < $num; $i++ ) {
        $tot++ if $q[$i] eq $a[$i];
    }

    return ($tot / $num);
}

=head2 dump_alignment_for_polyphen

  Arg[1]      : string $translation_stable_id - the stable of the Ensembl translation
                you want to run PolyPhen on
  Arg[2]      : string $file - the name of the file you want to write the alignment to
  Description : Fetches the Compara protein family containing the specified translation
                (if available), ungaps the alignment with respect to the translation, and 
                writes the alignment to the specified file in the format expected by PolyPhen
  Returntype  : none
  Exceptions  : throws if an alignment cannot be found, or if the file cannot be written 
  Status      : At Risk
  
=cut

sub dump_alignment_for_polyphen {

    my ($translation_stable_id, $file) = @_;

    # polyphen does not want the query included in the alignment

    my ($query_seq, $alignment) = _get_ungapped_alignment($translation_stable_id, 0);

    my @seqs = $alignment->each_seq;
    
    unless (scalar(@seqs)) {
        throw("No sequences in the alignment for $translation_stable_id");
    }

    if (length($seqs[0]->seq) > $MAX_PSIC_SEQLEN) {
        throw("$translation_stable_id sequence too long for PSIC");
    }

    # polyphen expects the alignment to be sorted descending by % id

    my $percent_id;

    for my $seq (@seqs) {
        $percent_id->{$seq->id} = _percent_id($query_seq, $seq);
    }

    my @sorted = sort { $percent_id->{$b->id} <=> $percent_id->{$a->id} } @seqs;

    # the format is similar to a clustalw .aln file, but it *must* have a
    # 70 character description beginning each line, and then the alignment string
    # on the rest of the line with no line breaks between sequences, PSIC
    # can also only deal with a fixed maximum sequence length and number
    # of sequences on the alignment, these are set in the constants at the
    # top of this file

    my $num_seqs = scalar(@seqs);

    open my $ALN, ">$file" or throw("Failed to open $file for writing");

    print $ALN "CLUSTAL $translation_stable_id (".($num_seqs > $MAX_PSIC_SEQS ? $MAX_PSIC_SEQS : $num_seqs).")\n\n";

    my $count = 0;

    for my $seq (@sorted) {
        if ($percent_id->{$seq->id} == 1) {
            #warn "Ignoring identical seq ".$seq->id."\n";
            #next;
        }

        my $id = $seq->id;
        my $extra = 70 - length($id);
        
        print $ALN $seq->id, ' ' x $extra, $seq->seq, "\n";                

        last if ++$count >= $MAX_PSIC_SEQS;
    }

    close $ALN;
}

=head2 dump_alignment_for_sift

  Arg[1]      : string $translation_stable_id - the stable id of the Ensembl translation
                you want to run SIFT on
  Arg[2]      : string $file - the name of the file you want to write the alignment to
  Description : Fetches the Compara protein family containing the specified translation
                (if available), ungaps the alignment with respect to the translation, and 
                writes the alignment to the specified file in the format expected by SIFT
  Returntype  : none
  Exceptions  : throws if an alignment cannot be found, or if the file cannot be written 
  Status      : At Risk
  
=cut

sub dump_alignment_for_sift {

    my ($translation_stable_id, $file) = @_;

    # SIFT is easy, we just fetch an ungapped alignment including the query 
    # and dump it in FASTA format

    my $aln = _get_ungapped_alignment($translation_stable_id, 1);

    open my $ALN, ">$file" or throw("Failed to open $file for writing");

    my $alignIO = Bio::AlignIO->newFh(
        -interleaved => 0,
        -fh          => $ALN,
        -format      => 'fasta',
        -idlength    => 20,
    );

    print $alignIO $aln;

    close $ALN;
}

1;

