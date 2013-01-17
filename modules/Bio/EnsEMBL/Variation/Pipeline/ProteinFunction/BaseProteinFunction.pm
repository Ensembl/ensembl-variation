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
package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::BaseProteinFunction;

use strict;
use warnings;

use Bio::SimpleAlign;
use Bio::LocatableSeq;

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');

sub get_protein_sequence {
    my ($self, $md5) = @_;

    my $tfa = $self->get_transcript_file_adaptor;

    my $fasta = $tfa->get_translation_fasta($md5);

    unless (length($fasta) > 34) {
        # sleep in case it's some race condition and then try again
        sleep(rand(5));
        $fasta = $tfa->get_translation_fasta($md5);

        unless (length($fasta) > 34) {
            die "$md5 looks wierdly short!";
        }
    }

    # strip out fasta header etc.

    $fasta =~ s/>.*\n//m;
    $fasta =~ s/\s//mg;
    
    die "No peptide for $md5?" unless length($fasta) > 0;

    return $fasta;
}

sub get_stable_id_for_md5 {
    my ($self, $md5) = @_;

    my $var_dba = $self->get_species_adaptor('variation');
    
    my $get_stable_id_sth = $var_dba->prepare(qq{
        SELECT  stable_id
        FROM    translation_mapping
        WHERE   md5 = ?
    });

    $get_stable_id_sth->execute($md5);

    my ($stable_id) = $get_stable_id_sth->fetchrow_array;

    return $stable_id;
}

1;
