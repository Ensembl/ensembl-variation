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
    
    my $get_stable_id_sth = $var_dba->dbc->prepare(qq{
        SELECT  stable_id
        FROM    translation_mapping
        WHERE   md5 = ?
    });

    $get_stable_id_sth->execute($md5);

    my ($stable_id) = $get_stable_id_sth->fetchrow_array;

    return $stable_id;
}

1;
