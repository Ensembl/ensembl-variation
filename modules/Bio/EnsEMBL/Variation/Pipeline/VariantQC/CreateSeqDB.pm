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



=head1 NAME 

Bio::EnsEMBL::Variation::Pipeline::VariantQC::CreateSeqDB

=head1 DESCRIPTION

Export and index fasta file of reference genome to avoid using database 
for species with large numbers of variants

=cut

package Bio::EnsEMBL::Variation::Pipeline::VariantQC::CreateSeqDB;


use strict;
use warnings;
use Bio::DB::Fasta;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);



sub run {

    my $self = shift;

    my $dir = $self->required_param('pipeline_dir');

    my $fasta_seq = $dir ."/genomic.fa";

    open my $fh, ">$fasta_seq" or die "can't open fasta file $! \n";

    ## get list of required seq ids
    ## not using fetch_all("toplevel"..) anymore as MHC haplotypes with variants are not top level
    my $var_dba = $self->get_species_adaptor('variation');
    my $id_ext_sth = $var_dba->dbc->prepare(qq [ select distinct seq_region_id from variation_feature ]);
    $id_ext_sth->execute()||die;
    my $seq_ids = $id_ext_sth->fetchall_arrayref();

    my $core_dba = $self->get_species_adaptor('core');
    my $slice_adaptor = $core_dba->get_SliceAdaptor();

     
    foreach my $id(@{$seq_ids}){
	next if $id->[0] ==0;

	my $slice = $slice_adaptor->fetch_by_seq_region_id($id->[0]);

	unless (defined $slice){
	    $self->warning('No slice found for seq_region_id '.$id->[0] );  
	    next;
	}
	print $fh ">".$slice->seq_region_name()."\n";
	my $seq = $slice->seq;
	## format & print seq
	fasta($seq, $fh);  

    }
    close $fh or die "failed to close sequence file: $!\n"; 

    ## create index file 
    my $fasta_db = Bio::DB::Fasta->new($fasta_seq);
}



## print sequence string in lines of 80 bases

sub fasta {

    my $seq = shift;
    my $fh  = shift;

    my $len = length($seq);
    my $start = 0;
    while($start < $len) {
        print $fh substr($seq, $start, 80), "\n";
        $start += 80;
    }
}




1;
