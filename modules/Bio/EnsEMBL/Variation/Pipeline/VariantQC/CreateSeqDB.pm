
=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.



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
