
=head1 LICENSE

 Copyright (c) 1999-2012 The European Bioinformatics Institute and
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

Bio::EnsEMBL::Variation::Pipeline::VariantQC::FlipPopulationGenotype

=head1 DESCRIPTION

Compliment alleles in population_genotype table where required
Runs alongside allele/ variation feature flipping & QC

=cut

package Bio::EnsEMBL::Variation::Pipeline::VariantQC::FlipPopulationGenotype;


use strict;
use warnings;

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);


sub run{

    my $self = shift;
   
    my $var_dba = $self->get_species_adaptor('variation');

    ## check there is something to do
    my $count_ext_sth  = $var_dba->dbc->prepare(qq[ select max(variation_id) from population_genotype]);
    $count_ext_sth->execute()||die "Failed to check population genotype \n";
    my $max = $count_ext_sth->fetchall_arrayref();

    return unless defined  $max->[0]->[0];

    ## get list of variations to flip
    my $flip = get_flip_list($var_dba);
   
    ## create new table, flipping where needed
    flip_genotypes( $var_dba, $max->[0]->[0], $flip );
}


sub get_flip_list{

    my $var_dba = shift;
    my %flip;

    my $flip_ext_sth = $var_dba->dbc->prepare(qq[ select  variation_id from variation_to_reverse_working ]);
    $flip_ext_sth->execute() ||die "Failed to extract flip list\n";;
    my $var_list = $flip_ext_sth->fetchall_arrayref();
    foreach my $l(@{$var_list}){
        $flip{$l->[0]} = 1;
    }

    return (\%flip);
}


sub flip_genotypes{

    my $var_dba    = shift;
    my $max_var_id = shift;
    my $flip       = shift;

    my %QUICK_COMP = ( "A" => "T",
                       "T" => "A",
                       "C" => "G",
                       "G" => "C"
    ); 


    my $data_ext_sth = $var_dba->dbc->prepare(qq[ select population_genotype_id, variation_id, subsnp_id, allele_1, allele_2, frequency, sample_id, count 
                                                  from population_genotype
                                                  where variation_id between ? and ?
                                                ]);

    my $data_ins_sth = $var_dba->dbc->prepare(qq[ insert into MTMP_population_genotype_working
                                                 (population_genotype_id, variation_id, subsnp_id, allele_1, allele_2, frequency, sample_id, count) 
                                                 values(?,?,?,?,?,?,?,?)
                                                ]);

   
    my $start = 1;
    my $batch = 10000;

    while (  $start < $max_var_id){

        my $end = $start  + $batch;
        $data_ext_sth->execute($start, $end)|| die "ERROR extracting population_genotype data\n";
  
        while (my $l = $data_ext_sth->fetchrow_arrayref() ){
            
            if($flip->{$l->[1]}){
                ## update genotypes for flipped variants                
                defined $QUICK_COMP{$l->[3]} ? $l->[3] = $QUICK_COMP{$l->[3]} : reverse_comp(\$l->[3]);
                defined $QUICK_COMP{$l->[4]} ? $l->[4] = $QUICK_COMP{$l->[4]} : reverse_comp(\$l->[4]);
            }
            ## write to Mart table
            $data_ins_sth->execute($l->[0], $l->[1], $l->[2], $l->[3], $l->[4], $l->[5], $l->[6], $l->[7]) ||die "Failed to enter pop gen data\n";
        }
        $start = $end +1;
    }

}



1;
