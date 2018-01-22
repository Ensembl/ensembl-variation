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

Bio::EnsEMBL::Variation::Pipeline::VariantQC::FlipPopulationGenotype

=head1 DESCRIPTION

Compliment alleles in population_genotype table where required
Runs after allele/ variation feature flipping & QC

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

    my $flip_ext_sth = $var_dba->dbc->prepare(qq[ select variation_id from variation_working where flipped = 1]);
    $flip_ext_sth->execute() ||die "Failed to extract flip list\n";;
    while(my $l = $flip_ext_sth->fetchrow_arrayref()){
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


    my $data_ext_sth = $var_dba->dbc->prepare(qq[ select population_genotype_id, variation_id, subsnp_id, allele_1, allele_2, frequency, population_id, count 
                                                  from population_genotype
                                                  where variation_id between ? and ?
                                                ]);

    my $data_ins_sth = $var_dba->dbc->prepare(qq[ insert into MTMP_population_genotype_working
                                                 (population_genotype_id, variation_id, subsnp_id, allele_1, allele_2, frequency, population_id, count) 
                                                 values(?,?,?,?,?,?,?,?)
                                                ]);

    $var_dba->dbc->do(qq{ ALTER TABLE MTMP_population_genotype_working DISABLE KEYS});
   
    my $start = 1;
    my $batch = 10000;

    while (  $start < $max_var_id){

        my $end = $start  + $batch;
        $data_ext_sth->execute($start, $end)|| die "ERROR extracting population_genotype data\n";
  
        while (my $l = $data_ext_sth->fetchrow_arrayref() ){
            
            if($flip->{$l->[1]} == 1){
                ## update genotypes for flipped variants                
                defined $QUICK_COMP{$l->[3]} ? $l->[3] = $QUICK_COMP{$l->[3]} : reverse_comp(\$l->[3]);
                defined $QUICK_COMP{$l->[4]} ? $l->[4] = $QUICK_COMP{$l->[4]} : reverse_comp(\$l->[4]);
            }
            ## write to Mart table
            $data_ins_sth->execute($l->[0], $l->[1], $l->[2], $l->[3], $l->[4], $l->[5], $l->[6], $l->[7]) ||die "Failed to enter pop gen data\n";
        }
        $start = $end +1;
    }
    $var_dba->dbc->do(qq{ ALTER TABLE MTMP_population_genotype_working ENABLE KEYS});

}



1;
