
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

=cut

package Bio::EnsEMBL::Variation::Pipeline::VariantQC::FinishVariantQC;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

my $DEBUG = 0;

## This runs some basic checks on new tables to produce a report
## swap tables names if they look correct  [TO DO]

sub run {
   
    my $self = shift;

    my $dir = $self->required_param('pipeline_dir');

    open my $report, ">", "$dir/QC_report.txt" || die "Failed to open QC_report.txt : $!\n";
    print $report "\n\tChecking post-QC tables \n\n";

    my $var_dba  = $self->get_species_adaptor('variation');


    ## Have all rows been processed
    my ($allele_number, $all_ok)  = check_all_processed($var_dba, $report);

     
    ## What are the failure rates for alleles and variants
    get_failure_rates($var_dba, $report, $allele_number );


    ## what are failure reasons for alleles and variants    
    check_failure_rates($var_dba, $report);


    ## rename tables
    if( $all_ok ==2 ){
	print $report "\n\tOK to rename post-QC tables \n\n";
	# rename_tables($var_dba);
    }

}



 ## Have all rows been processed
 sub check_all_processed{

    my $var_dba    = shift;
    my $report  = shift;

    my $old_varfeat_ext_sth  = $var_dba->dbc->prepare(qq[ select count(*) from variation_feature]);
    my $new_varfeat_ext_sth  = $var_dba->dbc->prepare(qq[ select count(*) from variation_feature_working]);
    my $old_allele_ext_sth   = $var_dba->dbc->prepare(qq[ select count(*) from allele]);
    my $new_allele_ext_sth   = $var_dba->dbc->prepare(qq[ select count(*) from allele_working]);
  

    $old_varfeat_ext_sth->execute() || die "Failed to extract old_varfeat count \n";
    $new_varfeat_ext_sth->execute() || die "Failed to extract new_varfeat count\n";
    $old_allele_ext_sth->execute()  || die "Failed to extract old_allele count\n";   
    $new_allele_ext_sth->execute()  || die "Failed to extract new_allele count\n";


    my $old_varfeat = $old_varfeat_ext_sth->fetchall_arrayref();
    my $new_varfeat = $new_varfeat_ext_sth->fetchall_arrayref();
    my $old_allele  = $old_allele_ext_sth->fetchall_arrayref();
    my $new_allele  = $new_allele_ext_sth->fetchall_arrayref();

    my $all_ok = 0; ## can we proceed to rename tables?

   if($old_varfeat->[0]->[0] == $new_varfeat->[0]->[0]){
	print $report "Variation_Feature: Correct number of entries seen in variation_feature: $new_varfeat->[0]->[0]\n\n";
	$all_ok++;
    }
    else{
	print $report "Variation_Feature: ERROR old table has :$old_varfeat rows, new table has $new_varfeat->[0]->[0]\n\n";
    }

    if($old_allele->[0]->[0] == $new_allele->[0]->[0]){
	print $report "Allele: Correct number of entries seen : $new_allele->[0]->[0]\n";
	$all_ok++;
    }
    else{
	print $report "Allele: ERROR old table has : $old_allele->[0]->[0] rows, new table has: $new_allele->[0]->[0]\n";
    }

    return ( $new_allele->[0]->[0], $all_ok );

}


 sub get_failure_rates{

     my ($var_dba, $report, $allele_number) = @_;
    
    my $variation_ext_sth    = $var_dba->dbc->prepare(qq[ select count(*) from variation]);
    my $varfail_ext_sth      = $var_dba->dbc->prepare(qq[ select count(distinct variation_id) from failed_variation_working]);
    my $allelefail_ext_sth   = $var_dba->dbc->prepare(qq[ select count(*) from failed_allele_working]);


    $variation_ext_sth->execute()   || die "Failed to extract variation count\n";  
    $varfail_ext_sth->execute()     || die "Failed to extract varfail count\n";    
    $allelefail_ext_sth->execute()  || die "Failed to extract allelefail count\n"; 


    my $variation   = $variation_ext_sth->fetchall_arrayref();
    my $varfail     = $varfail_ext_sth->fetchall_arrayref();
    my $allelefail  = $allelefail_ext_sth->fetchall_arrayref();
 

     my $var_fail_rate = substr((100 * $varfail->[0]->[0] / $variation->[0]->[0] ), 0,5 ); 
     
     print $report "\nVariation failure rate:  $var_fail_rate % [$varfail->[0]->[0] / $variation->[0]->[0] ]\n";
     
     
     my $allele_fail_rate = substr((100 * $allelefail->[0]->[0] / $allele_number ), 0,5 ); 

     print $report "\nAllele failure rate: $allele_fail_rate % [$allelefail->[0]->[0] /$allele_number ] \n\n";
     
     
}

## Breakdown of fails by reason
sub check_failure_rates{

    my ($var_dba, $report) = @_;
          
    my $vardesc_ext_sth      = $var_dba->dbc->prepare(qq[ select fd.description,count(*) 
                                                         from failed_description fd, failed_variation_working fv
                                                         where fv.failed_description_id = fd.failed_description_id 
                                                         group by fd.description]);

    my $alleledesc_ext_sth   = $var_dba->dbc->prepare(qq[ select fd.description,count(*) 
                                                         from failed_description fd, failed_allele_working fa
                                                         where fa.failed_description_id = fd.failed_description_id 
                                                         group by fd.description]);


    $vardesc_ext_sth->execute()     || die "Failed to extract variation fail reasons\n";  
    $alleledesc_ext_sth->execute()  || die "Failed to extract allele fail reasons\n"; 

  
    my $vardesc     = $vardesc_ext_sth->fetchall_arrayref();
    my $alleledesc  = $alleledesc_ext_sth->fetchall_arrayref();



    print $report "\nVariation Failure reasons:\n";
    foreach my $l (@{$vardesc}){
	print $report "\t$l->[1]\t$l->[0]\n";
    }

    print $report "\nAllele Failure reasons:\n";
    foreach my $l (@{$alleledesc}){
	print $report "\t$l->[1]\t$l->[0]\n";
    }
}


sub rename_tables{

    my ($var_dba) = shift;

    ## Keep orignal tables in short term
    #$var_dba->dbc->do(qq[ rename table allele to allele_before_pp ]);
    #$var_dba->dbc->do(qq[ rename table variation_feature to variation_feature_before_pp ]);
    #$var_dba->dbc->do(qq[ rename table failed_allele to failed_allele_before_pp ]);         ## Not needed post dev phase
    #$var_dba->dbc->do(qq[ rename table failed_variation to failed_variation_before_pp ]);   ## Not needed post dev phase


    ## Rename working tables 
    #$var_dba->dbc->do(qq[ rename table allele_working to allele ]);
    #$var_dba->dbc->do(qq[ rename table variation_feature_working to variation_feature ]);
    #$var_dba->dbc->do(qq[ rename table failed_allele_working to failed_allele ]);        ## Not needed post dev phase
    #$var_dba->dbc->do(qq[ rename table failed_variation_working to variation ]);         ## Not needed post dev phase


}
1;
