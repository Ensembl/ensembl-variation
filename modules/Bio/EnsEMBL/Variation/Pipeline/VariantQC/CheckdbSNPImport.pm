
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

package Bio::EnsEMBL::Variation::Pipeline::VariantQC::CheckdbSNPImport;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

my $DEBUG = 0;

## This runs some basic checks on the raw dbSNP import before full QC

sub run {
   
    my $self = shift;

    my $dir = $self->required_param('pipeline_dir');

    open my $report, ">", "$dir/preQC_report.txt" || die "Failed to open preQC_report.txt : $!\n";
    print $report "\n\tChecking pre-QC tables \n\n";


    my $varfeat_count   = $self->count_variation_feature();
    my $variation_count = $self->count_variation();
    my $allele_count    = $self->count_allele();
    my $fail_count      = $self->count_failed_variation();

    my $mean_vf   = substr(( $varfeat_count / $variation_count),0,5);
    my $mean_al   = substr(( $allele_count / $variation_count),0,5);
    my $fail_rate = substr((100 * $fail_count / $variation_count),0,5);

    my $var_no_ss_allele     = $self->count_no_ss_allele();
    my $var_no_allele_string = $self->count_no_allele_string();
     
    print $report "Post-import preQC check

Total Variation:        $variation_count
Total VariationFeature: $varfeat_count ( $mean_vf per variation )
Total Allele:           $allele_count ( $mean_al per variation )

Failed Variation:       $fail_count (failure rate: $fail_rate )

Variations without ss alleles:    $var_no_ss_allele  
Variations without allele_string: $var_no_allele_string 

\n";  

    if($var_no_ss_allele   > 0 || $var_no_allele_string > 0){

	print $report "Exiting - missing data to import\n"; 
	die;  ## rest of QC process does not start if this fails
    }
}


sub count_variation_feature{

    my $self = shift;
    my $var_dba   = $self->get_species_adaptor('variation');

    my $varfeat_ext_sth  = $var_dba->dbc->prepare(qq[ select count(*) from variation_feature]);
    $varfeat_ext_sth->execute() || die "Failed to extract varfeat count \n";
    my $varfeat_count = $varfeat_ext_sth->fetchall_arrayref();

    defined $varfeat_count->[0]->[0]  ?  return $varfeat_count->[0]->[0] : return 0;
}


sub count_variation{

    my $self = shift;
    my $var_dba   = $self->get_species_adaptor('variation');

    my $var_ext_sth   = $var_dba->dbc->prepare(qq[ select count(*) from variation]);
    $var_ext_sth->execute() || die "Failed to extract var count\n";
    my $var_count     = $var_ext_sth->fetchall_arrayref();

    defined $var_count->[0]->[0]  ?  return $var_count->[0]->[0] : return 0;
}

sub count_allele{

    my $self = shift;
    my $var_dba   = $self->get_species_adaptor('variation');

    my $allele_ext_sth = $var_dba->dbc->prepare(qq[ select count(*) from allele]);
    $allele_ext_sth->execute()|| die "Failed to extract allele count\n";   
    my $allele_count   = $allele_ext_sth->fetchall_arrayref();

    defined $allele_count->[0]->[0]  ?  return $allele_count->[0]->[0] : return 0;
}


sub count_failed_variation{

    my $self = shift;
    my $var_dba   = $self->get_species_adaptor('variation');


    my $fail_ext_sth  = $var_dba->dbc->prepare(qq[ select count(*) from failed_variation]);
    $fail_ext_sth->execute() || die "Failed to extract fail count\n";   
    my $fail_count    = $fail_ext_sth->fetchall_arrayref();

    defined $fail_count->[0]->[0]  ?  return $fail_count->[0]->[0] : return 0;
}
 
 sub count_no_ss_allele{

    my $self = shift;
    my $var_dba   = $self->get_species_adaptor('variation');

    my $no_allele_ext_sth      = $var_dba->dbc->prepare(qq[ select count(*) from variation 
                                                            where variation_id not in (select variation_id  from allele)]);

    $no_allele_ext_sth->execute() || die "Failed to extract no_allele count\n";   
    my $no_allele_count    = $no_allele_ext_sth->fetchall_arrayref();

    defined $no_allele_count->[0]->[0]  ?  return $no_allele_count->[0]->[0] : return 0;
}

sub count_no_allele_string{

    my $self = shift;
    my $var_dba   = $self->get_species_adaptor('variation');

    my $no_allele_str_ext_sth  = $var_dba->dbc->prepare(qq[ select count(*) from variation
                                                            where variation_id not in (select variation_id  from allele_string)]);
 
    $no_allele_str_ext_sth->execute() || die "Failed to extract no allele string count\n";
    my $no_allele_str_count    = $no_allele_str_ext_sth->fetchall_arrayref();

    defined $no_allele_str_count->[0]->[0]  ?  return $no_allele_str_count->[0]->[0] : return 0;

 
}




1;
