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

=head1 DESCRIPTION

This module is used in the ehive variant quality control process. 
It runs checks on a new dbSNP import to identify obvious missing data
before the full QC process starts.
If problems are discovered, the job dies killing the hive process

=cut

package Bio::EnsEMBL::Variation::Pipeline::VariantQC::CheckdbSNPImport;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);
use Bio::EnsEMBL::Variation::Utils::QCUtils qw( count_rows count_for_statement get_evidence_attribs);

my $DEBUG = 0;



sub run {
   
    my $self = shift;

    my $dir = $self->required_param('pipeline_dir');

    ## don't check final tmp_sample_genotype_single_bp for human as merge table
    my @genotype_tables_to_check = ("population_genotype", "sample_genotype_multiple_bp");
    if($self->required_param('species') =~/homo|human/){
	push @genotype_tables_to_check, "tmp_sample_genotype_single_bp_subind_ch22";
    }
    else{
	push @genotype_tables_to_check, "tmp_sample_genotype_single_bp";
    }


    open my $report, ">", "$dir/preQC_report.txt" || die "Failed to open preQC_report.txt : $!\n";
    print $report "\n\tChecking pre-QC tables \n\n";

    my $var_dba   = $self->get_species_adaptor('variation');
    my $varfeat_count   = count_rows($var_dba,'variation_feature');
    my $variation_count = count_rows($var_dba,'variation');
    my $allele_count    = count_rows($var_dba,'allele');
    my $fail_count      = count_rows($var_dba,'failed_variation');

    my $mean_vf   = substr(( $varfeat_count / $variation_count),0,5);
    my $mean_al   = substr(( $allele_count / $variation_count),0,5);
    my $fail_rate = substr((100 * $fail_count / $variation_count),0,5);

#    my $var_no_ss_allele     = $self->count_no_ss_allele();
    my $var_no_allele_string = $self->count_no_allele_string();

    my $geno_no_sample       = $self->count_sampleless_geno(\@genotype_tables_to_check);
    my $geno_no_subsnp       = $self->count_geno_ss_problem(\@genotype_tables_to_check);
    my $geno_no_allele       = $self->count_alleleless_geno(\@genotype_tables_to_check);

    my $varfeat_no_pos       = $self->count_bad_varfeat();
    my $varfeat_no_seqreg    = $self->count_seq_region_problem();

    my $complimented_desc    = $self->check_complimented_desc();
    my $bad_position         = $self->check_bad_position();

    my $attribs_loaded       = $self->check_attribs();

    my $maf_loaded = 1;
    if($self->required_param('species') =~/homo|human/){
       $maf_loaded = $self->count_maf();
    }


    print $report "Post-import preQC check

Total Variation:        $variation_count
Total VariationFeature: $varfeat_count ( $mean_vf per variation )
Total Allele:           $allele_count ( $mean_al per variation )

Failed Variation:       $fail_count (failure rate: $fail_rate )


Variations without allele_string:   $var_no_allele_string\n";

#Variations without ss alleles:      $var_no_ss_allele

print $report "
Genotypes without samples:          $geno_no_sample
Genotypes without real ss:          $geno_no_subsnp
Genotypes without alleles:          $geno_no_allele

VariationFeature without start/end: $varfeat_no_pos 
VariationFeature without seqregion: $varfeat_no_seqreg
VariationFeature where end+1<start: $bad_position
\n";  

    print $report "\nERROR: missing evidence attribs\n\n" if $attribs_loaded == 0;

    print $report "\nERROR: missing MAF \n\n" if $maf_loaded == 0;

    print $report "ERROR: $complimented_desc complimented descriptions found - to be fixed manually\n\n" if $complimented_desc >0;

    if(
       $var_no_allele_string > 0 || 
       $variation_count  == 0    ||
       $varfeat_count    == 0    ||
       $allele_count     == 0    ||
       $geno_no_sample     >0    ||
       $varfeat_no_pos     >0    ||
       $varfeat_no_seqreg  >0    ||
       $attribs_loaded   == 0    ||
       $maf_loaded       == 0
       ){

        print $report "Exiting - missing data to import\n"; 
        die;  ## rest of QC process does not start if this fails
    }


}


## Checks for missing data
## No longer a fail criteria - only holding alleles with freq's for human
 sub count_no_ss_allele{

    my $self = shift;
    my $var_dba   = $self->get_species_adaptor('variation');

    my $no_allele_ext_stat  = (qq[ select count(*) from variation 
                                   left outer join  allele on allele.variation_id = variation.variation_id
                                   where allele.allele is null
                                ]);

    return count_for_statement($var_dba, $no_allele_ext_stat );

}

=head2 count_no_allele_string

  Look for variants without allele data

=cut
sub count_no_allele_string{

    my $self = shift;
    my $var_dba   = $self->get_species_adaptor('variation');

    my $no_allele_str_ext_stat  = (qq[ select count(*) from variation
                                     where variation_id not in (select variation_id  from allele_string)]);
 
    return count_for_statement($var_dba, $no_allele_str_ext_stat );

}


=head2 count_bad_varfeat

  Look for variants without coordinates

=cut
sub count_bad_varfeat{

    my $self = shift;
    my $var_dba   = $self->get_species_adaptor('variation');

    my $varfeat_ext_stat = (qq[ select count(*) from variation_feature where seq_region_start =0 or seq_region_end =0  ]);
    return count_for_statement($var_dba, $varfeat_ext_stat);

}
=head2 count_seq_region_problem

  Look for variants without sequence locations data

=cut
sub count_seq_region_problem{

    my $self = shift;
    my $var_dba   = $self->get_species_adaptor('variation');

    my $varfeat_ext_sth = (qq[ select count(*) from variation_feature 
                               left outer join  seq_region on ( variation_feature.seq_region_id =  seq_region.seq_region_id)
                               where  seq_region.name is null  ]);

    return count_for_statement($var_dba, $varfeat_ext_sth );

}


### Checks on genotype tables

=head2 count_sampleless_geno

  Look for genotypes without sample data

=cut
sub count_sampleless_geno{

    my $self   = shift;
    my $tables = shift;

    my $var_dba   = $self->get_species_adaptor('variation');

    my $tot = 0;

   foreach my $table ( @{$tables} ){ 

       # indiviudal population???
       my $sample = '';
       my $sample_id = '';
       if ($table =~ /individual/) {
           $sample = 'individual';
           $sample_id = 'individual_id';
       }
       elsif ($table =~ /sample/) {
           $sample = 'sample';
           $sample_id = 'sample_id';
       }  
       else {
          $sample = 'population';
          $sample_id = 'population_id';
       }
       my $no_sample_ext_sth  = $var_dba->dbc->prepare(qq[select count(*) from $table where $sample_id not in (select $sample_id from $sample) ]);
       $self->warning("At sampleless_geno with table  $table");
       $no_sample_ext_sth->execute() || die "Failed to extract no sample count for $table\n";
       my $no_sample_count    = $no_sample_ext_sth->fetchall_arrayref();

       if (defined $no_sample_count->[0]->[0]){
	   $tot += $no_sample_count->[0]->[0];
       }
   }

   return $tot;

}

=head2 count_alleleless_geno

  Look for genotypes without allele data

=cut
sub count_alleleless_geno{

    my $self   = shift;
    my $tables = shift;

    my $var_dba   = $self->get_species_adaptor('variation');

    my $total_problem_genotypes = 0;

   foreach my $table ( @{$tables} ){ 

     my $no_allele_ext_stat  = (qq[select count(*) from $table where allele_1 is null or allele_2 is null ]);

     my $no_allele_count    = count_for_statement($var_dba, $no_allele_ext_stat );

     if (defined $no_allele_count){
        $total_problem_genotypes += $no_allele_count;
      }
  }

  return $total_problem_genotypes;

}

=head2 count_sampleless_geno

  Look for genotypes without subsnp ids

=cut
sub count_geno_ss_problem{

    my $self   = shift;
    my $tables = shift;

    my $var_dba   = $self->get_species_adaptor('variation');

   ## useful (slow) SQL:
   ## select count(*) from tmp_sample_genotype_single_bp  where subsnp_id  not in (select subsnp_id from allele)

    my $total_problem = 0;

    my $max_subsnp_ext_sth      = $var_dba->dbc->prepare(qq[ select max(subsnp_id) from variation_synonym]);
    $max_subsnp_ext_sth->execute()||die "Failed to find max subsnp_id\n";
    my $max_subsnp = $max_subsnp_ext_sth->fetchall_arrayref();
    unless(defined  $max_subsnp->[0]->[0]){ die "Failed to find max subsnp_id\n";}

   foreach my $table ( @{$tables} ){ 
    
    my $geno_ext_sth = $var_dba->dbc->prepare(qq[ select count(*) from $table where subsnp_id >  $max_subsnp->[0]->[0] ]);

      $geno_ext_sth->execute()|| die "Failed to extract failed genotype subsnp count for $table\n";   
      my $geno_count  = $geno_ext_sth->fetchall_arrayref();
      $total_problem += $geno_count->[0]->[0];
    }

    return $total_problem;
}


=head2 check_complimented_desc

  Look for described alleles complimented in error
     eg "(318 BP DELETION)" =>  ")NOIAELEH PV 813("
  These can be complex, so are flagged for manual fixing

=cut
sub check_complimented_desc{

    my $self = shift;

    my $var_dba   = $self->get_species_adaptor('variation');

    my $data_ext_stat = (qq[ select count(*) from allele_string 
                            where allele_string like '%NOIAELEH%'
                            or allele_string like '%NOIAYESNI%']);

    return count_for_statement($var_dba , $data_ext_stat);
}

=head2 check_bad_position

  Look for variation feature positions where end coord + 1 is less than start coord
  either start = end or start = end + 1 for insertions

=cut
sub check_bad_position{

    my $self = shift;

    my $var_dba   = $self->get_species_adaptor('variation');

    my $data_ext_stat = (qq[ select count(*) from variation_feature 
                           where seq_region_start > seq_region_end + 1
                           ]);

    return count_for_statement($var_dba , $data_ext_stat);
}

=head2 check_attribs

  Check the expected evidence attribs are available.
  Check early to avoid partial jobs - not future-proof but will catch some problems


=cut
sub check_attribs{

    my $self = shift;

    my $var_dba = $self->get_species_adaptor('variation');
    my $attribs = get_evidence_attribs( $var_dba->dbc() );

    my $found_everything = 1;
    foreach my $ev ( "1000Genomes", "Cited", "ESP", "ExAC", "Frequency", "HapMap", "Multiple_observations","1000Bull_Genomes", "WTSI_MGP"){
        $found_everything = 0 unless defined $attribs->{$ev};
    }
    return $found_everything;
}


=head2 count_maf

  Look low MAF counts in human (none loading or partial loading of 1KG table)

=cut
sub count_maf{

    my $self = shift;
    my $var_dba   = $self->get_species_adaptor('variation');

    my $maf_ext_stat  = (qq[ select count(*) from maf]);

    my $maf =  count_for_statement( $var_dba, $maf_ext_stat );
    $maf > 80000000 ? return 1 : return 0;
}



1;
