
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

    my $var_dba   = $self->get_species_adaptor('variation');
    my $core_dba  = $self->get_species_adaptor('core');


    ## Have all rows been processed in flipped/ re-ordered table
    my ($allele_number, $all_ok)  = $self->check_all_processed($var_dba, $report);

     
    ## What are the failure rates for alleles and variants
    my $suspiciously_poor = get_failure_rates($var_dba, $report, $allele_number );
    
    ## what are failure reasons for alleles and variants    
    check_failure_rates($var_dba, $report);
    
    ## crude check to ensure same MT sequence used
    check_MT_fails($var_dba, $report);


    if( defined $suspiciously_poor ){
        print $report "\nExiting due to high failure rates - not renaming tables or running updates\n"; 
        die;
    }
    
    unless( $all_ok == 1){        
        print $report "\nExiting - not all rows processed - not renaming tables or running updates\n";
        die; 
    }

    ### if the results of the checks look ok, update & rename tables 
    rename_tables($var_dba);

    update_failed_variation_set($var_dba, $report);

    update_meta_coord( $core_dba, $var_dba, 'variation_feature');

}



 ## Have all rows been processed 
 sub check_all_processed{

    my $self    = shift;
    my $var_dba = shift;
    my $report  = shift;

   
    my $old_varfeat   = count_rows($var_dba, 'variation_feature');
    my $new_varfeat   = count_rows($var_dba, 'variation_feature_working');

    my $old_allele    = count_rows($var_dba, 'allele');
    my $new_allele    = count_rows($var_dba, 'allele_working');
    my $mart_allele   = count_rows($var_dba, 'MTMP_allele_working');

    my $new_pop_geno  = count_rows($var_dba, 'population_genotype_working');
    my $old_pop_geno  = count_rows($var_dba, 'population_genotype');
    my $mart_pop_geno = count_rows($var_dba, 'MTMP_population_genotype_working');

    my $all_ok = 1; ## can we proceed to rename tables?

    print $report "Checking row counts between raw and post processed tables\n\n"; 

    if($old_varfeat == $new_varfeat){
       print $report "Variation_feature: OK ($new_varfeat  rows)\n";
    }
    else{
       $all_ok = 0;
       print $report "Variation_Feature: ERROR old table has :$old_varfeat rows, new table has $new_varfeat\n";
    }

    if($old_allele == $new_allele){
       print $report "Allele: OK ($new_allele)\n";
    }
    else{
       $all_ok = 0;
       print $report "Allele: ERROR old table has : $old_allele rows, new table has: $new_allele\n";
    }

    if($old_pop_geno == $new_pop_geno){
         print $report "Population_genotype: OK ($new_pop_geno  rows)\n";
      }
      else{
         $all_ok = 0;
         print $report "Population_genotype: ERROR old table has : $old_pop_geno rows, new table has: $new_pop_geno\n";
    }

    unless($self->required_param('species') =~/Homo|Human/i){ ## Mart tables not required for human
      if($old_allele == $mart_allele){
         print $report "Mart Allele: OK ($mart_allele rows)\n";
      }
      else{
         $all_ok = 0;
         print $report "Allele: ERROR old table has : $old_allele rows, mart table has: $mart_allele\n";
      }

      if($old_pop_geno == $mart_pop_geno){
         print $report "Mart population_genotype: OK ($mart_pop_geno rows)\n";
      }
      else{
         $all_ok = 0;
         print $report "Population_genotype: ERROR old table has : $old_pop_geno rows, mart table has: $mart_pop_geno\n";
      }

    }
    return ( $new_allele, $all_ok );

}

sub count_rows{

   my $var_dba     = shift;
   my $table_name  = shift;
   my $column_name = shift;

   my $row_count_ext_sth;

   if(defined $column_name){
     $row_count_ext_sth  = $var_dba->dbc->prepare(qq[ select count(distinct $column_name) from $table_name]);
   }
   else{
     $row_count_ext_sth  = $var_dba->dbc->prepare(qq[ select count(*) from $table_name]);
   }

   $row_count_ext_sth->execute();
   my $total_rows = $row_count_ext_sth->fetchall_arrayref();

   return $total_rows->[0]->[0]; 

}

 sub get_failure_rates{

    my ($var_dba, $report, $allele_number) = @_;
    
    my $variation_count        = count_rows($var_dba, 'variation');
    my $failed_variation_count = count_rows($var_dba, 'failed_variation_working','variation_id');
    my $failed_allele_count    = count_rows($var_dba, 'failed_allele_working');


    my $var_fail_rate = substr((100 * $failed_variation_count / $variation_count ), 0,5 ); 
     
    print $report "\nVariation failure rate:\t$var_fail_rate\% [$failed_variation_count / $variation_count ]\n";
     
     
    my $allele_fail_rate = substr((100 * $failed_allele_count / $allele_number ), 0,5 ); 

    print $report "\nAllele failure rate:\t$allele_fail_rate\% [$failed_allele_count /$allele_number ] \n\n";
     

    my $suspiciously_poor;
    if( $var_fail_rate >10 || $allele_fail_rate > 1){
         $suspiciously_poor = 1;
    }
    return  $suspiciously_poor;
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

sub check_MT_fails{

  my $var_dba= shift;
  my $report= shift;

  my $fail_rate = 0;
 
  my $all_MT_ext_sth      = $var_dba->dbc->prepare(qq[ select count(distinct variation_id) 
                                                        from  variation_feature,seq_region 
                                                        where seq_region.name = 'MT' 
                                                        and variation_feature.seq_region_id = seq_region.seq_region_id
                                                        ]);

  my $fail_MT_ext_sth      = $var_dba->dbc->prepare(qq[select count(distinct variation_feature.variation_id) 
                                                       from  variation_feature,seq_region,failed_variation_working
                                                       where seq_region.name = 'MT' 
                                                       and variation_feature.seq_region_id = seq_region.seq_region_id 
                                                       and failed_variation_working.variation_id = variation_feature.variation_id
                                                      ]);
  $all_MT_ext_sth->execute()||die;
  my $count_all_MT =  $all_MT_ext_sth->fetchall_arrayref();

  if( $count_all_MT->[0]->[0] > 0){

      $fail_MT_ext_sth->execute()||die;
      my $count_fail_MT =  $fail_MT_ext_sth->fetchall_arrayref();


      if( $count_fail_MT->[0]->[0] > 0){

          $fail_rate = substr((100 * $count_fail_MT->[0]->[0]/$count_all_MT->[0]->[0]),0,5);
      }

      print $report "\nMT failure rate:\t$fail_rate\% [$count_fail_MT->[0]->[0]/$count_all_MT->[0]->[0]]\n\n";
  }
  else{
      print $report "\nNo MT variants observed\n\n";
  }

}
  
sub rename_tables{

  my ($var_dba) = shift;

  ## Capture dbSNP 'suspect' variant fails 
  $var_dba->dbc->do(qq[ insert into failed_variation_working (variation_id, failed_description_id) select variation_id, failed_description_id from failed_variation ]); 

  ## Keep orignal tables in short term
  $var_dba->dbc->do(qq[ rename table allele to allele_before_pp ]) || die;
  $var_dba->dbc->do(qq[ rename table variation_feature to variation_feature_before_pp ]) || die;
  $var_dba->dbc->do(qq[ rename table failed_allele to failed_allele_before_pp ]) || die;         ## Not needed post dev phase
  $var_dba->dbc->do(qq[ rename table failed_variation to failed_variation_before_pp ]) || die;   ## Not needed post dev phase


  # Rename working tables 
  $var_dba->dbc->do(qq[ rename table allele_working to allele ]) || die;
  $var_dba->dbc->do(qq[ rename table MTMP_allele_working to MTMP_allele ]) || die;

  $var_dba->dbc->do(qq[ rename table variation_feature_working to variation_feature ]) || die;

  $var_dba->dbc->do(qq[ alter table variation_feature order by seq_region_id,seq_region_start,seq_region_end ]) || die;

  $var_dba->dbc->do(qq[ rename table failed_allele_working to failed_allele ]) || die;        ## Not needed post dev phase
  $var_dba->dbc->do(qq[ rename table failed_variation_working to failed_variation ]) || die;  ## Not needed post dev phase

  $var_dba->dbc->do(qq[ rename table population_genotype to population_genotype_before_pp]) || die;
  $var_dba->dbc->do(qq[ rename table MTMP_population_genotype_working to MTMP_population_genotype]) || die;
  $var_dba->dbc->do(qq[ rename table population_genotype_working to population_genotype]) || die;


  # does this need binning?
  $var_dba->dbc->do(qq[ update variation set flipped = 0  ]) || die;
  $var_dba->dbc->do(qq[ update variation set flipped = 1 where variation_id in (select variation_id from variation_to_reverse_working) ]) || die;

  

}

sub update_failed_variation_set{

  my $var_dba = shift;
  my $report  = shift;
    
   
  my $failed_var_ext_sth  = $var_dba->dbc->prepare(qq[ select distinct variation_id
                                                       from failed_variation
                                                     ]);


  ### sql to get variation set id
  my $fail_attrib_ext_sth  = $var_dba->dbc->prepare(qq[ select attrib_id
                                                         from attrib
                                                         where attrib_type_id = 9
                                                         and value = 'fail_all'

                                                       ]);
 
   my $variation_set_ext_sth  = $var_dba->dbc->prepare(qq[ select variation_set_id
                                                            from variation_set
                                                            where name = ?
                                                           ]);

   my $variation_set_ins_sth  = $var_dba->dbc->prepare(qq[ insert into variation_set
                                                            ( name, description, short_name_attrib_id)
                                                            values (?,?,?)
                                                          ]);


  my $var_set_ins_ext_sth  = $var_dba->dbc->prepare(qq[insert into variation_set_variation
                                                         (variation_id, variation_set_id)
                                                         values (?,?)
                                                        ]);


  $variation_set_ext_sth->execute('All failed variations')  || die "Failed to extract allele fail reasons\n"; 
  my $failed_set_id = $variation_set_ext_sth->fetchall_arrayref();

  unless(defined $failed_set_id->[0]->[0] ){
      #no set entered - look up attrib for short name and enter set

    $fail_attrib_ext_sth->execute() || die "Failed to extract failed set attrib reasons\n";
    my $attrib = $fail_attrib_ext_sth->fetchall_arrayref();

    if(defined $attrib->[0]->[0] ){
      $variation_set_ins_sth->execute('All failed variations',
                                       'Variations that have failed the Ensembl QC checks' ,
                                        $attrib->[0]->[0] )|| die "Failed to insert failed set\n";       
    }
    else{
        die "Exiting: Error - attribs not loaded. Load attribs then re-run FinishVariantQC\n";
    }
  }

  $variation_set_ext_sth->execute('All failed variations')  || die "Failed to extract allele fail reasons\n"; 
  $failed_set_id = $variation_set_ext_sth->fetchall_arrayref();


  ## extract list of variation ids to put in set
  $failed_var_ext_sth->execute() || die "Failed to extract variation fails\n";  
  my $all_failed_var = $failed_var_ext_sth->fetchall_arrayref();

  unless( defined $all_failed_var->[0]->[0]){
    print $report "ERROR: no variation fails available to update to variation set \n";
    return;
  }

  my $check_num = scalar @{ $all_failed_var};
  foreach my $var(@{ $all_failed_var}){

    $var_set_ins_ext_sth->execute( $var->[0], $failed_set_id->[0]->[0] )||die;
  }
  print $report "\n$check_num variants inserted into variation_set_variation table\n";
    
}




#
# updates the meta coord table 
# copied from old pipeline - $csname  not supplied => only ever updated as chromosome
sub update_meta_coord {

  my $core_dba   = shift;
  my $var_dba    = shift;
  my $table_name = shift;

  my $csname     = 'chromosome';

  my $csa = $core_dba->get_CoordSystemAdaptor();

  my $cs = $csa->fetch_by_name($csname);

  my $sth = $var_dba->dbc->prepare
    ('INSERT IGNORE INTO meta_coord set table_name = ?, coord_system_id = ?, max_length =500');

  $sth->execute($table_name, $cs->dbID());

  $sth->finish();

  return;
}


1;
