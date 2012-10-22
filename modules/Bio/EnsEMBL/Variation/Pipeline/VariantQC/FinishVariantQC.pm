
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

=head1 DESCRIPTION

This module is used in the ehive variant quality control process. 
It runs checks to identify large-scale process failure and investigate data quality.
If all checks are passed, final database updates take place to complete the process.

=cut


package Bio::EnsEMBL::Variation::Pipeline::VariantQC::FinishVariantQC;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);


my $DEBUG = 0;


sub run {
   
    my $self = shift;

    my $dir = $self->required_param('pipeline_dir');

    open my $report, ">", "$dir/QC_report.txt" || die "Failed to open QC_report.txt : $!\n";
    print $report "\n\tChecking post-QC tables \n\n";

    my $var_dba   = $self->get_species_adaptor('variation');
    my $core_dba  = $self->get_species_adaptor('core');


    ## Have all rows been processed in flipped/ re-ordered tables
    my ($allele_number, $process_error )  = $self->check_all_processed($var_dba, $report);

    ## check for consistency between new and old variation_feature tables
    print  $report "\nRunning checks variation_feature data\n";
    my $vf_consistancy_fail = check_variation_feature_consistency($var_dba, $report);
    print  $report "\t- variation_feature consistency checks passed\n" if  $vf_consistancy_fail == 0;
     
    print $report "\nChecking failure rates\n";
    ## What are the failure rates for alleles and variants
    my $suspiciously_poor = get_failure_rates($var_dba, $report, $allele_number );
    
    ## crude check to ensure same MT sequence used
    check_MT_fails($var_dba, $report);


    ## what are failure reasons for alleles and variants    
    fails_by_type($var_dba, $report);
    

    ## run checks of flipping process 
    print  $report "\nRunning checks on flipped data:\n";
    my ($popgen_fail)  = check_flipping_population_genotype($var_dba,'population_genotype',$report);
    my ($allele_fail)  = check_flipping_allele($var_dba,'allele',$report);
    my ($varfeat_fail) = check_flipping_variation_feature($var_dba,'variation_feature',$report);


    ## check all statuses and exit if there is a problem
    if( $suspiciously_poor   ==1 ||   #  high failure rates
        $popgen_fail         ==1 ||   #  flipping error in population_genotype_working
        $allele_fail         ==1 ||   #  flipping error in allele_working
        $varfeat_fail        ==1 ||   #  flipping error in variation_feature_working
        $process_error       ==1 ||   #  not all rows processed into *working tables
        $vf_consistancy_fail ==1      #  data missmatch between old and new variation_feature tables
       ){
        print $report "\n\nExiting due errors - not renaming tables or running updates\n"; 
        die;
    }
    
   

    ### if the results of the checks look ok, update & rename tables 
    print  $report "\nRunning updates and renaming working tables:\n";

    rename_tables($var_dba);

    update_failed_variation_set($var_dba, $report);

    update_meta_coord( $core_dba, $var_dba, 'variation_feature');

}


=head2 check_all_processed

  Compare the number of rows in raw and processed tables to check for 
  incomplete processing

=cut
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

    my $process_error = 0; ## can we proceed to rename tables?

    print $report "Checking row counts between raw and post processed tables:\n"; 

    if($old_varfeat == $new_varfeat){
	print $report "\tVariation_feature:\t\tOK ($new_varfeat rows)\n";
    }
    else{
	$process_error = 1;
	print $report "Variation_Feature: ERROR old table has :$old_varfeat rows, new table has $new_varfeat\n";
    }

    if($old_allele == $new_allele){
	print $report "\tAllele:\t\t\t\tOK ($new_allele  rows)\n";
    }
    else{
	$process_error = 1;
	print $report "Allele: ERROR old table has : $old_allele rows, new table has: $new_allele\n";
    }

    if($old_pop_geno == $new_pop_geno){
	print $report "\tPopulation_genotype:\t\tOK ($new_pop_geno rows)\n";
    }
    else{
	$process_error = 1;
	print $report "Population_genotype: ERROR old table has : $old_pop_geno rows, new table has: $new_pop_geno\n";
    }

    unless($self->required_param('species') =~/Homo|Human/i){ ## Mart tables not required for human
      if($old_allele == $mart_allele){
         print $report "\tMart Allele:\t\t\tOK ($mart_allele rows)\n";
      }
      else{
         $process_error = 1;
         print $report "Allele: ERROR old table has : $old_allele rows, mart table has: $mart_allele\n";
      }

      if($old_pop_geno == $mart_pop_geno){
         print $report "\tMart population_genotype:\tOK ($mart_pop_geno rows)\n";
      }
      else{
         $process_error = 1;
         print $report "Population_genotype: ERROR old table has : $old_pop_geno rows, mart table has: $mart_pop_geno\n";
      }

    }
    return ( $new_allele, $process_error );

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

sub count_group_by{

   my $var_dba     = shift;
   my $table_name  = shift;
   my $column_name = shift;

   return unless defined  $table_name  && defined $column_name;

   my %count;
   my  $row_count_ext_sth  = $var_dba->dbc->prepare(qq[ select $column_name, count(*) from $table_name group by $column_name]);   

   $row_count_ext_sth->execute();
   my $data = $row_count_ext_sth->fetchall_arrayref();
   foreach my $l(@{$data}){
       $count{$l->[0]} = $l->[1];
   }

   return \%count; 

}

=head2 get_failure_rates

  Calculate the failure rates for variants and alleles

=cut
sub get_failure_rates{

    my ($var_dba, $report, $allele_number) = @_;
    
    my $variation_count        = count_rows($var_dba, 'variation');
    my $failed_variation_count = count_rows($var_dba, 'failed_variation_working','variation_id');
    my $failed_allele_count    = count_rows($var_dba, 'failed_allele_working');


    my $var_fail_rate = substr((100 * $failed_variation_count / $variation_count ), 0,5 ); 
     
    print $report "\tVariation failure rate:\t$var_fail_rate\% [$failed_variation_count / $variation_count ]\n";
     
     
    my $allele_fail_rate = substr((100 * $failed_allele_count / $allele_number ), 0,5 ); 

    print $report "\tAllele failure rate:\t$allele_fail_rate\% [$failed_allele_count /$allele_number ]\n";
     

    my $suspiciously_poor;
    if( $var_fail_rate >10 || $allele_fail_rate > 3){
         $suspiciously_poor = 1;
    }
    return  $suspiciously_poor;
}

=head2 fails_by_type

  Report the number of variants and alleles failing for each reason

=cut
sub fails_by_type{

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

=head2 check_MT_fails

  There are occaisional differences between the mitochondrial  
  reference sequence used by  dbSNP and ensembl. When this happens, 
  MT variants will have a high failure rate due to the variant 
  alleles not matching the reference at the reported position. 
  This subroutine checks the MT failure rate for this reason
  to catch unexpected reference diferences.

=cut
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
                                                       and failed_variation_working.failed_description_id  = 2
                                                      ]);
  $all_MT_ext_sth->execute()||die;
  my $count_all_MT =  $all_MT_ext_sth->fetchall_arrayref();

  if( $count_all_MT->[0]->[0] > 0){
  ## if the species has any MT variants, count the number failing

      $fail_MT_ext_sth->execute()||die;
      my $count_fail_MT =  $fail_MT_ext_sth->fetchall_arrayref();

      if( $count_fail_MT->[0]->[0] > 0){
          $fail_rate = substr((100 * $count_fail_MT->[0]->[0]/$count_all_MT->[0]->[0]),0,5);
      }
      print $report "\tMT variant failure rate:$fail_rate\% [$count_fail_MT->[0]->[0]/$count_all_MT->[0]->[0]]\n\n";
  }
  else{
      print $report "\nNo MT variants observed\n";
  }

}

=head2 check_variation_feature_consistency

  Check for consistency between new and old variation_feature tables 
  by comparing the number of rows with each possible value for the 
  validation_status, somatic and seq_region_id columns

=cut
sub check_variation_feature_consistency{

  my ($var_dba, $report) = shift;

  my $consistency_fail = 0;

  my @columns = ('validation_status', 'somatic', 'seq_region_id');

  foreach my $column (@columns){

    my $new_status = count_group_by($var_dba, 'variation_feature_working', $column );
    my $old_status = count_group_by($var_dba, 'variation_feature', $column);

    foreach my $status (keys %{$old_status}){
      unless ($old_status->{$status} == $new_status->{$status}){
        print $report "Consistency error: variation_feature.$column $status ($old_status->{$status} => $new_status->{$status})\n";
        $consistency_fail = 1;
      }
    }
  }
  return $consistency_fail;
}
 


=head2 check_flipping_population_genotype

  Extract a subset of data from raw and processed population genotype
  tables to check genotypes have been flipped where expected.

=cut
sub check_flipping_population_genotype{

  my $var_dba = shift;
  my $type    = shift;
  my $report  = shift;

  my $flip_check_stmt  = qq[select v.name, CONCAT(old.allele_1,'/',old.allele_2), CONCAT(alc1.allele, '/',alc2.allele) 
                             from population_genotype old, population_genotype_working new, genotype_code gc1,
                             allele_code alc1,genotype_code gc2, allele_code alc2, variation v, variation_to_reverse_working vtr 
                             where v.variation_id = vtr.variation_id 
                             and v.variation_id = old.variation_id 
                             and old.population_genotype_id = new.population_genotype_id 
                             and new.genotype_code_id  = gc1.genotype_code_id 
                             and gc1.haplotype_id = 1 
                             and alc1.allele_code_id = gc1.allele_code_id 
                             and new.genotype_code_id  = gc2.genotype_code_id 
                             and gc2.haplotype_id = 2 
                             and alc2.allele_code_id = gc2.allele_code_id 
                             limit 3000
                           ];


  my $flip_check_sth = $var_dba->dbc->prepare($flip_check_stmt);
  $flip_check_sth->execute()||die "Failed to run pop_gen_flip check\n";
  my $pop_gen_flip_dat = $flip_check_sth->fetchall_arrayref();
  unless( defined $pop_gen_flip_dat->[0]->[0]){
      print  $report "No flipped population genotypes available for checking\n";
      return 0;
  }
 
  return  check_flipping($pop_gen_flip_dat, $type, $report );

}

=head2 check_flipping_allele

  Extract a subset of data from raw and processed allele
  tables to check alleles have been flipped where expected.

=cut
sub check_flipping_allele{

  my $var_dba = shift;
  my $type    = shift;
  my $report  = shift;

  my $flip_check_stmt = qq[select distinct old.subsnp_id, old.allele, v.name 
                             from allele old, variation v, variation_to_reverse_working vtr  
                             where v.variation_id = vtr.variation_id
                             and v.variation_id = old.variation_id                              
                             limit 3000
                           ];

 
  my $flip_check_sth = $var_dba->dbc->prepare($flip_check_stmt);
  $flip_check_sth->execute()||die "Failed to run old allele_flip check \n";
  my $allele_flip_dat = $flip_check_sth->fetchall_arrayref();

  unless( defined $allele_flip_dat->[0]->[0]){
      print  $report "No flipped alleles available for checking\n";
      return 0;
  }

  my $allele_new_stmt = qq[ select distinct ac.allele 
                            from allele_working al, allele_code ac
                            where al.subsnp_id = ?
                            and ac.allele_code_id = al.allele_code_id
                           ];

  my $allele_new_sth = $var_dba->dbc->prepare($allele_new_stmt);



  my %data;
  foreach my $l(@{$allele_flip_dat}){
      $data{$l->[0]}{old} .= $l->[1] . "/";
      $data{$l->[0]}{name} = $l->[2];

      $allele_new_sth->execute($l->[0])||die "Failed to run new allele_flip check \n";
      my $als = $allele_new_sth->fetchall_arrayref();
      foreach my $al(@{$als}){
	  $data{$l->[0]}{new} .= $als->[0] . "/";
      }
  }
  my @cleaned_data;
  foreach my $ss (keys %data){
      push @cleaned_data, [$data{$ss}{name},  $data{$ss}{old}, $data{$ss}{new}]; 
  }

  return  check_flipping(\@cleaned_data, $type, $report );
 

}

=head2 check_flipping_variation_feature

  Extract a subset of data from allele_string and processed variation_feature
  tables to check allele_strings have been flipped where expected.
      - don't check strings which are the same reverse complimented (eg. -/AT => -/AT)
      - don't check failed variants which may not have been flipped 
=cut
sub check_flipping_variation_feature{

  my $var_dba = shift;
  my $type    = shift;
  my $report  = shift;

  my $flip_check_stmt = qq[select distinct new.variation_name, old.allele_string, new.allele_string 
                             from allele_string old, variation_feature_working new, variation v,variation_to_reverse_working vtr  
                             where v.variation_id = vtr.variation_id 
                             and v.variation_id = old.variation_id 
                             and new.variation_id = v.variation_id
                             and new.variation_id not in (select variation_id from failed_variation_working)
                             and new.allele_string not like '%-%'
                             limit 3000
                           ];



  my $flip_check_sth = $var_dba->dbc->prepare($flip_check_stmt);
  $flip_check_sth->execute()||die "Failed to run variation_feature flip check\n";
  my $var_feat_flip_dat = $flip_check_sth->fetchall_arrayref();

  unless( defined $var_feat_flip_dat->[0]->[0]){
      print  $report "No flipped variation_features available for checking\n";
      return 0;
  }


  ## filter out G/C and A/T and sort alleles as switched to match reference
  my @cleaned_data;
  foreach my $l (@{$var_feat_flip_dat}){
      next if $l->[1] =~/G\/C|A\/T|T\/A|C\/G/;
      push @cleaned_data, [$l->[0], sort($l->[1]), sort ($l->[2])];  
  }

   unless( scalar @cleaned_data > 100){
      print  $report "ERROR checking variation_feature flipping - no data to check post clean up\n";
      return 1;
  }

  return check_flipping(\@cleaned_data, $type, $report );
  

}
=head2 check_flipping

  Compare bases to identify un-flipped data

=cut
sub check_flipping{

  my $data   = shift;
  my $type   = shift;
  my $report = shift;

  my $is_fail = 0;
  my $total_checked = scalar @{$data};

  foreach my $l (@{$data}){
    if($l->[1] eq $l->[2]){
      $is_fail = 1;
      print $report "Flipping error ($type) $l->[0] is $l->[1] pre and $l->[2] post\n";
    }
  }

  if($is_fail ==0){print  $report "\tOK: $type\t[$total_checked entries checked]\n";}
  return ($is_fail);    

}





=head2 rename_tables

  Rename original tables to tablename_before_pp
  Rename tablename_working tables to tablename
  Sort variation_feature table for rapid web access

=cut
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

=head2 update_failed_variation_set

  Add all failed variations to the failed_variation set
  (requires the attrib tables to be populated)

=cut
sub update_failed_variation_set{

  my $var_dba = shift;
  my $report  = shift;
    
   
  my $failed_var_ext_sth  = $var_dba->dbc->prepare(qq[ select distinct variation_id
                                                       from failed_variation
                                                     ]);



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
  print $report "\n$check_num failedvariants inserted into variation_set_variation table\n";
    
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
