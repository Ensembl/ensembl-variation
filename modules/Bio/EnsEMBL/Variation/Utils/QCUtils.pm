
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

=head1 DESCRIPTION

This module contains functions used in the variant quality control process. 

=cut


package Bio::EnsEMBL::Variation::Utils::QCUtils;

use strict;
use warnings;

use base qw(Exporter);


our @EXPORT_OK = qw(check_four_bases get_reference_base check_illegal_characters check_for_ambiguous_alleles remove_ambiguous_alleles find_ambiguous_alleles summarise_evidence count_rows count_group_by );


 
our %AMBIG_REGEXP_HASH = (
    'M' =>  '[AC]',
    'R' =>  '[AG]',
    'W' =>  '[AT]',
    'S' =>  '[CG]',
    'Y' =>  '[CT]',
    'K' =>  '[GT]',
    'V' =>  '[ACG]',
    'H' =>  '[ACT]',
    'D' =>  '[AGT]',
    'B' =>  '[CGT]',
    'X' =>  '[ACGT]',
    'N' =>  '[ACGT]'
    );

     
#ÊGet a string containing the possible ambiguity nucleotides
our $AMBIGUITIES = join("",keys(%AMBIG_REGEXP_HASH));

# Add the code for uracil in case some allele should have that
%AMBIG_REGEXP_HASH = (%AMBIG_REGEXP_HASH,('U' =>  'T'));



=head2 check_four_bases

  Checks if all 4 bases are present in allele string
  Returntype : 1 if present, 0 if not

=cut
sub check_four_bases{

  my $allele_string =shift;

  my @alleles = split /\//, $allele_string ;

  return 0 if scalar(@alleles) < 4;


  my %allele;
  foreach my $al(@alleles){
      $allele{$al} = 1;
  } 
  
  if( exists $allele{A} && 
    exists $allele{T} && 
    exists $allele{C} && 
    exists $allele{G} ){
    return 1;
  }
  else{
    return 0;
  }
}


=head2 get_reference_base

  Extract sequence from genomic reference at variant coordinates
  to check against supplied coordinates

=cut

sub get_reference_base{

   my ($var, $slice_ad) = @_;
   
   my $ref_seq;
   
   if( ($var->{end} +1) == $var->{start}){ ## convention for insertions to reference
     $ref_seq = "-";
   }

   elsif( $var->{end} < $var->{start}){ ## coordinate error
     warn "Incorrect coords $var->{start} - $var->{end}  for $var->{name} \n";    
   }    

   else{
   
     # retrieve the reference sequence at that mapping for deletion or substitution    

     my $slice = $slice_ad->fetch_by_region('toplevel', $var->{seqreg_name}, $var->{start}, $var->{end});
   
     unless (defined $slice){ die "ERROR Getting slice for $var->{seqreg_name}, $var->{start}, $var->{end}";}
     $ref_seq = $slice->seq();

     # correct for multi-mapping variants which may be on negative strand
     if($var->{strand} eq "-1"){ reverse_comp(\$ref_seq);}
   }

   return $ref_seq;

}
=head2 check_illegal_characters

  Checks for non ambiguity code/ATGC character to fail

=cut
sub check_illegal_characters{

    my $allele   = shift;

  if ($allele =~ m /[^ACGTU\-$AMBIGUITIES]/i){
    ## identify specfic alleles to flag
    my @fail;
    my @al = split/\//, $allele;
    foreach my $al(@al){ 
      if ($al =~ m /[^ACGTU\-$AMBIGUITIES]|\)\)$/i){
        push @fail, $al;
      }
    }
    return \@fail;
  }
  else{
    return undef;
  }
}

=head2 check_for_ambiguous_alleles

  Expand ambiguous alleles to A/T/C/G
  Returntype : expanded allele string

=cut
sub check_for_ambiguous_alleles{

  my $allele_string = shift;

  if($allele_string =~ m/[$AMBIGUITIES]/){
      return 1;
  }
  else{

      return undef;
  }
}

=head2 remove_ambiguous_alleles

  Expand ambiguous alleles to A/T/C/G
  Returntype : expanded allele string

=cut
sub remove_ambiguous_alleles{

  my $allele_string = shift;

  $allele_string =~ s/([U$AMBIGUITIES])/$AMBIG_REGEXP_HASH{$1}/ig;

  return $allele_string;
}

=head2 find_ambiguous_alleles

  Checks if any of the bases in an allele string are ambiguity codes
  Returntype : reference to array of ambiguous alleles

=cut
sub find_ambiguous_alleles{

  my $allele_string = shift;

  my @fail;
  my @al = split/\//, $allele_string;
  foreach my $al(@al){ 

    if ($al =~ m/[$AMBIGUITIES]/i){
      push @fail, $al;
    }
  }
  return \@fail;

}




=head2 summarise_evidence

  Description: Compile evidence of variant quality into a single comma seperated string
  Status     : At Risk
=cut
sub summarise_evidence{

    my $var_dbh = shift;
    my $species = shift;
    my $first   = shift;
    my $last    = shift;

    my %evidence;


    ## summarise ss information
    my $ss_variations =  get_ss_variations($var_dbh, $first, $last);


    ## extract list of variants with pubmed citations
    ## not using this to not fail variants assuming it is quicker to fail then revert in bulk
    my $pubmed_variations = get_pubmed_variations($var_dbh, $first, $last);


    ## get 1KG discovered variants (human only)
    my $kg_variations =  get_KG_variations($var_dbh, $first, $last) 
	if $species =~/Homo|Human/i ;


    foreach my $var(keys %$ss_variations ){

      ## dbSNP ss submissions
      $evidence{$var} = "Multiple_observations,"  if defined $ss_variations->{$var}->{count} && $ss_variations->{$var}->{count} > 1;

      $evidence{$var} .= "Frequency,"             if defined $ss_variations->{$var}->{'freq'};


      ## additional human evidence 
      $evidence{$var} .= "HapMap,"                if defined $ss_variations->{$var}->{'HM'}  ;
      
      $evidence{$var} .= "1000Genomes,"           if defined $kg_variations->{$var} ;

      ## pubmed citations
      $evidence{$var} .= "Cited,"                 if defined $pubmed_variations->{$var}; 

      ## remove trailing commas
      $evidence{$var} =~ s/\,$//                  if defined $evidence{$var} ;

    }
    return \%evidence;

}
=head2  get_ss_variations

  Summarise information from dbSNP ss submissions:

      - check the number of independant observations
             - ie. different submitter handle or different population
      - check for any allele frequency information
      - check for allele frequency information from Hapap if human

=cut

sub get_ss_variations{

    my $var_dbh = shift;
    my $first   = shift;
    my $last    = shift;
    
    my %evidence;

    my $obs_var_ext_sth  = $var_dbh->prepare(qq[ select al.variation_id,
                                                        h.handle,
                                                        al.population_id,
                                                        al.frequency,
                                                        p.name,
                                                        al.subsnp_id,
                                                        al.count
                                                 from   subsnp_handle h, allele al
                                                 left outer join population p on ( p.population_id = al.population_id)
                                                 where  al.variation_id between ? and ?
                                                 and    h.subsnp_id = al.subsnp_id ]);
    $obs_var_ext_sth->execute($first, $last );
    my $dat = $obs_var_ext_sth->fetchall_arrayref();

    my %save_by_var; #  gather ss entries by variant

    foreach my $l (@{$dat}){

	$evidence{$l->[0]}{'obs'}  = 1;  ## save default value for each variant - full set to loop through later.

	$l->[2] = "N" unless defined $l->[2];

        #save  submitter handle, population and ss id to try to discern independent submissions
	push  @{$save_by_var{$l->[0]}}, [  $l->[1], $l->[2], $l->[5] ];


       ## Save frequency evidence for variant by variant id - ensure at least 2 chromosomes assayed and variant poly

       if(defined $l->[3] && $l->[3] < 1 && $l->[3] > 0 && $l->[6] >1){
            ## flag if frequency data available
            $evidence{$l->[0]}{'freq'}  = 1;

	    ## special case for human only
	    $evidence{$l->[0]}{'HM'}  = 1 if defined $l->[4] && $l->[4]   =~/HapMap/i;

	}
    }


    foreach my $var (keys %save_by_var){

	$evidence{$var}{count} = count_ss(\@{$save_by_var{$var}} );
    }
    return \%evidence;
}

sub count_ss{


    my $submissions = shift;

    my %submitter;
    my %ss_ids;
    my $total_count;


    foreach my $submission (@{$submissions}){

        ## group ss ids
        $ss_ids{$submission->[2]} = 1;

        ## group populations seen by submitter
        push @{$submitter{$submission->[0]}{pop}} , $submission->[1];
        $submitter{$submission->[0]}{no_pop} = 1  if $submission->[1] eq "N";        
    }
    ## count ss ids
    my $ss_ids_held = scalar keys %ss_ids;

    ## count uniq submitter/population combinations  ( could be 2 submissions from re-analysis of same data)
    foreach my $submitter (keys %submitter){
        ## count by submitter ignoring submissions without populations if submissions with populations are present
        my %unique_sub;
        map { $unique_sub{$_} = 1; } @{$submitter{$submitter}{pop}};
        my $submitter_count = scalar(keys %unique_sub);
        $submitter_count-- if $submitter_count >1 && defined $submitter{$submitter}{no_pop} && $submitter{$submitter}{no_pop} ==1;
        $total_count += $submitter_count;

    }
    ## return lowest count
    $ss_ids_held  < $total_count ? return $ss_ids_held : return $total_count;
}

## cited variants are not failed

sub get_pubmed_variations{


    my $var_dbh = shift;
    my $first   = shift;
    my $last    = shift;
 
    my %pubmed_variations;

    my $pubmed_var_ext_sth  = $var_dbh->prepare(qq[ select variation_id from variation_citation ]);

    $pubmed_var_ext_sth->execute();
    my $data = $pubmed_var_ext_sth->fetchall_arrayref();
 
    foreach my $l (@{$data}){
	$pubmed_variations{$l->[0]} = 1;
    }

    return \%pubmed_variations;

}


## human only - extact ids for variants found in 1000 genomes project

sub get_KG_variations{


  my $var_dbh = shift;
  my $first   = shift;
  my $last    = shift;
    
  my %kg_variations;

  my $var_ext_sth  = $var_dbh->prepare(qq[ select variation.variation_id 
                                           from variation, tmp_1kg_var 
                                           where variation.snp_id =  tmp_1kg_var.rs_id 
                                           and variation.variation_id  between ? and ?
                                         ]);
  $var_ext_sth->execute($first, $last);
  my $data = $var_ext_sth->fetchall_arrayref();
 
  foreach my $l (@{$data}){
      $kg_variations{$l->[0]} = 1;
  }

  return \%kg_variations;

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

1;
