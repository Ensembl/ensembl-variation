
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

Bio::EnsEMBL::Variation::Pipeline::VariantQC::VariantQC

=head1 DESCRIPTION

Runs basic quality control on variant and allele data imported from extrenal sources

=cut

package Bio::EnsEMBL::Variation::Pipeline::VariantQC::VariantQC;


use strict;
use warnings;

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp expand);
use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);


our $DEBUG   = 1;

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



=head2 run

  Run checks on variant data; submit failed_variation and strand-corrected variation feature database updates then launch further allele checking

=cut
sub run {

  my $self = shift;
  
  ### start and end variation_id supplied
  my $first = $self->required_param('start_id'); 
  my $last  = $first + $self->required_param('batch_size') -1; 
  if($first ==1){$last--;} 

  if( $DEBUG == 1){$self->warning("Starting to run variantQC with $first & $last " );}

     
  #ÊGet a string containing the possible ambiguity nucleotides
  my $AMBIGUITIES = join("",keys(%AMBIG_REGEXP_HASH));

  # Add the code for uracil in case some allele should have that
  %AMBIG_REGEXP_HASH = (%AMBIG_REGEXP_HASH,('U' =>  'T'));


 
  my %fail_variant;   # hash of arrays of failed variation_ids
  my %fail_allele;    # hash of arrays of arrays failed variation_ids & alleles
  my %flip;           # hash of variation ids which have their strand changed in this process


  my $var_dba = $self->get_species_adaptor('variation');

  ## slice needed for ref check
  my $core_dba = $self->get_species_adaptor('core');
  my $slice_ad = $core_dba->get_SliceAdaptor;


  ## export current variation_feature data, adding allele_string if required 
  my $to_check;
  if($self->required_param('do_allele_string') == 1){
    $to_check = export_data_adding_allele_string($var_dba, $first, $last);
  }
  else{
    $to_check = export_data_with_allele_string($var_dba, $first, $last);
  }


  foreach my $var (@{$to_check}){
  

    ## Type 1 - fail variant if  >1 mapping seen
  
   if($var->{map} >1){       
      push @{$fail_variant{1}}, $var->{v_id};
    }


    ## Type 4 - novariation fails - flag variants & alleles as fails & don't run further checks

    if($var->{allele} =~ /NOVARIATION/){
      push @{$fail_variant{4}}, $var->{v_id};
      push @{$fail_allele{4}}, [$var->{v_id}, "NOVARIATION"];   ## unflipped allele failed throughout
      next;
    }
  
    if($var->{allele} =~ /HGMD_MUTATION/){
      ## only locations are available, so no further checking is possible.
      next;
    }


    # expand alleles if they contain brackets and numbers before other checks
    my $expanded = $var->{allele};
    expand(\$expanded);



  ##  Type 13 - non-nucleotide chars seen - flag variants & alleles as fails & don't run further checks

  my $illegal_alleles = check_illegal_characters($expanded,$AMBIGUITIES);
  if(defined $illegal_alleles->[0] ) {
    push @{$fail_variant{13}}, $var->{v_id};

    foreach my $ill_al( @{$illegal_alleles} ){
      push @{$fail_allele{13}},  [$var->{v_id}, $ill_al];   ## unflipped allele failed throughout
    }
    next;
  }


  ## Type 3  flag variation as fail if it has [A/T/G/C] allele string 

  my $all_possible_check = check_four_bases($expanded);
  if ($all_possible_check ==1){        
    push @{$fail_variant{3}}, $var->{v_id}
  }

  ## Type 14 resolve ambiguities before reference check - flag variants & alleles as fails

  if($expanded =~ m/[$AMBIGUITIES]/){
    $expanded = remove_ambiguous_alleles(\$expanded,$AMBIGUITIES);
    push @{$fail_variant{14}}, $var->{v_id};

    ## identify specific alleles to fail
    my $ambiguous_alleles = find_ambiguous_alleles($var->{allele},$AMBIGUITIES);
    foreach my $amb_al( @{$ambiguous_alleles} ){
      push @{$fail_allele{14}},  [$var->{v_id}, $amb_al]; ## unflipped allele failed throughout
    }

  }
  ## Further checks only run for variants with <3 map locations [Why not single]
  
  next if  $var->{map} > 1;

  
  ## flip allele string if on reverse strand and single mapping
  
    if( $var->{strand} eq "-1" ){
      reverse_comp(\$expanded );          ## for ref check
      if( $var->{allele}=~ /\(/){
          $var->{allele} = rev_tandem($var->{allele});
      }
      else{
          reverse_comp(\$var->{allele} );   ## for database storage
      }
      $var->{strand} = 1;

      ### store variation_id to use when flipping alleles
      $flip{$var->{v_id}} = 1;
    }
  
  
  
    # Extract reference sequence to run ref checks [ compliments for reverse strand multi-mappers]
  
    my $ref_seq = get_reference_base($var, $slice_ad) ;

    unless(defined $ref_seq){ 
      ## don't check further if obvious coordinate error
      push @{$fail_variant{15}}, $var->{v_id};
      next;
    }


    my $match_coord_length = 0; ## is either allele of compatible length with given coordinates?
  
    my @alleles = split/\//, $var->{allele} ; 
    foreach my $al(@alleles){
      my $ch = $al;
      expand(\$ch);  ## run ref check against ATATAT not (AT)3
      if($ch eq $ref_seq){
          $var->{ref} = $al;
          $match_coord_length = 1;
      }
      else{        
      $var->{alt} .= "/" . $al;
      ## if one of these suggests an insertion and the other a substitution, the sizes are not compatible
      if(length($ref_seq ) == length($ch) && $ch ne "-" && $ref_seq ne "-"){$match_coord_length = 1;}
      }
    }

    unless  ($match_coord_length == 1){
      push @{$fail_variant{15}}, $var->{v_id};
    }

    ## lengths ok - is actual base in agreement?
    if( defined $var->{ref}){
      ## re-order the allele_string field with the reference allele first 
      $var->{allele} = $var->{ref} . $var->{alt};            
    }
    else{
      ##  Type 2 - flag variants as fails if neither allele matches the reference
      push @{$fail_variant{2}}, $var->{v_id} ;
    }    
  }
  

  ## Database updates 
  
  ## write to new variation featues 
  insert_variation_features($var_dba, $to_check); 
 
  ## write to new failed_variation table
  write_variant_fails($var_dba, \%fail_variant);

  ## update variation table
  write_variant_flips($var_dba, \%flip);

   
  ## allele-specific checks - run on updated variation_feature table
  run_allele_checks($self, \%flip, \%fail_allele);  
  
}



=head2 run_allele_checks

  Run checks on allele table; submit submit failed_allele and strand-corrected allele updates 

=cut
sub run_allele_checks {

   my $self   = shift;
   my $flip   = shift; 
   my $fail   = shift;

   ### start and end variation_id supplied
   my $first = $self->required_param('start_id'); 
   my $last  = $first + $self->required_param('batch_size') - 1;
   if($first ==1){$last--;} 

   my %fail_all;

   my $var_dba = $self->get_species_adaptor('variation');


   my ($var_data, $allele_codes);
   ## export current allele & variation_feature data - alleles flipped where needed
   if( $self->required_param('schema') eq 'old'){ 
     $var_data = export_allele_data($var_dba, $first, $last, $flip);
   }
   else{
     $var_data = export_allele_data_new_schema($var_dba, $first, $last, $flip);
   }

   
   foreach my $var( keys %$var_data){

     ## cope with unmapped variants with alleles or named alleles, but don't check them
     next unless (defined $var_data->{$var}->{vf_allele} &&  $var_data->{$var}->{vf_allele} =~/\w+/); ## neaten

     ## isolate expected alleles 
     my %expected_alleles;
     my @expected_alleles = split/\//, $var_data->{$var}->{vf_allele};
     foreach my $exp(@expected_alleles){ 
       $expected_alleles{$exp} = 1;
     }


     ## check through allele submissions
     foreach my $submitted_data (@{$var_data->{$var}->{allele_data}}){

       my $check_allele = $submitted_data->[6];
      
       ## $submitted_data content: [ al.allele_id, al.subsnp_id, al.allele[code_id],  al.frequency, al.sample_id, al.count, al.allele ]

       unless( exists $expected_alleles{$check_allele} && $expected_alleles{$check_allele} ==1 ){ ## check expected on allele not allele_code

         if(defined $submitted_data->[4]){
           ###fail whole experimental set
           $fail_all{$submitted_data->[1]}{$submitted_data->[4]} = 1;   ## subsnp_id, sample_id
         }
         else{
           push @{$fail->{11}}, [$var, $submitted_data->[6] ];          ## var_id, allele
         }
       }
     }
   }
   ## write allele records to new table flipping strands as needed
   write_allele($var_dba, $var_data, $self->required_param('schema'));

   ### fail full experimental result for sample & subsnp_id
   write_allele_fails($var_dba , $fail, \%fail_all, $var_data,$self->required_param('schema') );
}


=head2 export_data_with_allele_string

  Extract variant data for checking taking allele_string from variation_feature.allele_string table 

=cut
sub export_data_with_allele_string{

   my ($var_dba, $first, $last) = @_;
   
   my @to_check;

   my $data_ext_sth = $var_dba->dbc->prepare(qq[SELECT v.variation_id,
                                                       v.name,
                                                       vf.variation_feature_id,
                                                       vf.seq_region_id,
                                                       vf.seq_region_start,
                                                       vf.seq_region_end,
                                                       vf.seq_region_strand,
                                                       vf.allele_string,
                                                       tmw.count,
                                                       vf.source_id,
                                                       vf.consequence_types,
                                                       vf.variation_set_id,
                                                       vf.somatic,
                                                       vf.class_attrib_id,
                                                       sr.name 
                                              FROM
                                                  variation v, 
                                                  variation_feature vf,
                                                  seq_region sr,
                                                  tmp_map_weight_working tmw
                                              WHERE
                                                  v.variation_id between ? and ?
                                                  and vf.variation_id = v.variation_id 
                                                  and vf.variation_id = tmw.variation_id  
                                                  and vf.seq_region_id = sr.seq_region_id
                                                ]);       
  
  $data_ext_sth->execute($first, $last)|| die "ERROR extracting variation feature info\n";
  
  my $data = $data_ext_sth->fetchall_arrayref();
  foreach my $l(@{$data}){
    my %save;

    $save{v_id}           = $l->[0];
    $save{name}           = $l->[1];
    $save{vf_id}          = $l->[2];
    $save{seqreg_id}      = $l->[3];
    $save{start}          = $l->[4];
    $save{end}            = $l->[5];
    $save{strand}         = $l->[6];
    $save{allele}         = $l->[7];
    $save{map}            = $l->[8];
    $save{source_id}      = $l->[9];
    $save{consequence_types} = $l->[10];
    $save{variation_set_id}  = $l->[11];
    $save{somatic}           = $l->[12];
    $save{class_attrib_id}   = $l->[13];
    $save{seqreg_name}       = $l->[14];

    push @to_check,\%save;;
  }

  return (\@to_check);  

}

=head2 export_data_adding_allele_string

  Extract variant data for checking taking allele_string from allele_string table 
  if required (dbSNP import pipeline does not populate variation_feature.allele_string)

=cut
sub export_data_adding_allele_string{

  my ($var_dba, $first, $last) = @_;
  
  my @to_check;

  my $variant_ext_sth = $var_dba->dbc->prepare(qq[SELECT v.variation_id,
                                                      v.name,
                                                      vf.variation_feature_id,
                                                      vf.seq_region_id,
                                                      vf.seq_region_start,
                                                      vf.seq_region_end,
                                                      vf.seq_region_strand,
                                                      tmw.count,
                                                      vf.source_id,
                                                      vf.consequence_types,
                                                      vf.variation_set_id,
                                                      vf.somatic,
                                                      vf.class_attrib_id,
                                                      sr.name 
                                             FROM
                                                 variation v, 
                                                 variation_feature vf,
                                                 seq_region sr,
                                                 tmp_map_weight_working tmw
                                             WHERE
                                                 v.variation_id between ? and ?
                                                 and vf.variation_id = v.variation_id 
                                                 and vf.variation_id = tmw.variation_id  
                                                 and vf.seq_region_id = sr.seq_region_id
                                               ]);       
 

  my $allele_ext_sth = $var_dba->dbc->prepare(qq[SELECT a.variation_id,
                                                        a.allele   
                                                FROM   allele_string a 
                                                WHERE  a.variation_id between ? and ?
                                               ]);       


  ### get allele information

  my %alleles;
  
  $allele_ext_sth->execute($first, $last)|| die "ERROR extracting allele feature info\n";

  my $allele_data = $allele_ext_sth->fetchall_arrayref();

  foreach my $l(@{$allele_data}){
    push @{$alleles{$l->[0]}},  $l->[1];  
  }
  
  ### get variant information
  
  $variant_ext_sth->execute($first, $last)|| die "ERROR extracting variation feature info\n";
  
  my $variant_data = $variant_ext_sth->fetchall_arrayref();
  foreach my $l(@{$variant_data}){
    my %save;
  
    $save{v_id}           = $l->[0];
    $save{name}           = $l->[1];
    $save{vf_id}          = $l->[2];
    $save{seqreg_id}      = $l->[3];
    $save{start}          = $l->[4];
    $save{end}            = $l->[5];
    $save{strand}         = $l->[6];


    $save{map}              = $l->[7];
    $save{source_id}        = $l->[8];
    $save{consequence_types}= $l->[9];
    $save{variation_set_id} = $l->[10];
    $save{somatic}          = $l->[11];

    $save{seqreg_name}      = $l->[13];
  
    if( defined $alleles{$l->[0]}->[0]){
      $save{allele}           = join '/', @{$alleles{$l->[0]}};
    }     
    else{
      $save{allele}           = "";
      warn "No alleles available for variant $l->[1]";
    } 

    push @to_check,\%save;;
  }

  return (\@to_check);  

}

=head2 export_allele_data

  Extract allele data and variation_feature allele string where available
  Compliments alleles if variant has unique map location and is on reverse strand

  Uses old variation schema for compatibility with import pipeline

=cut
sub export_allele_data{

   my ($var_dba, $first, $last, $flip) = @_;

   my %save;       
          
   my %done;
   ### Look up expected allele string from new variation_feature table with flipped allele strings
   my $data_ext_sth = $var_dba->dbc->prepare(qq[SELECT v.variation_id,
                                                       v.name,
                                                       vf.allele_string,
                                                       al.allele_id,
                                                       al.subsnp_id,
                                                       al.allele,
                                                       al.frequency,
                                                       al.sample_id,
                                                       al.count,
                                                       al.allele,
                                                       vf.variation_feature_id       
                                              FROM  variation v join allele al on(v.variation_id = al.variation_id )
                                                  left outer join variation_feature_working vf on (vf.variation_id = v.variation_id )
                                              WHERE
                                                  v.variation_id between ? and ?
                                                ]);       
  
   $data_ext_sth->execute($first, $last)|| die "ERROR extracting variation feature info\n";
  
   my $data = $data_ext_sth->fetchall_arrayref();
   
   foreach my $l(@{$data}){
   ## distinct on query causes database tmp to fill - handling hackily instead [only extracting variation_feature_id  for this purpose]
     next if $done{"$l->[3]\_$l->[10]"} ;
     $done{"$l->[3]\_$l->[10]"} = 1;
     $save{$l->[0]}{name}         = $l->[1];
     $save{$l->[0]}{vf_allele}    = $l->[2];
      
     if($flip->{$l->[0]}){

       ## update allele for flips
       reverse_comp(\$l->[5]);
       reverse_comp(\$l->[9]);   
     }

  push @{$save{$l->[0]}{allele_data}}, [$l->[3], $l->[4], $l->[5], $l->[6], $l->[7], $l->[8], $l->[9] ];

  }

  return (\%save);  

}

=head2 export_allele_data_new_schema

  Extract allele data and variation_feature allele string where available
  Compliments alleles if variant has unique map location and is on reverse strand

  Uses new variation schema 

=cut
sub export_allele_data_new_schema{

  my ($var_dba, $first, $last, $flip, $schema) = @_;
  
  my %save;

  
  # hack to use with old schema
  my %allele_codes;  ## save as looked up - should be quicker than holding all unusual alleles 

  #my $allele_code_ext_sth  = $var_dba->dbc->prepare(qq[ select allele_code_id from allele_code  where allele = ? ]);       
  
  #my $allele_code_ins_sth  = $var_dba->dbc->prepare(qq[ insert into allele_code (allele) values (?) ]);
      
  my %done;
  ### Look up expected allele string from new variation_feature table with flipped allele strings
  my $data_ext_sth = $var_dba->dbc->prepare(qq[SELECT v.variation_id,
                                                      v.name,
                                                      vf.allele_string,
                                                      al.allele_id,
                                                      al.subsnp_id,
                                                      al.allele_code_id,
                                                      al.frequency,
                                                      al.sample_id,
                                                      al.count,
                                                      alc.allele,
                                                      vf.variation_feature_id       
                                             FROM  variation v join allele al on(v.variation_id = al.variation_id )
                                                   join allele_code alc on (al.allele_code_id = alc.allele_code_id )
                                                   left outer join variation_feature_working vf on (vf.variation_id = v.variation_id )
                                             WHERE
                                                 v.variation_id between ? and ?
                                               ]);       
  
  $data_ext_sth->execute($first, $last)|| die "ERROR extracting allele info\n";
  
  my $data = $data_ext_sth->fetchall_arrayref();
  
  foreach my $l(@{$data}){
    ## distinct on query causes database tmp to fill - handling hackily instead [only extracting variation_feature_id  for this purpose]
    next if $done{"$l->[3]\_$l->[10]"} ;
    $done{"$l->[3]\_$l->[10]"} = 1;
    $save{$l->[0]}{name}         = $l->[1];
    $save{$l->[0]}{vf_allele}    = $l->[2];
        
    if($flip->{$l->[0]}){

        ## update allele for flips
        reverse_comp(\$l->[9]);
        unless ($allele_codes{$l->[9]}){
             $allele_codes{$l->[9]} = get_allele_code($var_dba, $l->[9] );
        }
        ## updated allele code for flip
        $l->[5] = $allele_codes{$l->[9]} ;
    }

    push @{$save{$l->[0]}{allele_data}}, [$l->[3], $l->[4], $l->[5], $l->[6], $l->[7], $l->[8], $l->[9] ];

    }

  return (\%save);  

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

  my $allele      = shift;
  my $AMBIGUITIES = shift;

  ## HGMD_MUTATION is permitted as an allele string
  $allele =~ s/\/|HGMD_MUTATION//g;
  if ($allele =~ m /[^ACGTU\-$AMBIGUITIES]/i){
    ## identify specfic alleles to flag
    my @fail;
    my @al = split/\//, $allele;
    foreach my $al(@al){ 
      if ($al =~ m /[^ACGTU\-$AMBIGUITIES]/i){
        push @fail, $al;
      }
    }
    return \@fail;
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
  my $AMBIGUITIES   = shift;

  $allele_string =~ s/([U$AMBIGUITIES])/$AMBIG_REGEXP_HASH{$1}/ig;

  return $allele_string;
}

=head2 find_ambiguous_alleles

  Checks if any of the bases in an allele string are ambiguity codes
  Returntype : reference to array of ambiguous alleles

=cut
sub find_ambiguous_alleles{

  my $allele_string = shift;
  my $AMBIGUITIES   = shift;

  my @fail;
  my @al = split/\//, $allele_string;
  foreach my $al(@al){ 
    if ($al =~ m /$AMBIGUITIES/i){
      push @fail, $al;
    }
  }
  return \@fail;

}
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


=head2 write_variant_fails

    Update failed_variation_working with all variation_id/reason combinations

=cut
sub write_variant_fails{


  my $var_dba   = shift; 
  my $fail_list = shift;

  my $fail_ins_sth = $var_dba->dbc->prepare(qq[insert into failed_variation_working
                                               (variation_id, failed_description_id)
                                               values (?,?)
                                               ]);       
          

  foreach my $reason (keys %{$fail_list}){ 

    ## duplicates arise due to running on variation_features not variations
    my @fails = unique(@{$fail_list->{$reason}});    

    foreach my $var( @fails  ){
      $fail_ins_sth->execute($var, $reason)|| die "ERROR inserting variation fails info\n";
    }
  }
}



=head2 write_allele_fails

    Update failed_allele_working with all variation_id/reason combinations


=cut
sub write_allele_fails{


  my $var_dba   = shift; 
  my $fail_list = shift;
  my $fail_all  = shift;  
  my $var_data  = shift;
  my $schema    = shift;


  my $allele_id_ext_sth ;
  ## find allele id to fail - cope with different schemas
  if($schema =~ /old/){
    $allele_id_ext_sth  = $var_dba->dbc->prepare(qq[select allele_id 
                                                   from allele_working
                                                   where allele =?
                                                   and allele_working.variation_id =?
                                                   ]);
  }
  else{
    $allele_id_ext_sth  = $var_dba->dbc->prepare(qq[select allele_id 
                                                   from allele_working, allele_code
                                                   where allele_code.allele =?
                                                   and allele_code.allele_code_id = allele_working.allele_code_id
                                                   and allele_working.variation_id =?
                                                   ]);
  }
 
                                                   

  ## those failing as a submitted set with frequencies
  my $allele_set_id_ext_sth = $var_dba->dbc->prepare(qq[select allele_id 
                                                       from  allele_working
                                                       where allele_working.variation_id =?
                                                       and allele_working.subsnp_id =?
                                                       and allele_working.sample_id = ?
                                                    ]);



  my $fail_ins_sth   = $var_dba->dbc->prepare(qq[insert into failed_allele_working
                                                (allele_id, failed_description_id)
                                                values (?,?)
                                                ]);       
          
  my %done;     ## save on old pk to drop out eariler
  my %done_new; ## save on new pk too

  foreach my $reason (keys %{$fail_list}){ 

    foreach my $var( @{$fail_list->{$reason}}  ){
      ## if the failure is a property of the reported allele, fail all entries with this allele

      if ($done{$reason}{$var->[0]}{$var->[1]}){ next; }  ## duplicates arise due to running on multiple variation_feature.allele_string's

      $done{$reason}{$var->[0]}{$var->[1]} = 1;
      
      $allele_id_ext_sth->execute($var->[1], $var->[0])|| die "ERROR extracting allele id info\n";
      my $allele_id = $allele_id_ext_sth->fetchall_arrayref();

       if ($allele_id->[0]->[0]){
         $done_new{$allele_id->[0]->[0]}{$reason} = 1; ## flag by reason and new id for screening later
        $fail_ins_sth->execute($allele_id->[0]->[0], $reason)|| die "ERROR inserting allele fails info\n";
      }
      else{
        warn "Error finding allele id for var :  $var->[0] & allele  $var->[1] - not failing \n";
      }        
  }
  
  
  foreach my $var( keys %$var_data){
      
      ## check through allele submissions & fail set of results for same subsnp & sample
    foreach my $submitted_data (@{$var_data->{$var}->{allele_data}}){
      
      next unless defined $submitted_data->[4];  ## only fails with sample info need be checked
      
      next unless defined $$fail_all{$submitted_data->[1]}{$submitted_data->[4]};  ## failed by ssid & sample        
      
      ## look up new pk for submission set on variation_id, subsnp_is, sample_id [quicker in bulk, but only update once]
      $allele_set_id_ext_sth->execute($var, $submitted_data->[1], $submitted_data->[4] )|| die "ERROR extracting allele id info\n";
      
      my $allele_set_ids = $allele_set_id_ext_sth->fetchall_arrayref();
      foreach my $allele_id (@{$allele_set_ids}){

        next if (exists $done_new{$allele_id->[0]}{11} && $done_new{$allele_id->[0]}{11} ==1);
        $done_new{$allele_id->[0]}{11} = 1;

        $fail_ins_sth->execute($allele_id->[0], 11 )|| die "ERROR inserting allele fails info\n";
      }
      undef $$fail_all{$submitted_data->[1]}{$submitted_data->[4]};    
      }
    }    
  }
}
=head2 insert_variation_features

  Create new version of variation feature tables with alleles ordered ref/alt where possible and complimented where necessary

=cut
sub insert_variation_features{
    
  my $var_dba  = shift;
  my $to_check = shift;

  my $varfeat_ins_sth = $var_dba->dbc->prepare(qq[insert into variation_feature_working
                                                (variation_id, variation_name, seq_region_id,  seq_region_start, seq_region_end, seq_region_strand,  
                                                 allele_string, map_weight,  source_id, consequence_types variation_set_id, somatic, class_attrib_id)
                                                 values (?,?,?,?,?,?,?,?,?,?,?,?,?)
                                                ]);       
  
  
  foreach my $data (@{$to_check}){ 
  
    $varfeat_ins_sth->execute($data->{v_id},
                            $data->{name},
                            $data->{seqreg_id},
                            $data->{start},
                            $data->{end},
                            $data->{strand},
                            $data->{allele},
                            $data->{map},
                            $data->{source_id},
                            $data->{consequence_types},
                            $data->{variation_set_id},
                            $data->{somatic},
                            $data->{class_attrib_id})|| die "ERROR importing variation feature info\n";


  }
}



sub write_variant_flips{
    
  my $var_dba = shift;
  my $flip    = shift;
  
  #my $flip_update_sth = $var_dba->dbc->prepare(qq[ update variation
  #                                                set flipped = 1
  #                                                where variation_id = ?
  #                                                ]);       

  my $flip_ins_sth = $var_dba->dbc->prepare(qq[ insert into variation_to_reverse_working
                                               (variation_id) values ( ?)
                                              ]);       
  

  foreach my $var (keys %{$flip} ){ 

    #### updating VARIATION table at end - no change to imported data unless all OK
    #$flip_update_sth->execute($var)|| die "ERROR updating variation flip status\n";   ## doing this at the end in bulk  
    $flip_ins_sth->execute($var)||    die "ERROR adding variation flip status\n";    

  }
}

=head2 insert_variation_features

  Create new version of allele table with alleles complimented where neccessary

=cut
sub write_allele{
    ### write data to new allele table flipping where needed; best way to recover partial updates caused by fails
    
  my ($var_dba, $var_data, $schema) = @_;
  my %done;

  my $allele_ins_sth ;
  if($schema =~ /old/){
    $allele_ins_sth     = $var_dba->dbc->prepare(qq[ insert into allele_working
                                                   (variation_id, subsnp_id, allele, frequency, sample_id, count)
                                                   values (?, ?, ?, ?, ?,? )
                                                    ]);
  }
  else{
    $allele_ins_sth     = $var_dba->dbc->prepare(qq[ insert into allele_working
                                                   (variation_id, subsnp_id, allele_code_id, frequency, sample_id, count)
                                                   values (?, ?, ?, ?, ?,? )
                                                    ]);
  }
  
  foreach my $var (keys %$var_data){

    foreach my $allele (@{$var_data->{$var}->{allele_data}}){  ## array of data from 1 row in table
      next if exists $done{$allele->[0]} ; #&& $done{$allele->[0]} ==1;  ## filtering on old pk
      $done{$allele->[0]} = 1;

      $allele_ins_sth->execute( $var, $allele->[1], $allele->[2], $allele->[3], $allele->[4], $allele->[5]) || die "ERROR inserting allele info\n";
    }
  }
}


sub get_allele_code{

  my ($var_dba, $current_al) = @_;

  my $allele_code;

  my $allele_code_ext_sth  =  $var_dba->dbc->prepare(qq [select allele_code_id from allele_code where allele =?]);
  my $allele_code_ins_sth  =  $var_dba->dbc->prepare(qq [insert into allele_code (allele) values (?) ]);


  $allele_code_ext_sth->execute($current_al) || die "Failed to look up allele code for $current_al\n";;

  my $al_code = $allele_code_ext_sth->fetchall_arrayref();
  if(defined $al_code->[0]->[0]){
    $allele_code = $al_code->[0]->[0];
  }
  else{
    ## insert new allele code - may be needed due to complimenting long odd things
    warn "Entering allele $current_al\n";

    $allele_code_ins_sth->execute($current_al) || die "Failed to enter allele code for $current_al\n";

    $allele_code_ext_sth->execute($current_al) || die "Failed to look up allele code for $current_al\n";;

    my $al_code = $allele_code_ext_sth->fetchall_arrayref();

    if(defined $al_code->[0]->[0]){
        $allele_code = $al_code->[0]->[0];
    }
    else{
        die "New allele code insertion failed oddly\n";
    }
  }
  return $allele_code;
}




##reverse_comp (AC)17/(AC)19  =>  91)GT(/71)GT(
sub rev_tandem{

    
  my $allele_string = shift;

  my $new_allele_string;

  my @parts = split/\//, $allele_string;
  
  foreach my $part (@parts){
  
    if( $part =~/\d+$/){
      my $num = $&;
      $part =~ s/\d+$|\(|\)//g;
      
      reverse_comp(\$part);
      print "doing $part\n";
      $new_allele_string .= "(" . $part .")" . $num . "/";
    }
    else{
      reverse_comp(\$part);
      $new_allele_string .= $part. "/";;
    }
  }
  $new_allele_string =~ s/\/$//;
  
  return $new_allele_string ;
}

### put somewhere sensible
sub unique {
  my %a;
  map { $a{$_} = 1; } @_;
  return sort keys %a;
}


1;
