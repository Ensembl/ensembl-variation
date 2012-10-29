
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
use Bio::EnsEMBL::Variation::Utils::dbSNP qw(get_alleles_from_pattern);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(revcomp_tandem);

use vars qw( @EXPORT_OK);
@EXPORT_OK = qw( check_four_bases check_illegal_characters remove_ambiguous_alleles find_ambiguous_alleles);

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

our %QUICK_COMP = ( "A" => "T",
                    "T" => "A",
                    "C" => "G",
                    "G" => "C"
    ); 

     
#ÊGet a string containing the possible ambiguity nucleotides
our $AMBIGUITIES = join("",keys(%AMBIG_REGEXP_HASH));

# Add the code for uracil in case some allele should have that
%AMBIG_REGEXP_HASH = (%AMBIG_REGEXP_HASH,('U' =>  'T'));



=head2 run

  Run checks on variant data; submit failed_variation and strand-corrected variation feature database updates then launch further allele checking

=cut
sub run {

  my $self = shift;
   
  

  ## variation_feature specific checks
  my ($var_data, $flip, $allele_string, $failed_variant, $failed_allele) = $self->run_variation_checks();

  ## allele specific checks 
  my ($allele_data, $failed_allele, $failed_ss) = $self->run_allele_checks( $flip, $failed_allele, $allele_string);  
 

  ## Database updates 
   my $var_dba = $self->get_species_adaptor('variation');

  ## write to new variation featues 
  insert_variation_features($var_dba, $var_data); 
 
  ## write to new failed_variation table
  write_variant_fails($var_dba, $failed_variant);

  ## write to new allele table
  $self->write_allele($var_dba, $allele_data);

  ### fail full experimental result for sample & subsnp_id
  write_allele_fails($var_dba , $failed_allele, $failed_ss );


  ## update variation table with flipped flag
  write_variant_flips($var_dba, $flip);

   
 
  
}

sub run_variation_checks{

  my $self = shift;

  ### start and end variation_id supplied
  my $first = $self->required_param('start_id');;
  my $last  = $first + $self->required_param('batch_size') -1;
  if($first ==1){$last--;}
  
  if( $DEBUG == 1){$self->warning("Starting to run variantQC with $first & $last " );}


  my %fail_variant;   # hash of arrays of failed variation_ids
  my %fail_allele;   # hash of arrays of arrays failed variation_ids & alleles
  my %flip;           # hash of variation ids which have their strand changed in this process
  my %allele_string;  # hash of expected allele for each variation_id strings saved for allele checking


  my $var_dba = $self->get_species_adaptor('variation');

  ## slice needed for ref check
  my $core_dba = $self->get_species_adaptor('core');
  my $slice_ad = $core_dba->get_SliceAdaptor;


  ## export current variation_feature data
  my $to_check = export_data_adding_allele_string($var_dba, $first, $last);
  


  foreach my $var (@{$to_check}){
  
    die "No allele string for $var->{name}\n" unless defined $var->{allele}; ## kill whole batch before any db updates
    ## save initial value of allele string for allele checking incase multi-mapping
    $allele_string{$var->{v_id}} = $var->{allele};

    ## Type 19 - fail variant if  >1 mapping seen
  
   if($var->{map} >1){       
      push @{$fail_variant{19}}, $var->{v_id};
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
    if( $var->{class_attrib_id} >0){
        ## named variants are given the SO term 'sequence alteration' on entry
        ## no further checking possible
	next;
    }

    # expand alleles if they contain brackets and numbers before other checks
    my $expanded = $var->{allele};
    expand(\$expanded);



  ##  Type 13 - non-nucleotide chars seen - flag variants & alleles as fails & don't run further checks

  my $illegal_alleles = check_illegal_characters($expanded);
  if(defined $illegal_alleles->[0] ) {
    push @{$fail_variant{13}}, $var->{v_id};

    foreach my $ill_al( @{$illegal_alleles} ){
      push @{$fail_allele{13}},  [$var->{v_id}, $ill_al];   ## unflipped allele failed throughout
    }
    next;  ## don't attempt to order variation_feature.allele_string or compliment
  }


  ## Type 3  flag variation as fail if it has [A/T/G/C] allele string 

  my $all_possible_check = check_four_bases($expanded);
  if ($all_possible_check ==1){        
    push @{$fail_variant{3}}, $var->{v_id}
  }

  ## Type 14 resolve ambiguities before reference check - flag variants & alleles as fails

  if($expanded =~ m/[$AMBIGUITIES]/){
    $expanded = remove_ambiguous_alleles(\$expanded);
    push @{$fail_variant{14}}, $var->{v_id};

    ## identify specific alleles to fail
    my $ambiguous_alleles = find_ambiguous_alleles($var->{allele});
    foreach my $amb_al( @{$ambiguous_alleles} ){
      push @{$fail_allele{14}},  [$var->{v_id}, $amb_al]; ## unflipped allele failed throughout
    }

  }
  ## Further checks only run for variants with 1 genomic location
  
  next if  $var->{map} > 1;

  
  ## flip allele string if on reverse strand and single mapping
  
    if( $var->{strand} eq "-1" ){
      reverse_comp(\$expanded );          ## for ref check
      if( $var->{allele}=~ /\(/){
          $var->{allele} = revcomp_tandem($var->{allele});
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
   ## save for allele checker if fully processed and possibly flipped
   $allele_string{$var->{v_id}} = $var->{allele};

  }
 

    return ($to_check, \%flip, \%allele_string, \%fail_variant, \%fail_allele) ;

}


=head2 run_allele_checks

  Check alleles in allele table match alleles in variation_feature table
      - catches frequency data with incorrect strand or no-calls counted in frequency calculations

  Submit failed_allele and strand-corrected allele updates 

=cut
sub run_allele_checks {

   my $self           = shift;
   my $flip           = shift; 
   my $fail           = shift;
   my $allele_string  = shift;

   ### start and end variation_id supplied
   my $first = $self->required_param('start_id'); 
   my $last  = $first + $self->required_param('batch_size') - 1;
   if($first ==1){$last--;} 

   my @fail_all;

   my $var_dba = $self->get_species_adaptor('variation');

   my %skip;  ## store sets to fail only once;

   ## export current allele & variation_feature data - alleles flipped where needed
   my $var_data = export_allele_data($var_dba, $first, $last, $flip );
   
   foreach my $var( keys %$var_data){
    
     my %expected_alleles;
     my @expected_alleles = split/\//, $allele_string->{$var};
     foreach my $exp(@expected_alleles){ 
       $expected_alleles{$exp} = 1;
     }
   

     ## check through allele submissions
     foreach my $submitted_data (@{$var_data->{$var}->{allele_data}}){

      
       ## $submitted_data content: [ al.allele_id, al.subsnp_id, al.allele,  al.frequency, al.sample_id, al.count, al.frequency_submitter_handle ]

       unless( exists $expected_alleles{$submitted_data->[2]} && $expected_alleles{$submitted_data->[2]} ==1 ){ ## check expected on allele not allele_code

         if(defined $submitted_data->[4]){
           ###fail whole experimental set if sample info supplied
           next if $skip{$submitted_data->[1]}{$submitted_data->[4]};    ## avoid duplicates
           push @fail_all, [$submitted_data->[1], $submitted_data->[4]]; ## subsnp_id, sample_id
         }
         else{
           push @{$fail->{11}}, [$var, $submitted_data->[2] ];          ## var_id, allele
         }
       }
     }
   }

   return ($var_data, $fail, \@fail_all);



}


=head2 export_data_adding_allele_string

  Extract variant data for checking taking allele_string from allele_string table 
  if required (dbSNP import pipeline does not populate variation_feature.allele_string)

=cut
sub export_data_adding_allele_string{

  my ($var_dba, $first, $last) = @_;
  
  my @to_check;

  my $variant_ext_sth = $var_dba->dbc->prepare(qq[SELECT vf.variation_id,
                                                      vf.variation_name,
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
                                                      sr.name,
                                                      vf.alignment_quality,
                                                      als.allele_string,
                                                      vf.validation_status
                                                 FROM variation_feature vf,
                                                      seq_region sr,
                                                      tmp_map_weight_working tmw,
                                                      allele_string als 
                                                 WHERE  vf.variation_id between ? and ? 
                                                 AND vf.variation_id = tmw.variation_id  
                                                 AND vf.seq_region_id = sr.seq_region_id
                                                 AND vf.variation_id = als.variation_id  
                                                 ]);       
 

 
  
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
    $save{class_attrib_id}  = $l->[12];
    $save{seqreg_name}      = $l->[13];
    $save{align_qual}       = $l->[14];
    $save{validation_status}= $l->[16];

    if($l->[15] =~ /^(\(.*\))\d+\/\d+/){## handle tandem
      my $expanded_alleles = get_alleles_from_pattern($l->[15]); 
      $save{allele} = join"/",@{$expanded_alleles};
    }
    else{
      $save{allele}         = $l->[15];
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

   my $data_ext_sth = $var_dba->dbc->prepare(qq[SELECT v.variation_id,
                                                       v.name,
                                                       al.allele_id,
                                                       al.subsnp_id,
                                                       al.allele,
                                                       al.frequency,
                                                       al.sample_id,
                                                       al.count,
                                                       al.frequency_submitter_handle
                                               FROM   variation v, allele al
                                               WHERE  v.variation_id between ? and ?
                                               AND    v.variation_id  = al.variation_id 
                                                ]);       
  
   $data_ext_sth->execute($first, $last)|| die "ERROR extracting variation feature info\n";
  
   my $data = $data_ext_sth->fetchall_arrayref();
   
   foreach my $l(@{$data}){  
    
       $save{$l->[0]}{name}   = $l->[1];
      
       if($flip->{$l->[0]}){
       ## update allele for flips
	   $l->[4]=~ /\(/   ?   $l->[4] = revcomp_tandem($l->[4]) :
	       defined $QUICK_COMP{$l->[4]}  ?  $l->[4] = $QUICK_COMP{$l->[4]} : 
	       reverse_comp(\$l->[4]);
	   
       }

       push @{$save{$l->[0]}{allele_data}}, [$l->[2], $l->[3], $l->[4], $l->[5], $l->[6], $l->[7], $l->[8] ];

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

  ## HGMD_MUTATION is permitted as an allele string
  $allele =~ s/\/|HGMD_MUTATION//g;
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

  my $fail_ins_sth = $var_dba->dbc->prepare(qq[insert  into failed_variation_working
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
  my $fail_list = shift;  ## has of arrays of alleles by fail class
  my $fail_all  = shift;  ## array of [subsnp_id, sample_id] 


  ## find allele id to fail by variant & allele  
  my $allele_id_ext_sth  = $var_dba->dbc->prepare(qq[ select al.allele_id 
                                                      from allele_working al, allele_code ac
                                                      where ac.allele =?
                                                      and  ac.allele_code_id = al.allele_code_id
                                                      and  al.variation_id =?
                                                   ]);
  
 
                                                   

  ## find allele id to fail as a submitted set with frequencies
  my $allele_set_id_ext_sth = $var_dba->dbc->prepare(qq[ select al.allele_id 
                                                         from   allele_working al
                                                         where  al.subsnp_id = ?
                                                         and    al.sample_id = ?
                                                       ]);



  my $fail_ins_sth   = $var_dba->dbc->prepare(qq[ insert into failed_allele_working
                                                  (allele_id, failed_description_id)
                                                  values (?,?)
                                                ]);       
          
  my %done;     ## save on old pk to drop out eariler
  my %done_new; ## save on new pk too

  ## deal with single allele fails
  foreach my $reason (keys %{$fail_list}){ 

    foreach my $var( @{$fail_list->{$reason}}  ){
      ## if the failure is a property of the reported allele, fail all entries with this allele

      if ($done{$reason}{$var->[0]}{$var->[1]}){ next; }  ## duplicates arise due to running on multiple variation_feature.allele_string's

      $done{$reason}{$var->[0]}{$var->[1]} = 1;
      
      $allele_id_ext_sth->execute($var->[1], $var->[0])|| die "ERROR extracting allele id info\n";
      my $allele_id = $allele_id_ext_sth->fetchall_arrayref();

       unless ($allele_id->[0]->[0]){
	   warn "Error finding allele id for var :  $var->[0] & allele  $var->[1] - not failing \n";
	   next;
       }
      foreach my $al(@{$allele_id}){
        $done_new{$al->[0]}{$reason} = 1; ## flag by reason and new id for screening later
	$fail_ins_sth->execute($al->[0], $reason)|| die "ERROR inserting allele fails info\n";
      }
    }
  }

  ### deal with sets of allele fails
  foreach my $set(@{$fail_all}){

      ## look up new pk for submission set on subsnp_id, sample_id 
      $allele_set_id_ext_sth->execute( $set->[0], $set->[1] )|| die "ERROR extracting allele id info\n";
      
      my $allele_set_ids = $allele_set_id_ext_sth->fetchall_arrayref();
      foreach my $allele_id (@{$allele_set_ids}){

        next if (exists $done_new{$allele_id->[0]}{11} && $done_new{$allele_id->[0]}{11} ==1);
        $done_new{$allele_id->[0]}{11} = 1;

        $fail_ins_sth->execute($allele_id->[0], 11 )|| die "ERROR inserting allele fails info\n";
      }
  }

}
=head2 insert_variation_features

  Create new version of variation feature tables with alleles ordered ref/alt where possible and complimented where necessary

=cut
sub insert_variation_features{
    
  my $var_dba  = shift;
  my $to_check = shift;

  my $varfeat_ins_sth = $var_dba->dbc->prepare(qq[insert  into variation_feature_working
                                                (variation_id, variation_name, seq_region_id,  seq_region_start, seq_region_end, seq_region_strand,  
                                                 allele_string, map_weight,  source_id, consequence_types, 
                                                 variation_set_id, somatic, class_attrib_id, alignment_quality,validation_status)
                                                 values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
                                                ]);       
  
  
  foreach my $data (@{$to_check}){ 

    die "No allele string for $data->{name}\n" unless defined  $data->{allele};
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
                            $data->{class_attrib_id},
                            $data->{align_qual},
			    $data->{validation_status}
                             )|| die "ERROR importing variation feature info\n";


  }
}


=head2 write_variant_flips

  Populate a tmp table with a list of variants whose data will be flipped

=cut
sub write_variant_flips{
    
  my $var_dba = shift;
  my $flip    = shift;
        

  my $flip_ins_sth = $var_dba->dbc->prepare(qq[ insert  into variation_to_reverse_working
                                               (variation_id) values ( ?)
                                              ]);       
  
  foreach my $var (keys %{$flip} ){ 

    #### updating VARIATION table at end - no change to imported data unless all OK
    $flip_ins_sth->execute($var)||    die "ERROR adding variation flip status\n";    

  }
}


=head2 write_allele

  Create new & old version of allele table with alleles complimented where neccessary
  Expects raw import format as input

=cut
sub write_allele{

  my $self      = shift;
  my $var_dba   = shift;
  my $var_data  = shift;

  my $code = read_allele_code($var_dba);

  my $allele_old_ins_sth = $var_dba->dbc->prepare(qq[ insert  into MTMP_allele_working
                                                    (variation_id, subsnp_id, allele, frequency, sample_id, count)
                                                    values (?, ?, ?, ?, ?,? )
                                                     ]);
 
  my $allele_new_ins_sth = $var_dba->dbc->prepare(qq[ insert  into allele_working
                                                     (variation_id, subsnp_id, allele_code_id, frequency, sample_id, count, frequency_submitter_handle)
                                                     values (?, ?, ?, ?, ?, ?, ? )
                                                    ]);
    

  foreach my $var (keys %$var_data){ ## look up missing allele codes before grabbing lock on table to update
    foreach my $allele (@{$var_data->{$var}->{allele_data}}){  ## array of data from 1 row in table
      ## allele data in coded format - insert if not available
      unless($code->{$allele->[2]}){
          $code->{$allele->[2]} = get_allele_code($var_dba, $allele->[2]);
      }
    }
  }

  foreach my $var (keys %$var_data){

    foreach my $allele (@{$var_data->{$var}->{allele_data}}){  ## array of data from 1 row in table

      ## keep allele data in un-coded format for Mart for non-human databases
       unless($self->required_param('species') =~/Homo|Human/i){
          $allele_old_ins_sth->execute( $var, $allele->[1], $allele->[2], $allele->[3], $allele->[4], $allele->[5]) || die "ERROR inserting allele info\n";
     }

      $allele_new_ins_sth->execute( $var, $allele->[1], $code->{$allele->[2]}, $allele->[3], $allele->[4], $allele->[5], $allele->[6]) || die "ERROR inserting allele info\n";
      #### DANGER - not restricting these to have the same PK
    }
  }
}




=head2 read_allele_code

  Description: Extract all allele codes held 
  Status     : At Risk
=cut

sub read_allele_code{


    my $var_dba   = shift; 
    my %allele_code;

    my $code_ext_sth = $var_dba->dbc->prepare(qq[select allele_code_id, allele
                                                 from allele_code          
                                                 ]);       
            
    $code_ext_sth->execute()||die "Error extracting allele_codes\n";
    my $list = $code_ext_sth->fetchall_arrayref();

    foreach my $code (@{$list}){ 
        $allele_code{$code->[1]} = $code->[0];
    }

    return \%allele_code ;
}

=head2 get_allele_code

  Description: Extract single allele code or enter if novel 
  Status     : At Risk
=cut
sub get_allele_code{


    my $var_dba   = shift; 
    my $allele    = shift;

    my $code_ext_sth = $var_dba->dbc->prepare(qq[select allele_code_id
                                                 from allele_code where allele =?          
                                                 ]);       
            
    $code_ext_sth->execute($allele )||die "Error extracting allele_codes\n";
    my $id = $code_ext_sth->fetchall_arrayref();

    unless(defined $id->[0]->[0]){  ## enter if unknown

       my $code_ins_sth = $var_dba->dbc->prepare(qq[insert ignore into allele_code (allele) values(?) ]);  # breaks with muliple processes without ignore     
       $code_ins_sth->execute($allele)||die "Error inserting allele_codes\n";;

       $code_ext_sth->execute($allele )||die "Error extracting allele_codes\n";
       $id = $code_ext_sth->fetchall_arrayref();
    }
    unless(defined $id->[0]->[0]){ die "Error fetching allele code for .$allele.\n";}
    return $id->[0]->[0] ;
}


### put somewhere sensible
sub unique {
  my %a;
  map { $a{$_} = 1; } @_;
  return sort keys %a;
}


1;
