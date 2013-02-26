
=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/legal/code_licence.html

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
use Data::Dumper ;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp expand);
use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);
use Bio::EnsEMBL::Variation::Utils::dbSNP qw(get_alleles_from_pattern);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(revcomp_tandem);
use Bio::EnsEMBL::Variation::Utils::QCUtils qw( check_four_bases get_reference_base check_illegal_characters check_for_ambiguous_alleles remove_ambiguous_alleles find_ambiguous_alleles summarise_evidence count_rows count_group_by);


our $DEBUG   = 1;


our %QUICK_COMP = ( "A" => "T",
                    "T" => "A",
                    "C" => "G",
                    "G" => "C"
    ); 



=head2 run

  Run checks on variant data; submit failed_variation and strand-corrected variation feature database updates then launch further allele checking

=cut
sub run {

  my $self = shift;
   

  ## variation_feature specific checks
  my ($var_data, $flip, $allele_string, $failed_variant, $failed_allele) = $self->run_variation_checks();

  ## allele specific checks 
  my ($allele_data, $failed_allele, $failed_ss) = $self->run_allele_checks( $flip, $failed_allele, $allele_string);  
 
  my $var_dba = $self->get_species_adaptor('variation');

  my $first = $self->required_param('start_id');
  my $last  = $first + $self->required_param('batch_size') -1;
  if($first ==1){$last--;}

  ## find supporting evidence & merge to single string
  my $evidence_summary = summarise_evidence($var_dba->dbc(),
                                             $self->required_param('species'),
                                             $first,
                                             $last);



  ## Database updates 

  ## write to new variation featues 
  write_variation_features($var_dba, $var_data, $evidence_summary ); 
 
  ## write to new failed_variation table
  write_variant_fails($var_dba, $failed_variant);

  ## write to new allele table
  my $allele_ids = $self->write_allele($var_dba, $allele_data);

  ### fail full experimental result for sample & subsnp_id
  write_allele_fails($var_dba , $allele_ids, $failed_allele, $failed_ss );


  ## write to new variation with evidence string & flip status
  my $minor_allele = $self->write_variation(  $evidence_summary , $flip );
  

  ## check dbSNP MAF against variation_feature allele string
  ## warn, but don't fix automatically (1KG genotypes loaded later)
  if($self->required_param('species') =~/homo/i){
     $self->check_minor_allele($minor_allele,  $allele_string); 
  }

}

sub run_variation_checks{

  my $self = shift;

  ### start and end variation_id supplied
  my $first = $self->required_param('start_id');;
  my $last  = $first + $self->required_param('batch_size') -1;
  if($first ==1){$last--;}
  
  if( $DEBUG == 1){$self->warning("Starting to run variantQC with $first & $last " );}


  my %fail_variant;   # hash of arrays of failed variation_ids
  my %fail_allele;    # hash of arrays of arrays failed variation_ids & alleles
  my %flip;           # hash of variation ids which have their strand changed in this process
  my %allele_string;  # hash of expected allele for each variation_id strings saved for allele checking

  


  my $var_dba = $self->get_species_adaptor('variation');

  ## slice needed for ref check
  my $core_dba = $self->get_species_adaptor('core');
  my $slice_ad = $core_dba->get_SliceAdaptor;


  ## export current variation_feature data
  my $to_check = export_data_adding_allele_string($var_dba, $first, $last);


  my $failed_set_id = find_failed_variation_set_id($var_dba);
  


  foreach my $var (@{$to_check}){
  
    die "No allele string for $var->{name}\n" unless defined $var->{allele}; ## kill whole batch before any db updates
    ## save initial value of allele string for allele checking incase multi-mapping
    $allele_string{$var->{v_id}} = $var->{allele};

   
    ## Type 19 - fail variant if  >1 mapping seen
  
   if($var->{map} >1){       
      push @{$fail_variant{19}}, $var->{v_id};
      $var->{variation_set_id} = $failed_set_id ;
    }


    ## Type 4 - novariation fails - flag variants & alleles as fails & don't run further checks

    if($var->{allele} =~ /NOVARIATION/){
      push @{$fail_variant{4}}, $var->{v_id};
      push @{$fail_allele{4}}, [$var->{v_id}, "NOVARIATION"];   ## unflipped allele failed throughout
      $var->{variation_set_id} = $failed_set_id ;
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

  my $illegal_alleles = check_illegal_characters($expanded);
  if(defined $illegal_alleles->[0] ) {

      ## named variants are given the SO term 'sequence alteration' on entry - do not fail these 
      next if $var->{class_attrib_id} eq 18;

      push @{$fail_variant{13}}, $var->{v_id};
      $var->{variation_set_id} = $failed_set_id ;

      foreach my $ill_al( @{$illegal_alleles} ){
         push @{$fail_allele{13}},  [$var->{v_id}, $ill_al];   ## unflipped allele failed throughout
    }
    next;  ## don't attempt to order variation_feature.allele_string or compliment
  }


  ## Type 3  flag variation as fail if it has [A/T/G/C] allele string 

  my $all_possible_check = check_four_bases($expanded);
  if ($all_possible_check ==1){        
      push @{$fail_variant{3}}, $var->{v_id};
      $var->{variation_set_id} = $failed_set_id ;
  }

  ## Type 14 resolve ambiguities before reference check - flag variants & alleles as fails

    my $is_ambiguous = check_for_ambiguous_alleles($expanded );
    if(defined $is_ambiguous){
	$expanded = remove_ambiguous_alleles(\$expanded);
	push @{$fail_variant{14}}, $var->{v_id};
	$var->{variation_set_id} = $failed_set_id ;

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
      $var->{variation_set_id} = $failed_set_id ;
      next;
    }


    my $match_coord_length = 1; ## is either allele of compatible length with given coordinates?
  
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
      }
    }

    unless  ($match_coord_length == 1){
	push @{$fail_variant{15}}, $var->{v_id};
	$var->{variation_set_id} = $failed_set_id ;
    }

    ## lengths ok - is actual base in agreement?
    if( defined $var->{ref}){
      ## re-order the allele_string field with the reference allele first 
      $var->{allele} = $var->{ref} . $var->{alt};            
    }
    else{
      ##  Type 2 - flag variants as fails if neither allele matches the reference
      push @{$fail_variant{2}}, $var->{v_id} ;
      $var->{variation_set_id} = $failed_set_id ;
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


   my @fail_all;

   my %skip;      ## store sets to fail only once;


   ## export current allele & variation_feature data - alleles flipped where needed
   my $var_data = $self->export_allele_data( $flip );
   
   foreach my $var( keys %$var_data){
    
     my %expected_alleles;
     my @expected_alleles = split/\//, $allele_string->{$var};
     foreach my $exp(@expected_alleles){ 
       $expected_alleles{$exp} = 1;
     }
   

     ## check through allele submissions
     foreach my $submitted_data (@{$var_data->{$var}->{allele_data}}){


       ## Apply QC checks
      
       ## $submitted_data content: [ al.allele_id, al.subsnp_id, al.allele,  al.frequency, al.sample_id, al.count, al.frequency_submitter_handle ]
       unless( exists $expected_alleles{$submitted_data->[2]} && $expected_alleles{$submitted_data->[2]} ==1 ){ ## check expected on allele not allele_code

       ## checking ss id here because of the way strains are handled in import
         if(defined $submitted_data->[4] && defined $submitted_data->[1] && $submitted_data->[1] >0  ){
           ###fail whole experimental set if sample info supplied
           next if defined $skip{$submitted_data->[1]}{$submitted_data->[4]} && $skip{$submitted_data->[1]}{$submitted_data->[4]} ==1;  ## avoid duplicates
	   $skip{$submitted_data->[1]}{$submitted_data->[4]} =1;
           push @fail_all, [$submitted_data->[1], $submitted_data->[4]]; ## subsnp_id, sample_id
         }
         else{
           push @{$fail->{11}}, [$var, $submitted_data->[2] ];          ## var_id, allele
         }
       }
     }
   }

   return ($var_data, $fail, \@fail_all );



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

   my $self           = shift;
   my $flip           = shift; 


   ### start and end variation_id supplied
   my $first = $self->required_param('start_id'); 
   my $last  = $first + $self->required_param('batch_size') - 1;
   if($first ==1){$last--;} 

   my $var_dba = $self->get_species_adaptor('variation');


   my %save;                 
   my %done;

   my $extraction_stmt = qq[SELECT v.variation_id,
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
                            AND    v.variation_id  = al.variation_id ];
  

   my $data_ext_sth = $var_dba->dbc->prepare($extraction_stmt);       
  
   $data_ext_sth->execute($first, $last)|| die "ERROR extracting allele info\n";
  
   my $data = $data_ext_sth->fetchall_arrayref();
   
   foreach my $l(@{$data}){  
    
       $save{$l->[0]}{name}   = $l->[1];
      
       if($flip->{$l->[0]}){
       ## update allele for flips
	   $l->[4]=~ /\(/   ?   $l->[4] = revcomp_tandem($l->[4]) :
	       defined $QUICK_COMP{$l->[4]}  ?  $l->[4] = $QUICK_COMP{$l->[4]} : 
	       reverse_comp(\$l->[4]);
	   
       }

       push @{$save{$l->[0]}{allele_data}}, [$l->[2], $l->[3], $l->[4], $l->[5], $l->[6], $l->[7], $l->[8], $l->[9] ];


  }

  return (\%save);  

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


  my $var_dba    = shift;
  my $allele_ids = shift;  ## allele_id look up
  my $fail_list  = shift;  ## hash of arrays of [variation_id, allele]'s by fail class
  my $fail_all   = shift;  ## array of [subsnp_id, sample_id] 

  my $fail_ins_sth   = $var_dba->dbc->prepare(qq[ insert into failed_allele_working
                                                  (allele_id, failed_description_id)
                                                  values (?,?)
                                                ]);       

  my %done;     ## save on old pk to drop out eariler
  my %done_new; ## save on new pk too

  ## deal with single allele fails (by variation and allele)
  foreach my $reason (keys %{$fail_list}){ 

      foreach my $var( @{$fail_list->{$reason}}  ){
      ## if the failure is a property of the reported allele, fail all entries with this allele

          if ($done{$reason}{$var->[0]}{$var->[1]}){ next; }  ## duplicates arise due to running on multiple variation_feature.allele_string's

	  $done{$reason}{$var->[0]}{$var->[1]} = 1;

	  unless ( defined $allele_ids->{va}{$var->[0]}{$var->[1]}->[0] ){
	      warn "Error finding allele ids to fail for var: $var->[0] & allele: $var->[1] - not failing on reason $reason\n";
	      next;
	  }
	  foreach my $al(@{$allele_ids->{va}{$var->[0]}{$var->[1]}}){
             $done_new{$al}{$reason} = 1; ## flag by reason and new id for screening later
	     $fail_ins_sth->execute($al, $reason)|| die "ERROR inserting allele fails info\n";
      }
    }
  }

  ### deal with submission sets of allele fails (by ss id and sample)
  foreach my $set(@{$fail_all}){
      unless (defined $allele_ids->{ss}{$set->[0]}{$set->[1]} ){
	  warn "Error finding allele ids to fail for var: $set->[0] & sample: $set->[1] - not failing \n";
	  next;
       }
      foreach my $allele_id (@{$allele_ids->{ss}{$set->[0]}{$set->[1]} }){
	  next if $done_new{$allele_id}{11} == 1;
	  $done_new{$allele_id}{11} = 1;
         $fail_ins_sth->execute($allele_id, 11 )|| die "ERROR inserting allele fails info\n";
      }
  }  

}
=head2 write_variation_features

  Create new version of variation feature tables with alleles ordered ref/alt where possible and complimented where necessary

=cut
sub write_variation_features{
    
  my $var_dba  = shift;
  my $to_load  = shift;
  my $evidence = shift;

  my $varfeat_ins_sth = $var_dba->dbc->prepare(qq[ insert  into variation_feature_working
                                                  (variation_id, variation_name, seq_region_id,  seq_region_start, 
                                                   seq_region_end, seq_region_strand,  allele_string, map_weight,  
                                                   source_id, consequence_types, variation_set_id, somatic,
                                                   class_attrib_id, alignment_quality, evidence)
                                                  values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
                                                ]);       
  
  
  foreach my $data (@{$to_load}){ 

    die "No allele string for $data->{name}\n" unless defined  $data->{allele};
    $varfeat_ins_sth->execute( $data->{v_id},
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
			       $evidence->{$data->{v_id}}
                             )|| die "ERROR importing variation feature info\n";


  }
}


=head2 write_variation

  Populate a variation_working table with new flipped & evidence statuses

=cut
sub write_variation{
    
  my $self     = shift;
  my $evidence = shift;
  my $flip     = shift;
        
  my $var_dba = $self->get_species_adaptor('variation');

  my $first = $self->required_param('start_id');
  my $last  = $first + $self->required_param('batch_size') -1;
  if($first ==1){$last--;}

  my %minor_allele;

  my $var_ext_sth = $var_dba->dbc->prepare(qq[ select variation_id,
                                                      source_id,
                                                      name, 
                                                      flipped,
                                                      class_attrib_id,
                                                      somatic,
                                                      minor_allele,
                                                      minor_allele_freq,
                                                      minor_allele_count,
                                                      clinical_significance_attrib_id
                                               from variation
                                               where variation_id between ? and ?
                                              ]);       


  my $var_ins_sth = $var_dba->dbc->prepare(qq[ insert into variation_working
                                               (variation_id,
                                               source_id,
                                               name, 
                                               flipped,
                                               class_attrib_id,
                                               somatic,
                                               minor_allele,
                                               minor_allele_freq,
                                               minor_allele_count,
                                               clinical_significance_attrib_id,
                                               evidence)
                                               values (?,?,?,?,?,?,?,?,?,?,?)
                                              ]); 

      
  

  $var_ext_sth->execute($first,$last)||die;
  my $data = $var_ext_sth->fetchall_arrayref();
  foreach my $v (@{$data} ){ 

    $flip->{$v->[0]} = 0 unless defined $flip->{$v->[0]}; ## update non-flipped with 0

    ## save minor allele for later checking
    $minor_allele{$v->[0]} = $v->[6] if defined $v->[6];

    ## substitute in new values as appropriate
    $var_ins_sth->execute( $v->[0],
                           $v->[1],
                           $v->[2],
                           $flip->{$v->[0]},
                           $v->[4],
                           $v->[5],
                           $v->[6],
                           $v->[7],
                           $v->[8],
                           $v->[9],
                           $evidence->{$v->[0]})||    die "ERROR writing variation table\n";    

  }
  return \%minor_allele;

}


=head2 write_allele

  Create new & old version of allele table with alleles complimented where neccessary
  Expects raw import format as input

=cut
sub write_allele{

  my $self      = shift;
  my $var_dba   = shift;
  my $var_data  = shift;

  my %save_id;

  my $code = read_allele_code($var_dba);

  my $allele_old_ins_sth = $var_dba->dbc->prepare(qq[ insert  into MTMP_allele_working
                                                    (allele_id, variation_id, subsnp_id, allele, frequency, sample_id, count)
                                                    values (?,?,?,?,?,?,?)
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
    
      $allele_new_ins_sth->execute( $var, $allele->[1], $code->{$allele->[2]}, $allele->[3], $allele->[4], $allele->[5], $allele->[6]) || die "ERROR inserting allele info\n";
      my $allele_id = $var_dba->db_handle->last_insert_id(undef, undef, qw(allele_working allele_id))|| die "no insert id for allele\n";


      ## save allele id for later fail update
      push @{$save_id{va}{$var}{$allele->[2]}},         $allele_id;
      push @{$save_id{ss}{$allele->[1]}{$allele->[4]}}, $allele_id;

        ## keep allele data in un-coded format for Mart for non-human databases
       unless($self->required_param('species') =~/Homo|Human|musculus|mouse/i){
          $allele_old_ins_sth->execute($allele_id, $var, $allele->[1], $allele->[2], $allele->[3], $allele->[4], $allele->[5]) || die "ERROR inserting allele info\n";
     }

    }
  }
  return \%save_id;
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


sub find_failed_variation_set_id{

    my $var_dba = shift;

    my $variation_set_ext_sth  = $var_dba->dbc->prepare(qq[ select variation_set_id
                                                            from variation_set
                                                            where name = ?
                                                           ]);

    ## check if present
    $variation_set_ext_sth->execute('All failed variations')  || die "Failed to extract failed variant set id\n";
    my $failed_set_id = $variation_set_ext_sth->fetchall_arrayref();

    die "Failed set not available\n" unless defined $failed_set_id->[0]->[0];

    return $failed_set_id->[0]->[0];
}
=head check_minor_allele

Description: Compare minor allele string to allele string for any variants with both and save suspect variants
Status     : At Risk

=cut
sub check_minor_allele{

    my $self          = shift;
    my $minor_allele  = shift; 
    my $allele_string = shift;
    
    my $mafail_ins_sth  = $var_dba->dbc->prepare(qq[ insert into failed_minor_allele_tmp (variation_id) values (?) ]);

    foreach my $var_id (keys %{$minor_allele}){

      next unless $allele_string->{$var_id};
      my @expected = split/\//, $allele_string->{$var_id};

      my $is_found = 0;
      foreach my $expected (@expected){
	  $is_found =1 if $minor_allele->{$var_id} eq $expected;
      } 
      unless ($is_found ===1){
	  $minor_allele_fail_ins_sth->execute($var_id});
    }

}
  

### put somewhere sensible
sub unique {
  my %a;
  map { $a{$_} = 1; } @_;
  return sort keys %a;
}


1;
