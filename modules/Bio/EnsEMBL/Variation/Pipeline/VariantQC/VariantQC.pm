=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Variation::Pipeline::VariantQC::VariantQC

=head1 DESCRIPTION

Runs basic quality control on variant and allele data imported from extrenal sources

=cut

package Bio::EnsEMBL::Variation::Pipeline::VariantQC::VariantQC;


use strict;
use warnings;
use Data::Dumper ;
use Bio::DB::Fasta;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp expand);
use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);
use Bio::EnsEMBL::Variation::Utils::dbSNP qw(get_alleles_from_pattern);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(revcomp_tandem);
use Bio::EnsEMBL::Variation::Utils::QCUtils qw( check_four_bases get_reference_base check_illegal_characters check_for_ambiguous_alleles remove_ambiguous_alleles find_ambiguous_alleles summarise_evidence count_rows count_group_by check_variant_size);


our $DEBUG   = 0;


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
   
  ## variation_feature specific checks + return bad allele
  my ($var_data, $flip, $allele_string, $failed_variant, $failed_varallele) = $self->run_variation_checks();

  ## allele specific checks (return by allele id & ss/sample)
  my ($allele_data, $allele_fails, $failed_ss) = $self->run_allele_checks( $flip, $failed_varallele, $allele_string);  
 
  my $var_dba = $self->get_species_adaptor('variation');

  my $first = $self->required_param('start_id');
  my $last  = $first + $self->required_param('batch_size') -1;
  if($first ==1){$last--;}

  ## find supporting evidence & merge to single string  
  my $evidence_attribs ;
  if($self->required_param('evidence_check') ==1){
      
       $evidence_attribs = summarise_evidence( $var_dba->dbc(),
                                               $self->required_param('species'),
                                               $first,
                                               $last);   

  }


  ## Database updates 

  ## write to new failed_variation table
  my $failed_var_ids = write_variant_fails($var_dba, $failed_variant);

  ## write to new variation featues 
  write_variation_features($var_dba, $var_data, $evidence_attribs, $failed_var_ids ); 
 
  ## write to new allele table & look up additional allele ids for failed experiments
  my $final_allele_fails = $self->write_allele($var_dba, $allele_data, $allele_fails, $failed_ss );

  ### fail full experimental result for sample & subsnp_id
  write_allele_fails($var_dba , $final_allele_fails );


  ## write to new variation with evidence string & flip status
  my $minor_allele = $self->write_variation(  $evidence_attribs , $flip, $failed_var_ids );
  

  ## check dbSNP MAF against variation_feature allele string
  ## warn, but don't fix automatically (1KG genotypes loaded later)
  if($self->required_param('species') =~/homo/i){
     $self->check_minor_allele($var_dba, $minor_allele,  $allele_string); 
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
  my %fail_varallele; # hash of arrays of arrays failed variation_ids & alleles
  my %flip;           # hash of variation ids which have their strand changed in this process
  my %allele_string;  # hash of expected allele for each variation_id strings saved for allele checking

  


  my $var_dba = $self->get_species_adaptor('variation');

  ## slice needed for ref check
  my $core_dba = $self->get_species_adaptor('core');
  my $slice_ad = $core_dba->get_SliceAdaptor;


  ## export current variation_feature data
  my ($to_check, $map_count, $strand_summary) = export_data_adding_allele_string($var_dba, $first, $last);


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

### Add post db patch
    ## Type 20 - It is unlikely a variant at the first base of a reference sequence is a good call
#    if( $var->{start}  ==1){
#      push @{$fail_variant{20}}, $var->{v_id};
#      $var->{variation_set_id} = $failed_set_id ;
#    }


    ## Type 4 - novariation fails - flag variants as fails & don't run further checks

    if($var->{allele} =~ /NOVARIATION/){
      push @{$fail_variant{4}}, $var->{v_id};
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
         push @{$fail_varallele{ $var->{v_id} }{ $ill_al }},13;   ## unflipped allele failed throughout
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
          push @{$fail_varallele{ $var->{v_id} }{ $amb_al }},14; 
	}
    }

    ## Compliment allele string only for variants with 1 genomic location or same strand for all locations 
    next if $map_count->{$var->{v_id}} > 1  && $strand_summary->{$var->{v_id}} eq '0';


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

      ### store variation_id to use when flipping alleles & allele stringfor multi-map flips
      $flip{$var->{v_id}} = 1;
      $allele_string{$var->{v_id}} = $var->{allele};
    }
  
    ## Reference match checks and re-ordering only run for variants with 1 genomic location
    next if $var->{map} > 1;
  
    # Extract reference sequence to run ref checks [ compliments for reverse strand multi-mappers]
  
    my $ref_seq ;
    if( $self->param('use_seqdb')  == 1){
        my $fasta_file = $self->required_param('pipeline_dir') . "/genomic.fa";
        my $db = Bio::DB::Fasta->new($fasta_file );
        $ref_seq = get_reference_base($var, $db, "fasta_seq") ;
    }
    else{
        $ref_seq = get_reference_base($var, $slice_ad, "coredb") ;
    }
    unless(defined $ref_seq){ 
      ## don't check further if obvious coordinate error
      push @{$fail_variant{15}}, $var->{v_id};
      $var->{variation_set_id} = $failed_set_id ;
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
       ## is the length ok even if the base doesn't match?
       my $size_ok = check_variant_size( $var->{start} , $var->{end} , $al);
       $match_coord_length = 1 if $size_ok ==1;
       ## save as alt allele if it doesn't match the reference
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
 

    return ($to_check, \%flip, \%allele_string, \%fail_variant, \%fail_varallele) ;

}


=head2 run_allele_checks

  Check alleles in allele table match alleles in variation_feature table
      - catches frequency data with incorrect strand or no-calls counted in frequency calculations
      - take in failed variation_id & allele string combination to look up allele_id 

  Input: flip status, failed allele by var and strand-corrected allele updates 
  Output: alleles, failed ids, failed ss/population combinations

=cut
sub run_allele_checks {

   my $self             = shift;
   my $flip             = shift; 
   my $failed_varallele = shift; ## hash of varid allele ->[array of fail codes] based on ref string checks
   my $allele_string    = shift;

   my %fail;     ## hash of fail reasons => array of allele ids
   my %fail_all; ## hash of var & ssids to fail


   ## export current allele & variation_feature data - alleles flipped where needed
   my $var_data = $self->export_allele_data( $flip );
   
   foreach my $var( keys %$var_data){
        
     my %expected_alleles;
     if (defined $allele_string->{$var} ){
         my @expected_alleles = split/\//, $allele_string->{$var};
         foreach my $exp(@expected_alleles){ 
           $expected_alleles{$exp} = 1;
         }
     }

       ## Apply QC checks
       foreach my $submitted_data (@{$var_data->{$var}->{allele_data}}){      
          ## $submitted_data content: [ al.allele_id, al.subsnp_id, al.allele,  al.frequency, al.population_id, al.count, al.frequency_submitter_handle ]

	   if ($submitted_data->[2] =~ /NOVARIATION/){
	       push @{$fail{4}}, $submitted_data->[0] ;
	       next;  ## no further checking
	   }
	   if( defined $failed_varallele->{ $var }{ $submitted_data->[2] }->[0] ){
	       foreach my $reason (@{ $failed_varallele->{ $var }{$submitted_data->[2] }}){
		   push @{$fail{$reason}}, $submitted_data->[0] ;
	       }
	       
	       next;  ## no further checking - must match ref
	   }

	## only check if strand corrected (mapped) reference variation available
        next unless defined $allele_string->{$var};
     

	   unless( exists $expected_alleles{$submitted_data->[2]} && $expected_alleles{$submitted_data->[2]} ==1 ){ 

	   ## checking ss id here because of the way strains are handled in import
	   if(defined $submitted_data->[4] && defined $submitted_data->[1] && $submitted_data->[1] >0  ){
           ###fail whole experimental set if sample info supplied
           $fail_all{$submitted_data->[1]}{$submitted_data->[4]} = 1; ## subsnp_id, sample_id
         }
         else{
           ##Additional submitted allele data from dbSNP does not agree with the dbSNP refSNP alleles (on allele_id)
           push @{$fail{11}}, $submitted_data->[0] ;       
         }
       }
     }
   }

   return ($var_data, \%fail, \%fail_all );



}


=head2 export_data_adding_allele_string

  Extract variant data for checking taking allele_string from allele_string table 
  if required (dbSNP import pipeline does not populate variation_feature.allele_string)

=cut
sub export_data_adding_allele_string{

  my ($var_dba, $first, $last) = @_;
  
  my @to_check;
  my %strand_summary;  ## hold strands seen for multi-mapping variants here rather than check later
  my %map_count;       ## count all locations including non-ref locations not used in map_weight calculation

  my %found_minor;     ## dbSNP flag 2 alleles as minor for some tri-allelic variants

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
                                                      v.class_attrib_id,
                                                      sr.name,
                                                      vf.alignment_quality,
                                                      als.allele_string,
                                                      vf.validation_status,
                                                      maf.allele,
                                                      maf.freq,
                                                      maf.count,
                                                      maf.is_minor_allele,
                                                      vf.flags
                                                 FROM variation_feature vf
                                                      left outer join tmp_map_weight_working  tmw on (vf.variation_id = tmw.variation_id ),
                                                      seq_region sr,
                                                      allele_string als,
                                                      variation v
                                                      left outer join maf on ( maf.snp_id = v.snp_id)
                                                 WHERE  v.variation_id between ? and ? 
                                                 AND vf.seq_region_id = sr.seq_region_id
                                                 AND vf.variation_id = als.variation_id  
                                                 AND vf.variation_id = v.variation_id
                                                 ]);       
 

 
  
  $variant_ext_sth->execute($first, $last)|| die "ERROR extracting variation feature info\n";
  
  my $variant_data = $variant_ext_sth->fetchall_arrayref();

  my %potentially_no_minor; ## if bi-allelic & equal frequencies assign one as minor at random

  foreach my $l(@{$variant_data}){
 
      next if defined $l->[20] &&  $l->[20] eq "0" && 
	  (! $potentially_no_minor{$l->[2]} && $l->[18] ==0.5);

      ## some variants are assigned 2 minor alleles
      next if defined $found_minor{$l->[2]} && $found_minor{$l->[2]} ==1;
      $found_minor{$l->[2]} =1; 

      $potentially_no_minor{$l->[2]} = 1 if defined $l->[18] &&  $l->[18] ==0.5;

      $map_count{$l->[0]}++;
      my %save;

      $save{v_id}           = $l->[0];
      $save{name}           = $l->[1];
      $save{vf_id}          = $l->[2];
      $save{seqreg_id}      = $l->[3];
      $save{start}          = $l->[4];
      $save{end}            = $l->[5];
      $save{strand}         = $l->[6];
      

      $save{map}              = $l->[7]  || 0 ;   ## variants mapping only to haplotypes/patches have map_weight 0
      $save{source_id}        = $l->[8];
      $save{consequence_types}= $l->[9];
      $save{variation_set_id} = $l->[10];
      $save{somatic}          = $l->[11];
      $save{class_attrib_id}  = $l->[12];
      $save{seqreg_name}      = $l->[13];
      $save{align_qual}       = $l->[14];
      $save{validation_status}= $l->[16];

      $save{min_allele}    = $l->[17];
      $save{min_af}        = $l->[18];
      $save{min_al_count}  = $l->[19];
      $save{flags}         = $l->[21];

    if($l->[15] =~ /^(\(.*\))\d+\/\d+/){## handle tandem
      my $expanded_alleles = get_alleles_from_pattern($l->[15]); 
      $save{allele} = join"/",@{$expanded_alleles};
    }
    else{
      $save{allele}         = $l->[15];
    }

    push @to_check,\%save;

    ## summarise strands seen - if multi-mapping but always reverse, safe to flip
    if(defined $strand_summary{$l->[0]} && $strand_summary{$l->[0]}  ne $save{strand} ){
	$strand_summary{$l->[0]} = 0;
    }
    else{
	$strand_summary{$l->[0]} = $save{strand};
    }
  }

  return (\@to_check, \%map_count, \%strand_summary);  

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
                                    al.population_id,
                                    al.count,
                                    al.frequency_submitter_handle,
                                    als.allele_string
                            FROM   variation v, allele al, allele_string als
                            WHERE  v.variation_id between ? and ?
                            AND    v.variation_id  = al.variation_id 
                            AND    als.variation_id = v.variation_id ];

  
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
    Return hash of failed ids for display setting

=cut
sub write_variant_fails{


  my $var_dba   = shift; 
  my $fail_list = shift;

  my %failed_ids;

  my $fail_ins_sth = $var_dba->dbc->prepare(qq[insert ignore into failed_variation_working
                                               (variation_id, failed_description_id)
                                               values (?,?)
                                               ]);       
          

  foreach my $reason (keys %{$fail_list}){ 

    ## duplicates arise due to running on variation_features not variations
    my @fails = unique(@{$fail_list->{$reason}});    

    foreach my $var( @fails  ){
      $fail_ins_sth->execute($var, $reason)|| die "ERROR inserting variation fails info\n";
      $failed_ids{$var} = 1;
    }
  }
  return \%failed_ids;
}



=head2 write_allele_fails

    Update failed_allele_working with all variation_id/reason combinations


=cut
sub write_allele_fails{


  my $var_dba       = shift;
  my $allele_fails  = shift;  ## hash of arrays of allele_id by fail class


  my $fail_ins_sth   = $var_dba->dbc->prepare(qq[ insert ignore into failed_allele_working
                                                  (allele_id, failed_description_id)
                                                  values (?,?)
                                                ]);       


  foreach my $reason (keys %{$allele_fails}){ 

      foreach my $allele_id( @{ $allele_fails->{$reason} }  ){
      
	  $fail_ins_sth->execute($allele_id, $reason)|| die "ERROR inserting allele fails info\n";
      }
    }
}
=head2 write_variation_features

  Create new version of variation feature tables with alleles ordered ref/alt where possible and complimented where necessary

=cut
sub write_variation_features{
    
  my $var_dba         = shift;
  my $to_load         = shift;
  my $evidence        = shift;
  my $failed_var_ids  = shift;

  my $varfeat_ins_sth = $var_dba->dbc->prepare(qq[ insert ignore into variation_feature_working
                                                  ( variation_feature_id, variation_id, variation_name, seq_region_id,  seq_region_start, 
                                                   seq_region_end, seq_region_strand,  allele_string, map_weight,  
                                                   source_id, consequence_types, variation_set_id, somatic,
                                                   class_attrib_id, alignment_quality, flags, evidence_attribs,
                                                   minor_allele, minor_allele_freq, minor_allele_count,validation_status, display )
                                                  values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
                                                ]);       
  
  
  foreach my $data (@{$to_load}){ 

    my $evidence_attribs = join(",", @{$evidence->{$data->{v_id}}} ) if defined $evidence->{$data->{v_id}};

    my $display = 1;
    $display = 0 if defined $failed_var_ids->{$data->{v_id}}; 

    die "No allele string for $data->{name}\n" unless defined  $data->{allele};
    $varfeat_ins_sth->execute( $data->{vf_id},
                               $data->{v_id},
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
                               $data->{flags},
			       $evidence_attribs,
                               $data->{min_allele},
                               $data->{min_af},
                               $data->{min_al_count},
			       $data->{validation_status},
                               $display
                             )|| die "ERROR importing variation feature info\n";


  }
}


=head2 write_variation

  Populate a variation_working table with new flipped & evidence statuses

=cut
sub write_variation{
    
  my $self           = shift;
  my $evidence       = shift;
  my $flip           = shift;
  my $failed_var_ids = shift;
        
  my $var_dba = $self->get_species_adaptor('variation');

  my $first = $self->required_param('start_id');
  my $last  = $first + $self->required_param('batch_size') -1;
  if($first ==1){$last--;}

  my $var_ext_sth = $var_dba->dbc->prepare(qq[ select v.variation_id,
                                                      v.source_id,
                                                      v.name, 
                                                      v.flipped,
                                                      v.class_attrib_id,
                                                      v.somatic,
                                                      maf.allele,
                                                      maf.freq,
                                                      maf.count,
                                                      maf.is_minor_allele,
                                                      clinical_significance,
                                                      validation_status
                                               from variation v left outer join maf on ( maf.snp_id = v.snp_id)
                                               where v.variation_id between ? and ?

                                              ]);       


  my $var_ins_sth = $var_dba->dbc->prepare(qq[ insert ignore into variation_working
                                               (variation_id,
                                               source_id,
                                               name, 
                                               flipped,
                                               class_attrib_id,
                                               somatic,
                                               minor_allele,
                                               minor_allele_freq,
                                               minor_allele_count,
                                               clinical_significance,
                                               evidence_attribs,
                                               validation_status,
                                               display)
                                               values (?,?,?,?,?,?,?,?,?,?,?,?,?)
                                              ]); 

      
  my %minor_allele;
  my %potentially_no_minor; ## if bi-allelic & equal frequencies assign one as minor at random

  $var_ext_sth->execute($first,$last)||die;
  my $data = $var_ext_sth->fetchall_arrayref();
  foreach my $v (@{$data} ){ 

    

    next if defined $v->[9] &&  $->[9] eq "0" && 
	  (! $potentially_no_minor{$v->[0]} && $v->[7] !=0.5);

    ## some variants are assigned 2 minor alleles
    next if defined $minor_allele{$v->[0]} ;

    $potentially_no_minor{$v->[0]} = 1 if defined  $v->[7] && $v->[7] ==0.5;
     
    $flip->{$v->[0]} = 0 unless defined $flip->{$v->[0]} ; ## update non-flipped with 0


    my $display = 1;
    $display = 0 if defined $failed_var_ids->{$v->[0]};


    ## save minor allele for later checking
    $minor_allele{$v->[0]} = $v->[6] if defined $v->[6];

    my $evidence_attribs = join(",", @{$evidence->{$v->[0]}}) if defined $evidence->{$v->[0]};

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
                           $v->[10],
                           $evidence_attribs,
                           $v->[11],
                           $display )||    die "ERROR writing variation table\n";    

  }
  return \%minor_allele;

}


=head2 write_allele

  Create new & old version of allele table with alleles complimented where neccessary
  Expects raw import format as input

=cut
sub write_allele{

  my $self         = shift;
  my $var_dba      = shift;
  my $var_data     = shift;
  my $allele_fails = shift; ## individual alleles to fail
  my $fail_exp     = shift; ## ss ids and samples to fail if one is suspect



  my $code = read_allele_code($var_dba);
 
  my $allele_new_ins_sth = $var_dba->dbc->prepare(qq[ insert ignore into allele_working
                                                     (allele_id, variation_id, subsnp_id, allele_code_id, frequency, population_id, count, frequency_submitter_handle)
                                                     values (?, ?, ?, ?, ?, ?, ?, ? )
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
    
      $allele_new_ins_sth->execute(  $allele->[0], $var, $allele->[1], $code->{$allele->[2]}, $allele->[3], $allele->[4], $allele->[5], $allele->[6]) || die "ERROR inserting allele info\n";

      push @{$allele_fails->{11}}, $allele->[0] if defined $allele->[4] &&
         defined  $fail_exp->{ $allele->[1]}->{$allele->[4]};

    }
  }
  ## return fail list with added data
  return $allele_fails;
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

    my $var_dba       = shift;
    my $minor_allele  = shift; 
    my $allele_string = shift;
    
    my $minor_allele_fail_ins_sth  = $var_dba->dbc->prepare(qq[ insert ignore into tmp_failed_minor_allele (variation_id) values (?) ]);

    foreach my $var_id (keys %{$minor_allele}){

	next unless $allele_string->{$var_id};
	my @expected = split/\//, $allele_string->{$var_id};
	
	my $is_found = 0;

	foreach my $expected (@expected){
	    $is_found = 1 if $minor_allele->{$var_id} eq $expected;
	} 
	
	$minor_allele_fail_ins_sth->execute($var_id) unless $is_found == 1;

    }
}
  

### put somewhere sensible
sub unique {
  my %a;
  map { $a{$_} = 1; } @_;
  return sort keys %a;
}


1;
