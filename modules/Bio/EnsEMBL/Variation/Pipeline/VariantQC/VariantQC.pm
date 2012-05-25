
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



=head2 run_allele_checks

  Run checks on variant data; submit failed_variation and strand-corrected variation feature database updates 

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
	
	## Type 1 - fail variant if  >3 mapping seen
	
 	if($var->{map} >3){
	   
	    push @{$fail_variant{1}}, $var->{v_id};
	}

	## Type 4 - novariation fails - flag variants & alleles as fails & don't run further checks

	if($var->{allele} =~ /NOVARIATION/){
	    push @{$fail_variant{4}}, $var->{v_id};
	    push @{$fail_allele{4}}, [$var->{v_id}, "NOVARIATION"];   ## unflipped allele failed throughout
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


	## Type 12  flag variation as fail if it has [A/T/G/C] allele string 
	my $all_possible_check = check_four_bases($expanded);
	if ($all_possible_check ==1){	    
	    push @{$fail_variant{12}}, $var->{v_id}
	}
	## Type 14 resolve ambiguities before reference check - flag variants & alleles as fails

	if($expanded =~ /$AMBIGUITIES/){
	    $expanded = remove_ambiguous_alleles(\$expanded,$AMBIGUITIES);
	    push @{$fail_variant{14}}, $var->{v_id};

	    ## identify specific alleles to fail
	    my $ambiguous_alleles = find_ambiguous_alleles($var->{allele},$AMBIGUITIES);
	    foreach my $amb_al( @{$ambiguous_alleles} ){
		push @{$fail_allele{14}},  [$var->{v_id}, $amb_al]; ## unflipped allele failed throughout
	    }

	}
	## Further checks only run for variants with single map locations
	
	next unless  $var->{map} == 1;

	
	## flip allele string if on reverse strand 
	## - store ids to update in allele table later
	
	if($var->{strand} eq "-1" ){
	    #warn "flipping: $var->{allele} strand: $var->{strand} mapwt: $var->{map}\n";
	    reverse_comp(\$var->{allele} );
	    $var->{strand} = 1;
	    ### store to flip everything else
	    $flip{$var->{v_id}} = 1;
	}
	
	
	
	# Extract reference sequence to run ref checks [ compliments for reverse strand multi-mappers]
	
	my $ref_seq = get_reference_base($var, $slice_ad) ;
	

	
	my @alleles = split/\//, $var->{allele};
	foreach my $al(@alleles){
	    if($al eq $ref_seq){
		$var->{ref} = $al;
	    }
	    else{$var->{alt} .= "/" . $al;}
	}
	
	if( defined $var->{ref}){
	    ## re-order the allele_string field with the reference allele first 
	    $var->{allele} = $var->{ref} . $var->{alt};
	    
	    ##  Type 15 - flag variants as fails if the reference allele does not match the expected length from the coordinates 
	    push @{$fail_variant{15}}, $var->{v_id} if check_variant_size($self, $var);

	}
	else{
	    ##  Type 2 - flag variants as fails if neither allele matches the reference
	    push @{$fail_variant{2}}, $var->{v_id} ;
	}	
    }
    

    ## Database updates 
    
    ## write to new variation featues 
    insert_variation_features($var_dba, $to_check);
 
   
    ## update failed_variation table
    write_variant_fails($var_dba, \%fail_variant);
   
 
     
    ## allele-specific checks - run on updated variation_feature table
    run_allele_checks($self, \%flip, \%fail_allele);  
  
}



=head2 run_allele_checks

  Run checks on allele data; submit submit failed_allele and strand-corrected allele database updates 

=cut
sub run_allele_checks {

    my $self   = shift;
    my $flip   = shift; 
    my $fail   = shift;

    ### start and end variation_id supplied
    my $first = $self->required_param('start_id'); 
    my $last  = $first + $self->required_param('batch_size') - 1;
    if($first ==1){$last--;} 

    open my $temp_log, ">", "al_log_$first";
    my %fail_all;


    my $var_dba = $self->get_species_adaptor('variation');


    ## export current allele & variation_feature data - alleles flipped where needed
    my $var_data = export_allele_data($var_dba, $first, $last, $flip);
      

    foreach my $var( keys %$var_data){

	## cope with unmapped variants with alleles or named alleles, but don't check them
	next if $var_data->{$var}->{vf_allele} =~ /NULL/;  

	## isolate expected alleles 
	my %expected_alleles;
	my @expected_alleles = split/\//, $var_data->{$var}->{vf_allele};
	foreach my $exp(@expected_alleles){ $expected_alleles{$exp} = 1;}


	## check through allele submissions
	foreach my $submitted_data (@{$var_data->{$var}->{allele_data}}){	   

	    unless( $expected_alleles{$submitted_data->[6]}){ ## check expected on allele not allele_code

		if(defined $submitted_data->[4]){
		    ###fail whole experimental set
		    $fail_all{$submitted_data->[1]}{$submitted_data->[4]} = 1;   ## subsnp_id, sample_id
		    print $temp_log "$submitted_data->[6] not in ". $var_data->{$var}->{vf_allele} . " Failing exp ss:$submitted_data->[1] sam:$submitted_data->[4]\n";
		}
		else{
		    push @{$fail->{11}}, [$submitted_data->[0], $submitted_data->[5] ];          ## var_id, allele
		    print $temp_log "$submitted_data->[6] not in ". $var_data->{$var}->{vf_allele} . " Failing one  var: $submitted_data->[0] allele: $submitted_data->[5]\n";
		}
	    }
	}
    }
    ## write allele records to new table flipping strands as needed
    write_allele($var_dba, $var_data);

    ### fail full experimental result for sample & subsnp_id
    write_allele_fails($var_dba , $fail, \%fail_all, $var_data);
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
                                                         vf.consequence_type,
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
	$save{consequence_type} = $l->[10];
	$save{variation_set_id} = $l->[11];
	$save{somatic}          = $l->[12];
	$save{class_attrib_id}  = $l->[13];
	$save{seqreg_name}      = $l->[14];

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
                                                         vf.consequence_type,
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
	 $save{consequence_type} = $l->[9];
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

=cut
sub export_allele_data{

     my ($var_dba, $first, $last, $flip) = @_;
     
     my %save;

     
     # hack to use with old schema
    # my %allele_codes;  ## save as looked up 

     #my $allele_code_ext_sth  = $var_dba->dbc->prepare(qq[ select allele_code_id from allele_code  where allele = ? ]);	   
     
     #my $allele_code_ins_sth  = $var_dba->dbc->prepare(qq[ insert into allele_code (allele) values (?) ]);
         
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
                                                         al.allele       
                                                FROM  variation v join allele al on(v.variation_id = al.variation_id )
                                                    left outer join variation_feature_working vf on (vf.variation_id = v.variation_id )
                                                WHERE
                                                    v.variation_id between ? and ?
                                                  ]);       
    
     $data_ext_sth->execute($first, $last)|| die "ERROR extracting variation feature info\n";
    
     my $data = $data_ext_sth->fetchall_arrayref();
     
     foreach my $l(@{$data}){
	 ## distinct on query causes database tmp to fill - handling hackily instead
	 next if $done{"$l->[0]\_$l->[2]\_$l->[3]"} ;
	 $done{"$l->[0]\_$l->[2]\_$l->[3]"} = 1;
	 $save{$l->[0]}{name}         = $l->[1];
	 $save{$l->[0]}{vf_allele}    = $l->[2];
		
	if($flip->{$l->[0]}){

	    ## update allele for flips
	    reverse_comp(\$l->[5]);
	    reverse_comp(\$l->[9]);
	   # unless ($allele_codes{$l->[9]}){
  	   #	$allele_codes{$l->[9]} = get_allele_code($l->[9], $allele_code_ext_sth, $allele_code_ins_sth);
	    #}
	    ## updated allele code for flip
	    #$l->[5] = $allele_codes{$l->[9]} ;
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
    
    if($var->{end} < $var->{start}){ ## convention for insertions to reference
	$ref_seq = "-";
    }
    else{
	
	# retrieve the reference sequence at that mapping 
	#warn "Getting slice for $var->{name} $var->{seqreg_name}, $var->{start}, $var->{end}\n";
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

sub remove_ambiguous_allele{

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

=head2 check_variant_size

  Compares reference allele string to coordinates for variation and checks length appropriate
  Returntype : 1 is OK, 0 if failed

=cut

sub check_variant_size{

    my $self     = shift;
    my $var_data = shift;
    
    my $ref_length = $var_data->{end} - $var_data->{start} +1;


    ### insertion to reference
    if( $var_data->{ref} eq "-" && $ref_length  == 0){ return 0;}

    ### deletion of reference or substitution- coordinates should reflect length of deleted string
    elsif( $ref_length  == length($var_data->{ref})){ return 0;}


    else{
	warn "Check here - $var_data->{name} inconsistant coords $var_data->{ref}  $var_data->{start} ->$var_data->{end}  \n";
	return 1;
    }
}

=head2 check_variant_size

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




sub write_allele_fails{


    my $var_dba   = shift; 
    my $fail_list = shift;
    my $fail_all  = shift;  
    my $var_data  = shift;

    ## find allele id to fail
   # my $allele_id_ext_sth  = $var_dba->dbc->prepare(qq[select allele_id 
   #                                                   from allele, allele_code
   #                                                   where allele_code.allele =?
   #                                                   and allele_code.allele_code_id = allele.allele_code_id
   #                                                   and allele.variation_id =?
   #                                                   ]);
    my $allele_id_ext_sth  = $var_dba->dbc->prepare(qq[select allele_id 
                                                      from allele
                                                      where allele =?
                                                      and allele.variation_id =?
                                                      ]);


    my $fail_ins_sth   = $var_dba->dbc->prepare(qq[insert into failed_allele_working
                                                  (allele_id, failed_description_id)
                                                  values (?,?)
                                                  ]);       
            
    my %done;
   
    foreach my $reason (keys %{$fail_list}){ 

	foreach my $var( @{$fail_list->{$reason}}  ){

	    if ($done{$reason}{$var->[0]}{$var->[1]}){ next; }  ## duplicates arise due to running on multiple variation_feature.allele_string's

	    $done{$reason}{$var->[0]}{$var->[1]} = 1;
	    
	    $allele_id_ext_sth->execute($var->[1], $var->[0])|| die "ERROR extracting allele id info\n";
	    my $allele_id = $allele_id_ext_sth->fetchall_arrayref();

	    if ($allele_id->[0]->[0]){
		
		$fail_ins_sth->execute($allele_id->[0]->[0], $reason)|| die "ERROR inserting allele fails info\n";
	    }
	    else{
		warn "Error finding allele id for var :  $var->[0] & allele  $var->[1] - not failing \n";
	    }	    

	}


	foreach my $var( keys %$var_data){
	    
	    ## check through allele submissions & fail set of results for same subsnp & sample
	    foreach my $submitted_data (@{$var_data->{$var}->{allele_data}}){
		
		next if $done{$submitted_data->[0]};       ## check for duplicates with 2 badly formed alleles in the set
		next unless defined $submitted_data->[4];  ## only fails with sample info need be checked
		
		if( defined  $$fail_all{$submitted_data->[1]}{$submitted_data->[4]}  ){ ## failed by ssid & sample
		    
		    $done{$submitted_data->[0]} = 1; ## old pk

		    ## look up new pk on base and variation_id
		    $allele_id_ext_sth->execute($submitted_data->[5] , $var)|| die "ERROR extracting allele id info\n";
		    my $allele_id = $allele_id_ext_sth->fetchall_arrayref();
		    
		    $fail_ins_sth->execute($allele_id->[0]->[0], 11 )|| die "ERROR inserting allele fails info\n";
		}
	    }
	}
	
    }
}

sub insert_variation_features{
    
    my $var_dba  = shift;
    my $to_check = shift;

    my $varfeat_ins_sth = $var_dba->dbc->prepare(qq[insert into variation_feature_working
                                                  (variation_id, variation_name, seq_region_id,  seq_region_start, seq_region_end, seq_region_strand,  
                                                   allele_string, map_weight,  source_id, consequence_type, variation_set_id, somatic, class_attrib_id)
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
				  $data->{consequence_type},
				  $data->{variation_set_id},
				  $data->{somatic},
				  $data->{class_attrib_id})|| die "ERROR importing variation feature info\n";


    }
}



sub write_variant_flips{
    
    my $var_dba = shift;
    my $flip    = shift;
    
    my $flip_ins_sth = $var_dba->dbc->prepare(qq[insert into variation_to_reverse_working
                                                 (variation_id)
                                                 values (?)
                                                 ]);       
    

    foreach my $var (keys %{$flip} ){ 

	$flip_ins_sth->execute($var)|| die "ERROR inserting variation flip info\n";	
    }
}


sub write_allele{
    ### write data to new allele table flipping where needed; best way to recover partial updates caused by fails
    
    my ($var_dba, $var_data) = @_;
    my %done;

    
    my $allele_ins_sth     = $var_dba->dbc->prepare(qq[ insert into allele_working
                                                       (variation_id, subsnp_id, allele, frequency, sample_id, count)
	                                                values (?, ?, ?, ?, ?,? )
                                                      ]);
  
    foreach my $var (keys %$var_data){

	foreach my $allele (@{$var_data->{$var}->{allele_data}}){  ## array of data from 1 row in table
	    next if $done{$allele->[0]} ==1;  ## filtering on old pk
	    $done{$allele->[0]} = 1;

	    $allele_ins_sth->execute( $var, $allele->[1], $allele->[2], $allele->[3], $allele->[4], $allele->[5]) || die "ERROR inserting allele info\n";
	}
    }
}


### put somewhere sensible
sub unique {
    my %a;
    map { $a{$_} = 1; } @_;
    return sort keys %a;
}


1;
