#
# Ensembl module for Bio::EnsEMBL::Variation::TranscriptVariation
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::TranscriptVariation

=head1 SYNOPSIS

  use Bio::EnsEMBL::Variation::TranscriptVariation;

  $tr_var = Bio::EnsEMBL::Variation::TranscriptVariation->new
    (-transcript        => $transcript,
     -pep_allele_string => 'N/K',
     -cdna_start        => 1127,
     -cdna_end          => 1127,
     -cds_start         => 558,
     -cds_end           => 558,
     -translation_start => 318,
     -translation_end   => 318,
     -consequence_type  => 'NON_SYNONYMOUS_CODING');


  print "variation: ", $tr_var->variation_feature()->variation_name(), "\n";
  print "transcript: ", $tr_var->transcript->stable_id(), "\n";
  print "consequence type: ", (join ",", @{$tr_var->consequence_type()}), "\n";
  print "cdna coords: ", $tr_var->cdna_start(), '-', $tr_var->cdna_end(), "\n";
  print "cds coords: ", $tr_var->cds_start(), '-', $tr_var->cds_end(), "\n";
  print "pep coords: ", $tr_var->translation_start(), '-',
        $tr_var->translation_end(), "\n";
  print "amino acid change: ", $tr_var->pep_allele_string(), "\n";


=head1 DESCRIPTION

A TranscriptVariation object represents a variation feature which is in close
proximity to an Ensembl transcript.  A TranscriptVariation object has several
attributes which define the relationship of the variation to the transcript.

=head1 AUTHOR - Graham McVicker

=head1 CONTACT

Post questions to the Ensembl development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::TranscriptVariation;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp expand);
use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Variation::ConsequenceType;

our @ISA = ('Bio::EnsEMBL::Storable');

my %CONSEQUENCE_TYPES = %Bio::EnsEMBL::Variation::ConsequenceType::CONSEQUENCE_TYPES;
my %AFFECTS_PEPTIDE = %Bio::EnsEMBL::Variation::ConsequenceType::AFFECTS_PEPTIDE;

=head2 new

  Arg [-ADAPTOR] :
    Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor

  Arg [-DBID] :
    int the unique internal identifier for this TranscriptVariation

  Arg [-TRANSCRIPT] :
    Bio::EnsEMBL::Transcript - The transcript affected by this
    TranscriptVariation

  Arg [-PEP_ALLELE_STRING] :
    string - A '/' delimited string representing amino acids altered by this
    variation

  Arg [-CDNA_START] :
    The start of this variation on the associated transcript in cdna
    coordinates

  Arg [-CDNA_END] :
    The end of this variation on the associated transcript in cdna coordinates

  Arg [-CDS_START] :
    The start of this variation on the associated transcript in cds
    coordinates

  Arg [-CDS_END] :
    The end of this variation on the associated transcript in cds coordinates

  Arg [-TRANSLATION_START] :
    The start of this variation on the translation of the associated transcript
    in peptide coordinates

  Arg [-TRANSLATION_END] :
    The end of this variation on the translation of the associated transcript
    in peptide coordinates

  Arg [-CONSEQUENCE_TYPE] :
    The type of this TranscriptVariation.  Must be one of:
    'INTRONIC', 'UPSTREAM', 'DOWNSTREAM', 'SYNONYMOUS_CODING',
    'NON_SYNONYMOUS_CODING', 'FRAMESHIFT_CODING', '5PRIME_UTR', '3PRIME_UTR'

  Example    : 
    $tr_var = Bio::EnsEMBL::Variation::TranscriptVariation->new
      (-transcript        => $transcript,
       -pep_allele_string => 'N/K',
       -cdna_start        => 1127,
       -cdna_end          => 1127,
       -cds_start         => 558,
       -cds_end           => 558,
       -translation_start => 318,
       -translation_end   => 318,
       -consequence_type  => 'NON_SYNONYMOUS_CODING');

  Description: Constructor. Instantiates a
               Bio::EnsEMBL::Variation::TranscriptVariation object
  Returntype : Bio::EnsEMBL::Variation::TranscriptVariation
  Exceptions : throw on bad argument
  Caller     : general, TranscriptVariationAdaptor
  Status     : At Risk

=cut

sub new {
  my $class = shift;

  my ($vf, $tr, $pep_allele, $cdna_start,$cdna_end, $cds_start, $cds_end, $tl_start,$tl_end, $consequence_type,
      $dbID, $adaptor, $transcript, $codons) =
    rearrange([qw(VARIATION_FEATURE TRANSCRIPT PEP_ALLELE_STRING CDNA_START
                  CDNA_END CDS_START CDS_END TRANSLATION_START TRANSLATION_END CONSEQUENCE_TYPE
                  DBID ADAPTOR TRANSCRIPT CODONS)], @_);

  if(defined($consequence_type)) {
      my @consequences = split /,/,@{$consequence_type};
      foreach my $consequence (@consequences){
	  $consequence = uc($consequence);  
	  if(!$CONSEQUENCE_TYPES{$consequence}) {
	      my $valid = join(',',map({"'$_'"} keys(%CONSEQUENCE_TYPES)));
	      throw("Type argument must be one of: $valid");
	  } 
      }   
  }

  if(defined($cdna_start) && ($cdna_start !~ /^\d+$/ || $cdna_start < 1)) {
    throw('CDNA start must be greater than or equal to 1');
  }

  if(defined($cdna_end) && ($cdna_end !~ /^\d+$/ || $cdna_start < 0)) {
    throw('CDNA end must be greater than or equal to 0');
  }

  if(defined($cds_start) && ($cds_start !~ /^\d+$/ || $cds_start < 1)) {
    throw('CDS start must be greater than or equal to 1');
  }

  if(defined($cds_end) && ($cds_end !~ /^\d+$/ || $cds_start < 0)) {
    throw('CDS end must be greater than or equal to 0');
  }

  if(defined($tl_start) && ($tl_start !~ /^\d+$/ || $tl_start < 1)) {
    throw('Translation start must be greater than or equal to 1');
  }

  if(defined($tl_end) && ($tl_end !~ /^\d+$/ || $tl_start < 0)) {
    throw('Translation end must be greater than or equal to 0');
  }

  return bless {'dbID'              => $dbID,
                'adaptor'           => $adaptor,
                'variation_feature' => $vf,
                'transcript'        => $tr,
                'pep_allele_string' => $pep_allele,
                'cdna_start'        => $cdna_start,
                'cdna_end'          => $cdna_end,
		'cds_start'         => $cds_start,
                'cds_end'           => $cds_end,
                'translation_start' => $tl_start,
                'translation_end'   => $tl_end,
                'consequence_type'  => $consequence_type,
				'codons'				=> $codons,}, $class;
}

sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
}


=head2 transcript_stable_id

  Example    : print $trvar->transcript_stable_id(), "\n";
  Description: Get the stable_id of the Transcript that is affected by this
               TranscriptVariation. This will NOT trigger lazy-loading of
	       the transcript.
  Returntype : string
  Exceptions : throw on bad argument
  Caller     : general
  Status     : At Risk

=cut

sub transcript_stable_id {
  my $self = shift;

  # If the transcript object has been loaded, get the stable_id from it. Otherwise, return the stable_id that is stored in this TranscriptVariation  
  return $self->{'transcript'}->stable_id() if (defined($self->{'transcript'}));
  return $self->{'_transcript_stable_id'};
}

=head2 transcript

  Arg [1]    : (optional) Bio::EnsEMBL::Transcript $transcript
  Example    : print $trvar->transcript()->stable_id(), "\n";
  Description: Getter/Setter for the Transcript that is affected by this
               TranscriptVariation.
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : throw on bad argument
  Caller     : general
  Status     : At Risk

=cut

sub transcript {
  my $self = shift;

  if(@_) {
    my $tr = shift;
    if(defined($tr) && (!ref($tr) || !$tr->isa('Bio::EnsEMBL::Transcript'))) {
      throw('Bio::EnsEMBL::Transcript argument expected');
    }
    $self->{'transcript'} = $tr;
  }
  else{
      #lazy-load the transcript object into the transcript_variation, and return it
      if (!defined $self->{'transcript'} && $self->{'adaptor'} && defined($self->{'_transcript_stable_id'})){
	  my $transcript_adaptor = $self->{'adaptor'}->db()->dnadb()->get_TranscriptAdaptor();
	  $self->transcript($transcript_adaptor->fetch_by_stable_id($self->{'_transcript_stable_id'}));
	  delete $self->{'_transcript_stable_id'};
      }
  }

  return $self->{'transcript'};
}


=head2 variation_feature

  Args       : none
  Example    : print $trvar->variation_feature()->variation_name(), "\n";
  Description: Getter for the VariationFeature associated with this
               transcript variation.
  Returntype : Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub variation_feature {
  my $self = shift;

  if(defined($self->{'_vf_id'}) && $self->{'adaptor'}){
      #lazy-load  from database on demand
      my $vf = $self->{'adaptor'}->db()->get_VariationFeatureAdaptor();
      return $vf->fetch_by_dbID($self->{'_vf_id'});
  }
  return $self->{'variation_feature'};
}




=head2 pep_allele_string

  Arg [1]    : string $allele (optional) 
               The new value to set the pep_allele_string attribute to
  Example    : $pep_allele_string = $obj->pep_allele_string()
  Description: Getter/Setter for the pep_allele_string attribute.  The
               pep allele string is a '/' delimited string of amino acid
               codes representing the change to the peptide made by this
               variation.  A '-' represents one half of a insertion/deletion.
               The reference allele (the one found in the ensembl peptide)
               should be first.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub pep_allele_string{
  my $self = shift;
  return $self->{'pep_allele_string'} = shift if(@_);
  return $self->{'pep_allele_string'};
}



=head2 cdna_start

  Arg [1]    : (optional) int $start
  Example    : $cdna_start = $trvar->cdna_start();
  Description: Getter/Setter for the start position of this variation on the
               transcript in CDNA coordinates.
  Returntype : int
  Exceptions : throw if $start is not an int
               throw if $start < 1
  Caller     : general
  Status     : Stable

=cut

sub cdna_start {
  my $self = shift;

  if(@_) {
    my $cdna_start = shift;
    if(defined($cdna_start) && ($cdna_start !~ /^\d+$/ || $cdna_start < 1)) {
      throw('cdna start must be an integer greater than 0');
    }
    $self->{'cdna_start'} = $cdna_start;
  }
  
  return $self->{'cdna_start'};
}



=head2 cdna_end

  Arg [1]    : (optional) int $end
  Example    : $cdna_end = $trvar->cdna_end();
  Description: Getter/Setter for the end position of this variation on the
               transcript in cdna coordinates.
  Returntype : int
  Exceptions : throw if $end is not an int
               throw if $end < 0
  Caller     : general
  Status     : Stable

=cut

sub cdna_end {
  my $self = shift;

  if(@_) {
    my $cdna_end = shift;
    if(defined($cdna_end) && ($cdna_end !~ /^\d+$/ || $cdna_end < 0)) {
      throw('cdna end must be an integer greater than or equal to 0');
    }
    $self->{'cdna_end'} = $cdna_end;
  }

  return $self->{'cdna_end'};
}


=head2 cds_start

  Arg [1]    : (optional) int $start
  Example    : $cds_start = $trvar->cds_start();
  Description: Getter/Setter for the start position of this variation on the
               transcript in CDS coordinates.
  Returntype : int
  Exceptions : throw if $start is not an int
               throw if $start < 1
  Caller     : general
  Status     : Stable

=cut

sub cds_start {
  my $self = shift;

  if(@_) {
    my $cds_start = shift;
    if(defined($cds_start) && ($cds_start !~ /^\d+$/ || $cds_start < 1)) {
      throw('cds start must be an integer greater than 0');
    }
    $self->{'cds_start'} = $cds_start;
  }
  
  return $self->{'cds_start'};
}



=head2 cds_end

  Arg [1]    : (optional) int $end
  Example    : $cds_end = $trvar->cds_end();
  Description: Getter/Setter for the end position of this variation on the
               transcript in cds coordinates.
  Returntype : int
  Exceptions : throw if $end is not an int
               throw if $end < 0
  Caller     : general
  Status     : Stable

=cut

sub cds_end {
  my $self = shift;

  if(@_) {
    my $cds_end = shift;
    if(defined($cds_end) && ($cds_end !~ /^\d+$/ || $cds_end < 0)) {
      throw('cds end must be an integer greater than or equal to 0');
    }
    $self->{'cds_end'} = $cds_end;
  }

  return $self->{'cds_end'};
}


=head2 translation_start

  Arg [1]    : (optional) int $tl_start
  Example    : $tl_start = $trvar->translation_start();
  Description: Getter/Setter for the start position of this variation on the
               translation of the associated transcript in peptide coordinates.
  Returntype : int
  Exceptions : throw if $start is not an int
               throw if $start < 0
  Caller     : general
  Status     : Stable

=cut

sub translation_start {
  my $self = shift;

  if(@_) {
    my $tl_start = shift;
    if(defined($tl_start) && ($tl_start !~ /^\d+$/ || $tl_start < 1)) {
      throw('translation start must be an integer greater than or equal to 1');
    }
    $self->{'translation_start'} = $tl_start;
  }
  
  return $self->{'translation_start'};
}


=head2 translation_end

  Arg [1]    : (optional) int $tl_end
  Example    : $tl_end = $trvar->translation_end();
  Description: Getter/Setter for the end position of this variation on the
               translation of the associated transcript in peptide coordinates.
  Returntype : int
  Exceptions : throw if $end is not an int
               throw if $end < 0
  Caller     : general
  Status     : Stable

=cut

sub translation_end {
  my $self = shift;

  if(@_) {
    my $tl_end = shift;
    if(defined($tl_end) && ($tl_end !~ /^\d+$/ || $tl_end < 0)) {
      throw('translation end must be an integer greater than or equal to 0');
    }
    $self->{'translation_end'} = $tl_end;
  }

  return $self->{'translation_end'};
}


=head2 consequence_type

  Arg [1]    : (optional) string $consequence_type
  Example    : if($tr_var->consequence_type()->[0] eq 'INTRONIC') { do_something(); }
  Description: Getter/Setter for the consequence type of this transcript variation.
               Allowed values are: 'ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST','FRAMESHIFT_CODING',
		  'NON_SYNONYMOUS_CODING','SPLICE_SITE','SYNONYMOUS_CODING','REGULATORY_REGION',
		  '5PRIME_UTR','3PRIME_UTR','INTRONIC','UPSTREAM','DOWNSTREAM'
  Returntype : ref to array of strings
  Exceptions : throw if provided argument is not one of the valid strings
  Caller     : general
  Status     : At Risk

=cut

sub consequence_type {
  my $self = shift;
 
  if(@_) {
      my $consequence_type = shift;
      if(defined($consequence_type)) {
	  $consequence_type = uc($consequence_type);  
	  if(!$CONSEQUENCE_TYPES{$consequence_type}) {
	      my $valid = join(',',map({"'$_'"} keys(%CONSEQUENCE_TYPES)));
	      throw("Type argument must be one of: $valid");
	  } 
	    
      }
      push @{$self->{'consequence_type'}}, $consequence_type;
  }
  
  return $self->{'consequence_type'};
}

=head2 display_consequence

  Args       : none
  Example    : $display_consequence = $tv->display_consequence();
  Description: Getter for the consequence type to display,
               when more than one
  Returntype : string
  Exceptions : throw on incorrect argument
  Caller     : webteam
  Status     : At Risk

=cut

sub display_consequence{
    my $self = shift;

    my $highest_priority;

    #get the value to display from the consequence_type attribute
    $highest_priority = 'INTERGENIC';
    foreach my $ct (@{$self->consequence_type}){
	if ($CONSEQUENCE_TYPES{$ct} < $CONSEQUENCE_TYPES{$highest_priority}){
	    $highest_priority = $ct;
	}
    }
    
    
    return $highest_priority;
}



=head2 codons

  Arg [1]    : (optional) string $codons
  Example    : $codons = $consequence_type->codons
  Description: Getter/Setter for the possible codons created in this transcript
  Returntype : "/"-separated string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub codons {
  my $self = shift;

  if(@_) {
    $self->{'codons'} = shift;
  }
  
  elsif(
	!(defined($self->{'codons'}))
	&& defined($self->transcript)
	&& defined($self->translation_start)
	&& defined($self->variation_feature)
	&& ((join ',', @{$self->consequence_type}) =~ /SYN|STOP|PAR/)
  ) {
	
	my @codons;
	
	# get the transcript and its translateable sequence
	my $tr = $self->transcript;
	my $cds = $tr->translateable_seq();
	
	$cds = "\L$cds";
	
	# get the VF
	my $vf = $self->variation_feature;
	
	# get codon start end and length in CDS terms
	my $codon_cds_start = $self->translation_start * 3 - 2;
	my $codon_cds_end   = $self->translation_end   * 3;
	my $codon_len = $codon_cds_end - $codon_cds_start + 1;
	
	# add reference codon to codon list
	#push @codons, substr($cds, $codon_cds_start-1, $codon_len);
	
	my $var_len = $vf->length;
	
	# fetch and expand the allele string
	my $allele_string = $vf->allele_string;
	expand(\$allele_string);
	
	my @alleles = split /\//, $allele_string;
	my $strand = $vf->strand;
	
	# if we need to flip strand to match the transcript
	if ($strand != $tr->strand()) {
		
		# flip feature onto same strand as transcript
		for (my $i = 0; $i < @alleles; $i++) {
		  reverse_comp(\$alleles[$i]);
		}
		
		$strand = $tr->strand();
	}
	
	# first allele is reference
	#shift @alleles;
	
	# iterate through remaining alleles
	foreach my $a(@alleles) {
	  
	  # make a copy of the cds
	  my $tmp_cds = $cds;
	  
	  # work out which coord to edit
	  my $exon_phase = $tr->start_Exon->phase;
	  
	  my $edit_coord = ($self->cdna_start - $tr->cdna_coding_start) + ($exon_phase > 0 ? $exon_phase : 0);
	  
	  # do the edit
	  substr($tmp_cds, $edit_coord, $var_len) = $a;
	  
	  # get the new codon
	  my $codon_str = substr($tmp_cds, $codon_cds_start-1, $codon_len + length($a)-$var_len);
	  
	  # add it to the list
	  push @codons, $codon_str;
	}
	
	# create slash-separated codon string
	$self->{'codons'} = join '/', @codons;
  }
  
  return $self->{'codons'};
}

=head2 affects_transcript

  Example    : if($tv->affects_transcript) { ... }
  Description: Indicates if this transcript variation affects the protein
               sequence of the affected transcript
  Returntype : boolean
  Exceptions : none
  Caller     : Internal
  Status     : At Risk

=cut

sub affects_transcript {
  my $self = shift;
  return defined($AFFECTS_PEPTIDE{$self->display_consequence});
}



=head2 _calc_transcript_coords

  Example    : $tv->_calc_transcript_coords
  Description: Internal method for calculating absent transcript-relative coords
  Returntype : none
  Exceptions : none
  Caller     : Internal
  Status     : At Risk

=cut

sub _calc_transcript_coords {
  my $self = shift;
  
  my $tr = $self->transcript;
  
  if(defined($tr)) {
	
	# get a transcript mapper
	my $tm = $self->{'_transcript_mapper'} || $tr->get_TranscriptMapper();
	$self->{'_transcript_mapper'} ||= $tm;
	
	# get the variation feature object
	my $vf = $self->variation_feature;
	
	# do pep coords
	unless(defined $self->{'translation_start'} && defined $self->{'translation_end'}) {
	  my @pep_coords = $tm->genomic2pep($vf->start,$vf->end,$vf->strand);
	  if(scalar @pep_coords == 1 && !$pep_coords[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
		$self->{'translation_start'} ||= $pep_coords[0]->start;
		$self->{'translation_end'} ||= $pep_coords[0]->end;
	  }
	}
	
	# do cdna coords
	unless(defined $self->{'cdna_start'} && defined $self->{'cdna_end'}) {
	  my @cdna_coords = $tm->genomic2cdna($vf->start,$vf->end,$vf->strand);
	  if(scalar @cdna_coords == 1 && !$cdna_coords[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
		$self->{'cdna_start'} ||= $cdna_coords[0]->start;
		$self->{'cdna_end'} ||= $cdna_coords[0]->end;
	  }
	}
	
	# do cds coords
	unless(defined $self->{'cds_start'} && defined $self->{'cds_end'}) {
	  my @cds_coords = $tm->genomic2cds($vf->start,$vf->end,$vf->strand);
	  if(scalar @cds_coords == 1 && !$cds_coords[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
		$self->{'cds_start'} ||= $cds_coords[0]->start;
		$self->{'cds_end'} ||= $cds_coords[0]->end;
	  }
	}
  }
}


1;
