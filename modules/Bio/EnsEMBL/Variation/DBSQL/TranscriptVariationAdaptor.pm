#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor

=head1 SYNOPSIS

  # connect to variation database
  $vdb = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(...);

  # connect to core database
  $db  = Bio::EnsEMBL::DBSQL::DBAdaptor->new(...);


  # tell the variation database where to obtain core database data
  $vdb->dnadb($db);

  $trva = $vdb->get_TranscriptVariationAdaptor();
  $vfa  = $vdb->get_VariationFeatureAdaptor();
  $tra  = $db->get_TranscriptAdaptor();


  # retrieve a TranscriptVariation by its internal identifier

  $trv = $tra->fetch_by_dbID(552);
  print $trv->transcript()->stable_id(), ' ',
        $trv->variation_feature()->variation_name(), ' ',
        $trv->consequence_type(), "\n";

  # retrieve all TranscriptVariations associated with a Transcript

  $tr = $tra->fetch_by_stable_id('ENST00000278995');
  foreach $trv (@{$trva->fetch_all_by_Transcripts([$tr])}) {
    print $trv->variation_feature->variation_name(), ' ', $trv->consequence_type(), "\n";
  }

  # retrieve all TranscriptVariations associated with a VariationFeature

  $vf = $vfa->fetch_by_dbID(99123);
  foreach $trv (@{$trva->fetch_all_by_VariationFeatures($vf)}) {
    print $trva->transcript->stable_id(), ' ', $trva->consequence_type(), "\n";
  }


=head1 DESCRIPTION

This adaptor provides database connectivity for TranscriptVariation objects.
TranscriptVariations which represent an association between a variation and
a Transcript may be retrieved from the Ensembl variation database via several
means using this module.

=head1 AUTHOR - Graham McVicker

=head1 CONTACT

Post questions to the Ensembl development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp expand);
use Bio::EnsEMBL::Utils::TranscriptAlleles qw(type_variation);

use Bio::EnsEMBL::Variation::TranscriptVariation;

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');


=head2 fetch_all_by_Transcripts

    Arg[1]      : listref of Bio::EnsEMBL::Transcript
    Example     : $tr = $ta->fetch_by_stable_id('ENST00000278995');
                  @tr_vars = @{$tr_var_adaptor->fetch_all_by_Transcripts([$tr])};
    Description : Retrieves all TranscriptVariation objects associated with
                  provided Ensembl Transcript. Attaches them to the TranscriptVariation
    ReturnType  : ref to list of Bio::EnsEMBL::Variation::TranscriptVariations
    Exceptions  : throw on bad argument
    Caller      : general
    Status      : Stable

=cut

sub fetch_all_by_Transcripts{
    my $self = shift;
    my $transcript_ref = shift;

    if (ref($transcript_ref) ne 'ARRAY' or ! $transcript_ref->[0]->isa('Bio::EnsEMBL::Transcript')){
      throw('Array Bio::EnsEMBL::Transcript expected');
    }
    
    my %tr_by_id;

    foreach my $tr (@{$transcript_ref}) {
      if (!$tr->isa('Bio::EnsEMBL::Transcript')){
	throw('Bio::EnsEMBL::Transcript is expected');
      }
      $tr_by_id{$tr->dbID()} = $tr;
    }
    my $instr = join (",", keys( %tr_by_id));
    my $transcript_variations = $self->generic_fetch( "tv.transcript_id in ( $instr )");
    for my $tv (@{$transcript_variations}){
	#add to the TranscriptVariation object all the Transcripts
	$tv->{'transcript'} = $tr_by_id{$tv->{'_transcript_id'}};
	delete $tv->{'_transcript_id'}; #remove the transcript_id from the transcript_variation object
    }
    return $transcript_variations;
}

=head2 fetch_all_by_VariationFeatures

  Arg [1]    : Bio::EnsEMBL::Variation::VariationFeature $vf
  Example    :
     $vf = $vf_adaptor->fetch_by_dbID(1234);
     @tr_vars = @{$tr_var_adaptor->fetch_all_by_VariationFeatures([$vf])});
  Description: Retrieves all TranscriptVariation objects associated with
               provided Ensembl variation features. Attaches them to the given variation
               features.
  Returntype : ref to list of Bio::EnsEMBL::Variation::TranscriptVariations
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_VariationFeatures {
  my $self = shift;
  my $vf_ref = shift;

  if(ref($vf_ref) ne 'ARRAY') {
    throw('Array Bio::EnsEMBL::Variation::VariationFeature expected');
  }

  my %vf_by_id;
  my $tvs;
  
  my @have_dbID;
  my @no_dbID;
  
  foreach my $vf(@$vf_ref) {
	if(defined $vf->dbID) {
		push @have_dbID, $vf;
	}
	
	else {
		push @no_dbID, $vf;
	}
  }
  
  if(scalar @have_dbID) {

	%vf_by_id = map {$_->dbID(), $_ } @have_dbID;
	my $instr = join (",",keys( %vf_by_id));
	$tvs = $self->generic_fetch( "tv.variation_feature_id in ( $instr )" );
	for my $tv ( @$tvs ) {
	  #add to the variation feature object all the transcript variations
	  $vf_by_id{ $tv->{'_vf_id'} }->add_TranscriptVariation( $tv );
	  delete $tv->{'_vf_id'}; #remove the variation_feature_id from the transcript_variation object    
	}
  }
  
  if(scalar @no_dbID) {
	
	foreach my $vf(@no_dbID) {
		
		my $slice = $vf->feature_Slice;
		
		my @transcripts = @{$slice->get_all_Transcripts()};
		
		while(my $transcript = shift @transcripts) {
			
			# expand the allele string
			my $allele_string = $vf->allele_string;
			expand(\$allele_string);
			
			my @alleles = split /\//, $allele_string;
			my $strand = $vf->strand;
			
			# if we need to flip strand to match the transcript
			if ($strand != $transcript->strand()) {
				
				# flip feature onto same strand as transcript
				for (my $i = 0; $i < @alleles; $i++) {
				  reverse_comp(\$alleles[$i]);
				}
				
				$strand = $transcript->strand();
			}
			
			# shift off the reference allele
			shift @alleles;
			
			# convert the start and end to slice coords
			my $start = $vf->start - $slice->start + 1;
			my $end = $vf->end - $slice->start + 1;
			
			# create a consequence type object
			my $consequence_type = Bio::EnsEMBL::Variation::ConsequenceType->new(
				$transcript->dbID,
				undef,				# this is the VF dbID, we don't have one
				$start,
				$end,
				$strand,
				\@alleles
			);
			
			my $consequences;
			
			if($allele_string !~ /\+/) {
				$consequences = type_variation($transcript, "", $consequence_type);
			}
			
			foreach my $con(@$consequences) {
				my $trv = Bio::EnsEMBL::Variation::TranscriptVariation->new_fast( {
					'dbID' 				=> undef,
					'adaptor' 			=> $self,
					'cdna_start'		=> $con->cdna_start,
					'cdna_end'			=> $con->cdna_end,
					'translation_start' => $con->aa_start,
					'translation_end'	=> $con->aa_end,
					'pep_allele_string' => join("/", @{$con->aa_alleles || []}),
					'consequence_type'	=> $con->type
				} );
		  
				$trv->{'_vf_id'} = undef; #add the variation feature
				$trv->{'_transcript_id'} = $transcript->dbID; #add the transcript id
				
				push @{$tvs}, $trv;
				$vf->add_TranscriptVariation($trv);
			}
		}
	}
  }
  
  
  return $tvs;
}


#
# internal method responsible for constructing transcript variation objects
# from an executed statement handle.  Ordering of columns in executed statement
# must be consistant with function implementation.
#
sub _objs_from_sth {
  my $self = shift;
  my $sth = shift;

  my ($trv_id, $tr_id, $vf_id, $cdna_start, $cdna_end, $tl_start, $tl_end,
      $pep_allele, $consequence_type);

  $sth->bind_columns(\$trv_id, \$tr_id, \$vf_id, \$cdna_start, \$cdna_end,
                     \$tl_start, \$tl_end, \$pep_allele, \$consequence_type);


  my %tr_hash;
  my %vf_hash;

  my @results;

  # construct all of the TranscriptVariation objects

  while($sth->fetch()) {

      my @consequences = split /,/,$consequence_type;
    
      my $trv = Bio::EnsEMBL::Variation::TranscriptVariation->new_fast
	  ( { 'dbID' => $trv_id,
	      'adaptor' => $self,
	      'cdna_start' => $cdna_start,
	      'cdna_end'   => $cdna_end,
	      'translation_start' => $tl_start,
	      'translation_end' => $tl_end,
	      'pep_allele_string' => $pep_allele,
	      'consequence_type' => \@consequences} );

    $trv->{'_vf_id'} = $vf_id; #add the variation feature
    $trv->{'_transcript_id'} = $tr_id; #add the transcript id
    push @results, $trv;
  }


  return \@results;
}

sub _tables {return ['transcript_variation','tv'];}

sub _columns {
    return qw (tv.transcript_variation_id tv.transcript_id
	       tv.variation_feature_id tv.cdna_start tv.cdna_end
	       tv.translation_start tv.translation_end
	       tv.peptide_allele_string tv.consequence_type
	       );
}
1;
