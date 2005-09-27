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

  # retrieve all TranscriptsVariations associated with a Variation

  $v = $v->fetch_by_name('rs1445');
  foreach $trv (@{$trva->fetch_all_by_Variation($v)}) {
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

=cut

sub fetch_all_by_Transcripts{
    my $self = shift;
    my $transcript_ref = shift;
    
    if (ref($transcript_ref) ne 'ARRAY'){
	throw('Array Bio::EnsEMBL::Transcript expected');
    }
    
    my %tr_by_id;
    %tr_by_id = map {$_->dbID(), $_} @{$transcript_ref};
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

=cut

sub fetch_all_by_VariationFeatures {
  my $self = shift;
  my $vf_ref = shift;

  if(ref($vf_ref) ne 'ARRAY') {
    throw('Array Bio::EnsEMBL::Variation::VariationFeature expected');
  }

  my %vf_by_id;

  %vf_by_id = map {$_->dbID(), $_ } @$vf_ref;
  my $instr = join (",",keys( %vf_by_id));
  my $tvs = $self->generic_fetch( "tv.variation_feature_id in ( $instr )" );
  for my $tv ( @$tvs ) {
      #add to the variation feature object all the transcript variations
      $vf_by_id{ $tv->{'_vf_id'} }->add_TranscriptVariation( $tv );
      delete $tv->{'_vf_id'}; #remove the variation_feature_id from the transcript_variation object    
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

    #added new attribute, splice_site
    my %splice_sites = %Bio::EnsEMBL::Variation::ConsequenceType::SPLICE_SITES; #get the hash with the valid SPLICE_SITES
    my $splice_site;
    my @types = split(',',$consequence_type); #get the different consequence types and separate in 2 different attributes: slice_site and
    #there is a splice_site type
    if (@types > 1){
	$splice_site = shift @types; #relying in the order of the SET column in the database
    }
    $consequence_type = shift @types if (@types > 0);

    my $trv = Bio::EnsEMBL::Variation::TranscriptVariation->new_fast
      ( { 'dbID' => $trv_id,
	  'adaptor' => $self,
	  'cdna_start' => $cdna_start,
	  'cdna_end'   => $cdna_end,
	  'translation_start' => $tl_start,
	  'translation_end' => $tl_end,
	  'pep_allele_string' => $pep_allele,
	  'consequence_type' => $consequence_type,
          'splice_site'      => $splice_site || ''} );

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
