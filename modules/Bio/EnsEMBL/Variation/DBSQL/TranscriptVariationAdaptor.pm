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
        $trv->type(), "\n";

  # retrieve all TranscriptVariations associated with a Transcript

  $tr = $tra->fetch_by_stable_id('ENST00000278995');
  foreach $trv (@{$trva->fetch_all_by_Transcript($tr)}) {
    print $trv->variation_feature->variation_name(), ' ', $trv->type(), "\n";
  }

  # retrieve all TranscriptVariations associated with a VariationFeature

  $vf = $vfa->fetch_by_dbID(99123);
  foreach $trv (@{$trva->fetch_all_by_VariationFeature($vf)}) {
    print $trva->transcript->stable_id(), ' ', $trva->type(), "\n";
  }

  # retrieve all TranscriptsVariations associated with a Variation

  $v = $v->fetch_by_name('rs1445');
  foreach $trv (@{$trva->fetch_all_by_Variation($v)}) {
    print $trva->transcript->stable_id(), ' ', $trva->type(), "\n";
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


=head2 fetch_by_dbID

  Arg [1]    : int $dbID
  Example    : $var = $var_adaptor->fetch_by_dbID(5526);
  Description: Retrieves a TranscriptVariation object via its internal
               identifier. If no such TranscriptVariation exists in the 
               database undef is returned.
  Returntype : Bio::EnsEMBL::Variation::TranscriptVariation
  Exceptions : throw if dbID arg is not defined
  Caller     : general

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  throw('dbID argument expected') if(!defined($dbID));

  my $sth = $self->prepare
    (q{SELECT tv.transcript_variation_id, tv.transcript_id,
              tv.variation_feature_id, tv.cdna_start, tv.cdna_end,
              tv.translation_start, tv.translation_end,
              tv.peptide_allele_string, tv.type
       FROM   transcript_variation tv
       WHERE  tv.transcript_variation_id = ?});
  $sth->execute($dbID);

  my $result = $self->_objs_from_sth($sth);
  $sth->finish();

  return undef if(!@$result);

  return $result->[0];
}



=head2 fetch_all_by_Transcript

  Arg [1]    : Bio::EnsEMBL::Transcript $tr
  Example    :
     $tr = $tr_adaptor->fetch_by_stable_id('ENST00000278995');
     @tr_vars = @{$tr_var_adaptor->fetch_all_by_Transcript($tr)});
  Description: Retrieves all TranscriptVariation objects associated with
               a provided Ensembl Transcript.
  Returntype : ref to list of Bio::EnsEMBL::Variation::TranscriptVariations
  Exceptions : throw on bad argument
  Caller     : general

=cut

sub fetch_all_by_Transcript {
  my $self = shift;
  my $tr = shift;

  if(!ref($tr) || !$tr->isa('Bio::EnsEMBL::Transcript')) {
    throw('Bio::EnsEMBL::Transcript argument expected');
  }

  if(!$tr->dbID()) {
    warning('Can not retrieve TranscriptVariations ' .
            'for transcript without dbID');
    return [];
  }

  my $sth = $self->prepare
    (q{SELECT tv.transcript_variation_id, tv.transcript_id,
              tv.variation_feature_id, tv.cdna_start, tv.cdna_end,
              tv.translation_start, tv.translation_end,
              tv.peptide_allele_string, tv.type
       FROM   transcript_variation tv
       WHERE  tv.transcript_id = ?});
  $sth->execute($tr->dbID());

  my $result = $self->_objs_from_sth($sth);

  $sth->finish();

  return $result;
}



=head2 fetch_all_by_VariationFeature

  Arg [1]    : Bio::EnsEMBL::Variation::VariationFeature $vf
  Example    :
     $vf = $vf_adaptor->fetch_by_dbID(1234);
     @tr_vars = @{$tr_var_adaptor->fetch_all_by_VariationFeature($vf)});
  Description: Retrieves all TranscriptVariation objects associated with
               a provided Ensembl variation feature.
  Returntype : ref to list of Bio::EnsEMBL::Variation::TranscriptVariations
  Exceptions : throw on bad argument
  Caller     : general

=cut

sub fetch_all_by_VariationFeature {
  my $self = shift;
  my $vf = shift;

  if(!ref($vf) || !$vf->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
    throw('Bio::EnsEMBL::Variation::VariationFeature argument expected');
  }

  if(!$vf->dbID()) {
    warning('Can not retrieve TranscriptVariations ' .
            'for variation feature without dbID');
    return [];
  }

  my $sth = $self->prepare
    (q{SELECT tv.transcript_variation_id, tv.transcript_id,
              tv.variation_feature_id, tv.cdna_start, tv.cdna_end,
              tv.translation_start, tv.translation_end,
              tv.peptide_allele_string, tv.type
       FROM   transcript_variation tv
       WHERE  tv.variation_feature_id = ?});
  $sth->execute($vf->dbID());

  my $result = $self->_objs_from_sth($sth);

  $sth->finish();

  return $result;
}



=head2 fetch_all_by_Variation

  Arg [1]    : Bio::EnsEMBL::Variation::Variation $v
  Example    :
     $v = $var_adaptor->fetch_by_name('rs234');
     @tr_vars = @{$tr_var_adaptor->fetch_all_by_Variation($v)});
  Description: Retrieves all TranscriptVariation objects associated with
               a provided Ensembl variation.
  Returntype : ref to list of Bio::EnsEMBL::Variation::TranscriptVariations
  Exceptions : throw on bad argument
  Caller     : general

=cut

sub fetch_all_by_Variation {
  my $self = shift;
  my $v = shift;

  if(!ref($v) || !$v->isa('Bio::EnsEMBL::Variation::Variation')) {
    throw('Bio::EnsEMBL::Variation::Variation argument expected');
  }

  if(!$v->dbID()) {
    warning('Can not retrieve TranscriptVariations ' .
            'for variation without dbID');
    return [];
  }

  my $sth = $self->prepare
    (q{SELECT tv.transcript_variation_id, tv.transcript_id,
              tv.variation_feature_id, tv.cdna_start, tv.cdna_end,
              tv.translation_start, tv.translation_end,
              tv.peptide_allele_string, tv.type
       FROM   transcript_variation tv, variation_feature vf
       WHERE  tv.variation_feature_id = vf.variation_feature_id
       AND    vf.variation_id = ?});
  $sth->execute($v->dbID());

  my $result = $self->_objs_from_sth($sth);

  $sth->finish();

  return $result;
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
      $pep_allele, $type);

  $sth->bind_columns(\$trv_id, \$tr_id, \$vf_id, \$cdna_start, \$cdna_end,
                     \$tl_start, \$tl_end, \$pep_allele, \$type);


  my %tr_hash;
  my %vf_hash;

  my @results;

  # construct all of the TranscriptVariation objects

  while($sth->fetch()) {
    my $trv = Bio::EnsEMBL::Variation::TranscriptVariation->new
      (-dbID => $trv_id,
       -adaptor => $self,
       -cdna_start => $cdna_start,
       -cdna_end   => $cdna_end,
       -translation_start => $tl_start,
       -translation_end => $tl_end,
       -pep_allele_string => $pep_allele,
       -type => $type);

    $tr_hash{$tr_id} ||= [];
    $vf_hash{$vf_id} ||= [];
    push @{$tr_hash{$tr_id}}, $trv;
    push @{$vf_hash{$vf_id}}, $trv;

    push @results, $trv;
  }

  # load all transcripts and variation features with one query -
  # much faster than individual queries

  my $tra = $self->db()->get_TranscriptAdaptor();
  my $vfa = $self->db()->get_VariationFeatureAdaptor();

  my @tr_ids = keys %tr_hash;
  my @vf_ids = keys %vf_hash;

  my @trs = @{$tra->fetch_all_by_dbID_list(\@tr_ids)};
  my @vfs = @{$vfa->fetch_all_by_dbID_list(\@vf_ids)};


  # add the transcripts and variation features to the
  # already constructed transcript variation objects

  foreach my $tr (@trs) {
    foreach my $trv (@{$tr_hash{$tr->dbID()}}) {
      $trv->transcript($tr);
    }
  }
  foreach my $vf (@vfs) {
    foreach my $trv (@{$vf_hash{$vf->dbID()}}) {
      $trv->variation_feature($vf);
    }
  }

  return \@results;
}




1;
