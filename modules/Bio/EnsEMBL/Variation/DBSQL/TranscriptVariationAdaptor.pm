=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

=cut

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor

=head1 SYNOPSIS
  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_registry_from_db(-host => 'ensembldb.ensembl.org', -user => 'anonymous');

  my $ta  = $registry->get_adaptor('human', 'core', 'Transcript');
  my $tva = $registry->get_adaptor('human', 'variation', 'TranscriptVariation');
  my $va  = $registry->get_adaptor('human', 'variation', 'Variation');
  my $vfa = $registry->get_adaptor('human', 'variation', 'VariationFeature');

  # fetch all TranscriptVariations related to a Transcript
  my $transcript = $ta->fetch_by_stable_id('ENST00000380152');
  for my $tv (@{ $tva->fetch_all_by_Transcripts([$transcript]) }) {
  print $tv->display_consequence, "\n";
  }

  # fetch all TranscriptVariations related to a VariationFeature
  my $vf = $vfa->fetch_all_by_Variation($va->fetch_by_name('rs669'))->[0];
  for my $tv (@{ $tva->fetch_all_by_VariationFeatures([$vf]) }) {
  print $tv->display_consequence, "\n";
  }

  # fetch all TranscriptVariations related to a Translation
  for my $tv (@{ $tva->fetch_all_by_translation_id('ENSP00000447797') }) {
  foreach my $allele (keys %{$tv->hgvs_protein}) {
    my $hgvs_notation = $tv->hgvs_protein->{$allele} || 'hgvs notation is NA';
    print "$allele $hgvs_notation\n";
  }
  }

  # fetch all TranscriptVariations related to a Translation with given SO terms
  for my $tv (@{ $tva->fetch_all_by_translation_id_SO_terms('ENSP00000447797', ['missense_variant']) }) {
  foreach my $allele (keys %{$tv->hgvs_protein}) {
    my $hgvs_notation = $tv->hgvs_protein->{$allele} || 'hgvs notation is NA';
    print "$allele $hgvs_notation\n";
  }
  }

=head1 DESCRIPTION

This adaptor allows you to fetch TranscriptVariation objects either by the Transcripts
the associated VariationFeature falls in, or by VariationFeature directly. Storing
TranscriptVariation objects in a variation schema database is also supported. In the
database there will a separate row for each alternative allele of a TranscriptVariation, 
but the methods here will fetch all alleles associated with the TranscriptVariation
at once.

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor;

use Bio::EnsEMBL::Variation::TranscriptVariation;
use Bio::EnsEMBL::Variation::TranscriptVariationAllele;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(MAX_DISTANCE_FROM_TRANSCRIPT);
use Bio::EnsEMBL::Variation::Utils::Constants qw(%OVERLAP_CONSEQUENCES);
use Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix;

use base qw(Bio::EnsEMBL::Variation::DBSQL::VariationFeatureOverlapAdaptor);
use Scalar::Util qw(weaken);

our $DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME  = 1;

=head2 store

  Arg [1]  : Bio::EnsEMBL::Variation::TranscriptVariation $tv
  Description: Store the TranscriptVariation in the database
  Status   : At risk

=cut

sub store {
  my ($self, $tv, $mtmp) = @_;

  my $write_data = $self->_get_write_data($tv);
  $self->_store_write_data($write_data);

  $self->_store_mtmp_write_data($self->_get_mtmp_write_data_from_tv_write_data($write_data)) if $mtmp;
}

=head2 fetch_all_by_Transcripts_SO_terms

  Arg [1]  : listref of Bio::EnsEMBL::Transcripts
  Arg [2]  : listref of SO terms
  Description: Fetch all germline TranscriptVariations associated with the
         given list of Transcripts with consequences with given SO terms
  Returntype : listref of Bio::EnsEMBL::Variation::TranscriptVariations
  Status   : At risk

=cut

sub fetch_all_by_Transcripts_SO_terms {
  my ($self, $transcripts, $terms) = @_;
  my $constraint = $self->_get_consequence_constraint($terms);
  return $self->fetch_all_by_Transcripts_with_constraint($transcripts, $constraint.' AND somatic = 0');
}

=head2 fetch_all_somatic_by_Transcripts_SO_terms

  Arg [1]  : listref of Bio::EnsEMBL::Transcripts
  Arg [2]  : listref of SO terms
  Description: Fetch all somatic TranscriptVariations associated with the
         given list of Transcripts with consequences with given SO terms
  Returntype : listref of Bio::EnsEMBL::Variation::TranscriptVariations
  Status   : At risk

=cut

sub fetch_all_somatic_by_Transcripts_SO_terms {
  my ($self, $transcripts, $terms) = @_;
  my $constraint = $self->_get_consequence_constraint($terms);
  return $self->fetch_all_by_Transcripts_with_constraint($transcripts, $constraint.' AND somatic = 1');
}

=head2 fetch_all_by_VariationFeatures_SO_terms

  Arg [1]    : listref of Bio::EnsEMBL::Variation::VariationFeatures
  Arg [2]    : listref of SO terms
  Description: Fetch all germline TranscriptVariations associated with the
               given list of VariationFeatures with consequences with given
               SO terms
  Returntype : listref of Bio::EnsEMBL::Variation::TranscriptVariations
  Status   : At risk

=cut

sub fetch_all_by_VariationFeatures_SO_terms {
  my ($self, $vfs, $transcripts, $terms, $without_children, $included_so) = @_;
  my $constraint = $self->_get_consequence_constraint($terms, $without_children, $included_so);
  if (!$constraint) {
    return [];
  }
  return $self->SUPER::fetch_all_by_VariationFeatures_with_constraint($vfs, $transcripts, $constraint);
}

=head2 count_all_by_VariationFeatures_SO_terms

  Arg [1]    : listref of Bio::EnsEMBL::Variation::VariationFeatures
  Arg [2]    : listref of SO terms
  Description: Count TranscriptVariations associated with given
               VariationFeatures with consequences with given SO terms
  Returntype : int
  Status     : At risk

=cut

sub count_all_by_VariationFeatures_SO_terms {
  my ($self, $vfs, $transcripts, $terms, $included_so) = @_;
  my $constraint = $self->_get_consequence_constraint($terms, 1, $included_so);
  if (!$constraint) {
    return 0;
  }
  return $self->SUPER::count_all_by_VariationFeatures_with_constraint($vfs, $transcripts, $constraint);
}

=head2 fetch_all_by_Transcripts

  Arg [1]  : listref of Bio::EnsEMBL::Transcripts
  Description: Fetch all germline TranscriptVariations associated with the
         given list of Transcripts
  Returntype : listref of Bio::EnsEMBL::Variation::TranscriptVariations
  Status   : Stable

=cut

sub fetch_all_by_Transcripts {
  my ($self, $transcripts) = @_;
  return $self->fetch_all_by_Transcripts_with_constraint($transcripts, 'somatic = 0');
}

=head2 fetch_all_somatic_by_Transcripts

  Arg [1]  : listref of Bio::EnsEMBL::Transcripts
  Description: Fetch all somatic TranscriptVariations associated with the
         given list of Transcripts
  Returntype : listref of Bio::EnsEMBL::Variation::TranscriptVariations
  Status   : Stable

=cut

sub fetch_all_somatic_by_Transcripts {
  my ($self, $transcripts) = @_;
  return $self->fetch_all_by_Transcripts_with_constraint($transcripts, 'somatic = 1');
}

=head2 fetch_all_by_translation_id

  Arg[1]   : String $translation_id
         The stable identifier of the translation
  Description: Fetch all germline TranscriptVariations associated with the given Translation
  Returntype : listref of Bio::EnsEMBL::Variation::TranscriptVariation
  Status   : At Risk
=cut

sub fetch_all_by_translation_id {
  my ($self, $translation_id) = @_;
  my $transcript = $self->_transcript($translation_id);
  my $all_tvs = $self->fetch_all_by_Transcripts([$transcript]);
  return $self->_transcript_variations_on_protein($all_tvs);
}

=head2 fetch_all_somatic_by_translation_id

  Arg[1]   : String $translation_id
         The stable identifier of the translation.
  Description: Fetch all somatic TranscriptVariations associated with the given Translation
  Returntype : listref of Bio::EnsEMBL::Variation::TranscriptVariation
  Status   : At Risk
=cut

sub fetch_all_somatic_by_translation_id {
  my ($self, $translation_id) = @_;
  my $transcript = $self->_transcript($translation_id);
  my $all_tvs = $self->fetch_all_somatic_by_Transcripts([$transcript]);
  return $self->_transcript_variations_on_protein($all_tvs);
}

=head2 fetch_all_by_translation_id_SO_terms

  Arg[1]    : String $translation_id
          The stable identifier of the translation
  Arg[2]    : listref of SO terms
  Description : Fetch all germline TranscriptVariations associated with the given Translation
          and having consequence types as given in the input list of SO terms
  Returntype  : listref of Bio::EnsEMBL::Variation::TranscriptVariation
  Status    : At Risk
=cut

sub fetch_all_by_translation_id_SO_terms {
  my ($self, $translation_id, $terms) = @_;
  my $transcript = $self->_transcript($translation_id);
  my $all_tvs = $self->fetch_all_by_Transcripts_SO_terms([$transcript], $terms);
  return $self->_transcript_variations_on_protein($all_tvs);
}

=head2 fetch_all_somatic_by_translation_id_SO_terms

  Arg[1]   : String $translation_id
         The stable identifier of the translation
  Arg[2]   : listref of SO terms
  Description: Fetch all somatic TranscriptVariations associated with the given Translation
         and having consequence types as given in the input list of SO terms
  Returntype : listref of Bio::EnsEMBL::Variation::TranscriptVariation
  Status   : At Risk 
=cut

sub fetch_all_somatic_by_translation_id_SO_terms {
  my ($self, $translation_id, $terms) = @_;
  my $transcript = $self->_transcript($translation_id);
  my $all_tvs = $self->fetch_all_somatic_by_Transcripts_SO_terms([$transcript], $terms);
  return $self->_transcript_variations_on_protein($all_tvs);
}

# Returns the associated Transcript for a given translation id 
sub _transcript {
  my ($self, $translation_id) = @_;
  my $transcript_adaptor = $self->db()->dnadb()->get_TranscriptAdaptor(); 
  my $transcript = $transcript_adaptor->fetch_by_translation_stable_id($translation_id);
  return $transcript;
}

# Returns listref of TranscriptVariations whose coordinates can be mapped to the protein sequence 
sub _transcript_variations_on_protein {
  my ($self, $all_tvs) = @_;  
  my @tvs;
  foreach my $tv (@$all_tvs) {
    if ($tv->translation_start && $tv->translation_end) {
      push(@tvs, $tv);
    }
  }
  return \@tvs;
}

=head2 fetch_all_by_Transcripts_with_constraint

  Arg [1]  : listref of Bio::EnsEMBL::Transcripts
  Arg [2]  : extra SQL constraint for the query
  Description: Fetch all TranscriptVariations associated with the
         given list of Transcripts
  Returntype : listref of Bio::EnsEMBL::Variation::TranscriptVariations
  Status   : At risk

=cut

sub fetch_all_by_Transcripts_with_constraint {
  my ($self, $transcripts, $constraint) = @_;
  
  return $self->SUPER::fetch_all_by_Features_with_constraint($transcripts, $constraint);
}

sub _fetch_all_by_VariationFeatures_no_DB {
  my ($self, $vfs, $features, $constraint, $dont_add_to_vf) = @_;

  # get features?
  if(!$features || !@$features) {
    my $slices = $self->_get_ranged_slices_from_VariationFeatures($vfs);
    @$features = map {@{$_->get_all_Transcripts(1)}} @$slices;
  }
  
  my @return;
  
  foreach my $f(@$features) {

    my $f_slice = $f->slice;
    
    foreach my $vf(@$vfs) {
      my $vfo = Bio::EnsEMBL::Variation::TranscriptVariation->new(
        -variation_feature  => $vf,
        -transcript         => $f,
        -adaptor            => $self,
        -no_ref_check      => 1,
        -no_transfer       => ($vf->slice + 0) == ($f_slice + 0)
      );
      
      $vf->add_TranscriptVariation($vfo) unless $dont_add_to_vf;

      push @return, $vfo;
    }
  }
  
  return \@return;
}

sub _objs_from_sth {
  my ($self, $sth) = @_;

  #warn $sth->sql;

  my (
    $transcript_variation_id,
    $variation_feature_id, 
    $feature_stable_id, 
    $allele_string,
    $consequence_types,
    $cds_start,
    $cds_end,
    $cdna_start,
    $cdna_end,
    $translation_start,
    $translation_end,
    $distance_to_transcript,
    $codon_allele_string,
    $pep_allele_string,
    $hgvs_genomic,
    $hgvs_transcript,
    $hgvs_protein,
    $polyphen_prediction,
    $polyphen_score,
    $sift_prediction,
    $sift_score,
    $display
  );
  
  $sth->bind_columns(
    \$transcript_variation_id,
    \$variation_feature_id, 
    \$feature_stable_id, 
    \$allele_string,
    \$consequence_types,
    \$cds_start,
    \$cds_end,
    \$cdna_start,
    \$cdna_end,
    \$translation_start,
    \$translation_end,
    \$distance_to_transcript,
    \$codon_allele_string,
    \$pep_allele_string,
    \$hgvs_genomic,
    \$hgvs_transcript,
    \$hgvs_protein,
    \$polyphen_prediction,
    \$polyphen_score,
    \$sift_prediction,
    \$sift_score,
    \$display
  );
  
  my %tvs;
  
  while ($sth->fetch) {
    
    my ($ref_allele, $alt_allele)   = split /\//, $allele_string;
    my ($ref_codon, $alt_codon)   = split /\//, $codon_allele_string || '';
    my ($ref_pep, $alt_pep)     = split /\//, $pep_allele_string || '';
    
    # for HGMD mutations etc. just set the alt allele to the ref allele
    $alt_allele ||= $ref_allele;
    
    # for synonymous mutations the peptides are the same and 
    # there is no / in the string
    $alt_pep ||= $ref_pep;
    
    # for TranscriptVariations with multiple alternative alleles
    # there will be multiple rows in the database, so we construct
    # the TV object and the reference allele object when we see 
    # the first row, but then only add extra allele objects when 
    # we see further rows, we track existing TVs in the %tvs hash, 
    # keyed by variation_feature_id and feature_stable_id

    my $key = $variation_feature_id.'_'.$feature_stable_id;
    
    my $tv = $tvs{$key};
    
    unless ($tv) {
      $tv = Bio::EnsEMBL::Variation::TranscriptVariation->new_fast({
        _variation_feature_id   => $variation_feature_id,
        _feature_stable_id    => $feature_stable_id,
        cds_start         => $cds_start,
        cds_end         => $cds_end,
        cdna_start        => $cdna_start,
        cdna_end        => $cdna_end,
        translation_start     => $translation_start,
        translation_end     => $translation_end,
        distance_to_transcript  => $distance_to_transcript,
        display         => $display,
        adaptor         => $self,
      });
      
      $tvs{$key} = $tv;
      
      my $ref_allele = Bio::EnsEMBL::Variation::TranscriptVariationAllele->new_fast({
        is_reference        => 1,
        variation_feature_seq     => $ref_allele,
        transcript_variation    => $tv, 
        codon             => $ref_codon,
        peptide           => $ref_pep, 
        dbID            => $transcript_variation_id,
      });

      $tv->add_TranscriptVariationAllele($ref_allele);
    }
     
    #my $overlap_consequences = $self->_transcript_variation_consequences_for_set_number($consequence_types);

    my $overlap_consequences = [ map { $OVERLAP_CONSEQUENCES{$_} } split /,/, ($consequence_types || 'sequence_variant') ];
    
    my $allele = Bio::EnsEMBL::Variation::TranscriptVariationAllele->new_fast({
      is_reference        => 0,
      variation_feature_seq     => $alt_allele,
      transcript_variation    => $tv, 
      codon             => $alt_codon,
      peptide           => $alt_pep,
      hgvs_genomic        => $hgvs_genomic,
      hgvs_transcript       => $hgvs_transcript,
      hgvs_protein        => $hgvs_protein,
      overlap_consequences    => $overlap_consequences, 
      polyphen_prediction     => $polyphen_prediction,
      polyphen_score        => $polyphen_score,
      sift_prediction       => $sift_prediction, 
      sift_score          => $sift_score, 
      dbID            => $transcript_variation_id,
    });
    
    $tv->add_TranscriptVariationAllele($allele);
  }
  
  return [values %tvs];
}

sub _tables {
  return (
    ['transcript_variation', 'tv']
  );
}

sub _columns {
  return qw(
    transcript_variation_id 
    variation_feature_id 
    feature_stable_id 
    allele_string 
    consequence_types 
    cds_start 
    cds_end 
    cdna_start 
    cdna_end 
    translation_start 
    translation_end 
    distance_to_transcript 
    codon_allele_string 
    pep_allele_string 
    hgvs_genomic 
    hgvs_transcript 
    hgvs_protein 
    polyphen_prediction 
    polyphen_score 
    sift_prediction
    sift_score
    display
  );
}

sub _write_columns {
  return qw(
    variation_feature_id 
    feature_stable_id 
    allele_string
    somatic 
    consequence_types 
    cds_start 
    cds_end 
    cdna_start 
    cdna_end 
    translation_start 
    translation_end 
    distance_to_transcript 
    codon_allele_string 
    pep_allele_string 
    hgvs_genomic 
    hgvs_transcript 
    hgvs_protein 
    polyphen_prediction 
    polyphen_score 
    sift_prediction
    sift_score
    display
  );
}

sub _mtmp_write_columns {
  my $self = shift;

  if(!exists($self->{_mtmp_write_columns})) {
    my %mtmp_exclude = map {$_ => 1} $self->_mtmp_remove_columns();

    $self->{_mtmp_write_columns} = [grep {!$mtmp_exclude{$_}} $self->_write_columns()];
  }

  return @{$self->{_mtmp_write_columns}};
}

sub _mtmp_remove_columns {
  return qw(transcript_variation_id hgvs_genomic hgvs_protein hgvs_transcript somatic codon_allele_string);
}

sub _default_where_clause {

  my $self = shift;
  
  my $clause = ' tv.display = 1 ' unless $self->db->include_failed_variations();
  return $clause;

}

## helper routines used by store
################################

# returns an arrayref of arrayrefs
# each one corresponding to a row of data in the table
sub _get_write_data {
  my ($self, $tv) = @_;

  my @return;

  my $vf = $tv->base_variation_feature;

  # cache up front, quicker when multiple alleles
  my (
    $vfID,
    $trID,
    $somatic,
    $cds_start,
    $cds_end,
    $cdna_start,
    $cdna_end,
    $translation_start,
    $translation_end,
    $distance_to_transcript,
    $display
  ) = (
    $vf->dbID,
    $tv->feature->stable_id,
    $vf->is_somatic,
    undef,
    undef,
    undef,
    undef,
    undef,
    undef,
    $tv->distance_to_transcript,
    $vf->display
  );

  foreach my $allele(@{$tv->get_all_alternate_TranscriptVariationAlleles}) {

    # use pre-predicate data to avoid running costly subs
    # we also need the HGVS tva in case shifting has changed things (this might be the original tva anyway, no extra cost)
    my $pre = $allele->_pre_consequence_predicates;
    my $hgvs_pre = $allele->_hgvs_tva ? $allele->_hgvs_tva->_pre_consequence_predicates : {};

    my (
      $codon_allele_string, $pep_allele_string,
      $hgvs_transcript, $hgvs_protein,
      $polyphen_prediction, $polyphen_score,
      $sift_prediction, $sift_score
    );

    # exon
    if($pre->{exon}) {
      $cdna_start ||= $tv->cdna_start;
      $cdna_end   ||= $tv->cdna_end;
    }

    # coding
    if($pre->{coding}) {
      $cds_start         ||= $tv->cds_start;
      $cds_end           ||= $tv->cds_end;
      $translation_start ||= $tv->translation_start;
      $translation_end   ||= $tv->translation_end;

      $codon_allele_string = $allele->codon_allele_string;
      $pep_allele_string   = $allele->pep_allele_string;
      $polyphen_prediction = $allele->polyphen_prediction;
      $polyphen_score      = $allele->polyphen_score;
      $sift_prediction     = $allele->sift_prediction;
      $sift_score          = $allele->sift_score;
    }

    # HGVS-specific
    if($hgvs_pre->{within_feature}) {
      $hgvs_transcript = $allele->hgvs_transcript;
    }

    if($hgvs_pre->{coding}) {
      $hgvs_protein = $allele->hgvs_protein;
    }

    push @return, [
      $vfID,
      $trID,
      $allele->allele_string,
      $somatic,
      (join ',', map { $_->SO_term } @{ $allele->get_all_OverlapConsequences }),
      $cds_start, 
      $cds_end,
      $cdna_start,
      $cdna_end,
      $translation_start,
      $translation_end,
      $distance_to_transcript,
      $codon_allele_string,
      $pep_allele_string,
      $allele->hgvs_genomic,
      $hgvs_transcript,
      $hgvs_protein,
      $polyphen_prediction,
      $polyphen_score,
      $sift_prediction,
      $sift_score,
      $display
    ];
  }

  return \@return;
}

# stores the data generated by _get_write_data
sub _store_write_data {
  my ($self, $write_data, $no_delay) = @_;

  $self->_generic_get_sth('', $no_delay)->execute(@$_) for @$write_data;
}

# gets a cached statement handle
sub _generic_get_sth {
  my ($self, $prefix, $no_delay) = @_;

  $prefix ||= '';

  my $key = '_'.$prefix.'store_sth';

  if(!exists($self->{$key})) {
    my $method = '_'.lc($prefix).'write_columns'; 
    my @write_columns = $self->$method;
    my $write_columns_str = join(",\n", @write_columns);
    my $values_str = '?'.(',?' x ((scalar @write_columns) - 1));
    
    my $delayed = $no_delay ? '' : 'DELAYED';
    my $table = $prefix.'transcript_variation';

    $self->{$key} ||= $self->dbc->db_handle->prepare_cached(qq{
      INSERT $delayed INTO $table ($write_columns_str) VALUES ($values_str)
    });
  }

  return $self->{$key};
}

# creates an arrayref of arrayrefs
# each one corresponding to a row of data in the MTMP table
# requires arrayref of row arrayrefs as output from _get_write_data
sub _get_mtmp_write_data_from_tv_write_data {
  my ($self, $all_data) = @_;

  my (@return, %hash);
  my @write_columns = $self->_write_columns;
  my @mtmp_write_columns = $self->_mtmp_write_columns;

  foreach my $data(@$all_data) {
    for my $i(0..$#write_columns) {
      $hash{$write_columns[$i]} = $data->[$i];
    }

    # this is the key difference between transcript_variation and MTMP_
    # MTMP_ has one row per consequence, whereas transcript_variation can
    # have multiple cons separated by "," in a single row
    foreach my $term(split(',', $hash{consequence_types})) {
      my %copy = %hash;
      $copy{consequence_types} = $term;

      push @return, [map {$copy{$_}} @mtmp_write_columns];
    }
  }

  return \@return;
}

# stores the data generated by _get_mtmp_write_data_from_tv_write_data
sub _store_mtmp_write_data {
  my ($self, $write_data, $no_delay) = @_;

  $self->_generic_get_sth('MTMP_', $no_delay)->execute(@$_) for @$write_data;
}

1;
