
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


# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor

=head1 SYNOPSIS

  $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_registry_from_db(-host => 'ensembldb.ensembl.org', -user => 'anonymous');

  my $vfa = $registry->get_adaptor('human', 'variation', 'variationfeature');
  my $sa = $registry->get_adaptor('human', 'core', 'slice');
  my $va = $registry->get_adaptor('human', 'variation', 'variation');

  # Get a VariationFeature by its internal identifier
  my $vf = $va->fetch_by_dbID(145);

  # Include the variations that have been flagged as failed
  $vfa->db->include_failed_variations(1);

  # Get all VariationFeatures in a region
  my $slice = $sa->fetch_by_region('chromosome', 'X', 1e6, 2e6);
  foreach my $vf ( @{ $vfa->fetch_all_by_Slice($slice) } ) {
    print $vf->start(), '-', $vf->end(), ' ', $vf->allele_string(), "\n";
  }

  # Fetch all genome hits for a particular variation
  my $v = $va->fetch_by_name('rs56');

  foreach my $vf ( @{ $vfa->fetch_all_by_Variation($v) } ) {
    print $vf->seq_region_name(), $vf->seq_region_start(), '-',
          $vf->seq_region_end(), "\n";
  }
  
=head1 DESCRIPTION

This adaptor provides database connectivity for VariationFeature objects.
Genomic locations of variations can be obtained from the database using this
adaptor. See the base class BaseFeatureAdaptor for more information.
By default, the 'fetch_all_by_...'-methods will not return variations
that have been flagged as failed in the Ensembl QC. This behaviour can be modified
by setting the include_failed_variations flag in Bio::EnsEMBL::Variation::DBSQL::DBAdaptor.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;

use Bio::EnsEMBL::Variation::Allele;
use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);
use Bio::EnsEMBL::Utils::Iterator;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Variation::Utils::Constants qw(%OVERLAP_CONSEQUENCES);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_hgvs_alleles trim_sequences);
use Bio::SeqUtils;
use Scalar::Util qw(looks_like_number);

our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor', 'Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor');
our $MAX_VARIATION_SET_ID = 64;
our $DEBUG =0;


## Used for the itterator function
my $DEFAULT_ITERATOR_CACHE_SIZE = 10_000;


sub store {
    my ($self, $vf) = @_;
    
    my $dbh = $self->dbc->db_handle;
    
    # look up source_id
    if(!defined($vf->{_source_id})) {
        my $sth = $dbh->prepare(q{
            SELECT source_id FROM source WHERE name = ?
        });
        $sth->execute($vf->{source}->name);
        
        my $source_id;
        $sth->bind_columns(\$source_id);
        $sth->fetch();
        $sth->finish(); 
        $vf->{_source_id} = $source_id;
    }
    if( defined $vf->{evidence}){
	## store these by attrib id to allow differnt values in different species
	my $aa = $self->db->get_AttributeAdaptor;

	foreach my $ev_term( @{$vf->{evidence}} ){

	    my $ev_class_id = $aa->attrib_id_for_type_value('evidence',$ev_term);
	    push @{$vf->{evidence_attribs}},  $ev_class_id;
	}
    }
    throw("No source ID found for source name ", $vf->{source}->name)
        unless defined($vf->{_source_id});
    
    my $sth = $dbh->prepare(q{
        INSERT INTO variation_feature (
            seq_region_id,
            seq_region_start,
            seq_region_end,
            seq_region_strand,
            variation_id,
            allele_string,
            variation_name,
            map_weight,
            flags,
            source_id,
            consequence_types,
            variation_set_id,
            class_attrib_id,
            somatic,
            minor_allele,
            minor_allele_freq,
            minor_allele_count,
            alignment_quality,
            evidence_attribs
        ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
    });
    
    $sth->execute(
        defined($vf->{seq_region_id}) ? $vf->{seq_region_id} : $vf->slice->get_seq_region_id,
        $vf->{slice} ? $vf->seq_region_start : $vf->{start},
        $vf->{slice} ? $vf->seq_region_end : $vf->{end},
        $vf->strand,
        $vf->variation ? $vf->variation->dbID : $vf->{_variation_id},
        $vf->allele_string,
        $vf->variation_name,
        $vf->map_weight || 1,
        $vf->{flags},
        $vf->{source} ? $vf->{source}->dbID : $vf->{_source_id},
        $vf->{slice} ? (join ",", @{$vf->consequence_type('SO')}) : 'intergenic_variant',
        $vf->{variation_set_id} || '',
        $vf->{class_attrib_id} || ( $vf->{class_SO_term} && $self->db->get_AttributeAdaptor->attrib_id_for_type_value('SO_term', $vf->{class_SO_term}) ) || 18,
        $vf->is_somatic,
        $vf->minor_allele,
        $vf->minor_allele_frequency,
        $vf->minor_allele_count,
        $vf->flank_match,
        $vf->{evidence_attribs} ? (join ",", @{$vf->{evidence_attribs}}) : undef
    );
    
    $sth->finish;
    
    # get dbID
	my $dbID = $dbh->last_insert_id(undef, undef, 'variation_feature', 'variation_feature_id');
    $vf->{dbID}    = $dbID;
    $vf->{adaptor} = $self;
}


=head2 update

  Arg [1]    : Bio::EnsEMBL:VariationFeature
  Example    : $variation_feature_adaptor->update( $variation_feature )};
  Description: Updates the given $variation_feature in the database,
               storing ancy changes that may have been made.
  Returntype : Bool
  Exceptions : 
  Caller     : 
  Status     : Experimental

=cut

sub update {
    my ($self, $vf) = @_;

    my $dbh = $self->dbc->db_handle;

    # look up source_id
    if(!defined($vf->{_source_id})) {
        my $sth = $dbh->prepare(q{
            SELECT source_id FROM source WHERE name = ?
        });
        $sth->execute($vf->{source}->name);

        my $source_id;
        $sth->bind_columns(\$source_id);
        $sth->fetch();
        $sth->finish();
        $vf->{_source_id} = $source_id;
    }

    throw("No source ID found for source name ", $vf->{source}->name)
        unless defined($vf->{_source_id});

    my $sth = $dbh->prepare(q{
        UPDATE variation_feature SET
            seq_region_id = ?,
            seq_region_start = ?,
            seq_region_end = ?,
            seq_region_strand = ?,
            variation_id = ?,
            allele_string = ?,
            variation_name = ?,
            map_weight = ?,
            flags = ?,
            source_id = ?,
            consequence_types = ?,
            variation_set_id = ?,
            class_attrib_id = ?,
            somatic = ?,
            minor_allele = ?,
            minor_allele_freq = ?,
            minor_allele_count = ?,
            alignment_quality = ?
         WHERE variation_feature_id = ?
    });
    
    $sth->execute(
        defined($vf->{seq_region_id}) ?
          $vf->{seq_region_id} : $vf->slice->get_seq_region_id,
        $vf->{slice} ? $vf->seq_region_start : $vf->{start},
        $vf->{slice} ? $vf->seq_region_end   : $vf->{end},
        $vf->strand,
        $vf->variation ? $vf->variation->dbID : $vf->{_variation_id},
        $vf->allele_string,
        $vf->variation_name,
        $vf->map_weight || 1,
        $vf->{flags},
        $vf->{source} ? $vf->{source}->dbID : $vf->{_source_id},
        $vf->{slice} ?
          (join ",", @{$vf->consequence_type('SO')}) : 'intergenic_variant',
        $vf->{variation_set_id} || '',
        $vf->{class_attrib_id} ||
          $vf->adaptor->db->get_AttributeAdaptor->
            attrib_id_for_type_value('SO_term', $vf->{class_SO_term}) || 18,
        $vf->is_somatic,
        $vf->minor_allele,
        $vf->minor_allele_frequency,
        $vf->minor_allele_count,
        $vf->flank_match,
        $vf->dbID
    );

    $sth->finish;
}



=head2 fetch_all

  Description: Returns a listref of all germline variation features
  Returntype : listref of VariationFeatures
  Status     : Stable

=cut

sub fetch_all {
    my $self = shift;
    my $constraint = 'vf.somatic = 0';
    return $self->generic_fetch($constraint);
}

=head2 fetch_all_somatic

  Description: Returns a listref of all somatic variation features
  Returntype : listref of VariationFeatures
  Status     : Stable

=cut

sub fetch_all_somatic {
    my $self = shift;
    my $constraint = 'vf.somatic = 1';
    return $self->generic_fetch($constraint);
}

=head2 fetch_all_by_Slice_constraint

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Arg [2]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Description: Returns a listref of germline variation features created 
               from the database which are on the Slice defined by $slice 
               and fulfill the SQL constraint defined by $constraint.
  Returntype : listref of VariationFeatures
  Exceptions : thrown if $slice is not defined
  Caller     : Bio::EnsEMBL::Slice
  Status     : Stable

=cut

sub fetch_all_by_Slice_constraint {
    my ($self, $slice, $constraint, $no_cons) = @_;
    
    # by default, filter out somatic mutations
    my $somatic_constraint = 'vf.somatic = 0';

    if ($constraint) {
        $constraint .= " AND $somatic_constraint";
    }
    else {
        $constraint = $somatic_constraint;
    }
    
    my $use_vcf = $self->db->use_vcf();
    my @vfs;
    
    if($use_vcf) { # && !$constraint) {
      push @vfs,
        map {@{$_->get_all_VariationFeatures_by_Slice($slice, $no_cons)}}
        grep {$_->use_as_source}
        @{$self->db->get_VCFCollectionAdaptor->fetch_all() || []};
    }
    if($use_vcf <= 1) {
      push @vfs, @{$self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint)};
    }
        
    return \@vfs;
}

sub fetch_all_by_Slice_constraint_with_Variations {
    my $self = shift;
    $self->{_get_variations} = 1;
    my $vfs = $self->fetch_all_by_Slice_constraint(@_);
    $self->{_get_variations} = 0;
    return $vfs;
}

sub fetch_all_by_Slice_constraint_with_TranscriptVariations {
    my $self = shift;
    $self->{_get_transcript_variations} = 1;
    my $vfs = $self->fetch_all_by_Slice_constraint(@_);
    $self->{_get_transcript_variations} = 0;
    return $vfs;
}

sub fetch_all_somatic_by_Slice_constraint_with_TranscriptVariations {
    my $self = shift;
    $self->{_get_transcript_variations} = 1;
    my $vfs = $self->fetch_all_somatic_by_Slice_constraint(@_);
    $self->{_get_transcript_variations} = 0;
    return $vfs;
}

=head2 fetch_all_somatic_by_Slice_constraint

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Arg [2]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Description: Returns a listref of somatic variation features created 
               from the database which are on the Slice defined by $slice 
               and fulfill the SQL constraint defined by $constraint.
  Returntype : listref of VariationFeatures
  Exceptions : thrown if $slice is not defined
  Caller     : Bio::EnsEMBL::Slice
  Status     : Stable

=cut

sub fetch_all_somatic_by_Slice_constraint {
    my ($self, $slice, $constraint) = @_;
    
    my $somatic_constraint = 'vf.somatic = 1';
    
    if ($constraint) {
        $constraint .= " AND $somatic_constraint";
    }
    else {
        $constraint = $somatic_constraint;
    }
    
    return $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);
}

=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Example    : my $vfs = $vfa->fetch_all_by_Slice($slice);
  Description: Retrieves all variation features on the given Slice.
               NOTE: only germline variations will be returned, if you want 
               somatic mutations use the fetch_all_somatic_by_Slice method.
  Returntype : listref of Bio::EnsEMBL::VariationFeatures
  Exceptions : none
  Caller     : Bio::EnsEMBL::Slice
  Status     : Stable

=cut

sub fetch_all_by_Slice {
  my ($self, $slice, $no_cons) = @_;
  return $self->fetch_all_by_Slice_constraint($slice, '', $no_cons);
}

=head2 fetch_all_somatic_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Example    : my $svfs = $vfa->fetch_all_somatic_by_Slice($slice);
  Description: Retrieves a list of variation features representing somatic mutations on the given Slice.
  Returntype : listref of Bio::EnsEMBL::VariationFeatures
  Exceptions : none
  Caller     : Bio::EnsEMBL::Slice
  Status     : Stable

=cut

sub fetch_all_somatic_by_Slice {
  my ($self, $slice, $no_cons) = @_;
  return $self->fetch_all_somatic_by_Slice_constraint($slice, '', $no_cons);
}

=head2 fetch_all_somatic_by_Slice_Source
  Arg [1]    : Bio::EnsEMBL::Slice $slice the slice from which to obtain features
  Arg [2]    : Bio::EnsEMBL::Variation::Source $source only return somatic mutations for the given source
  Example    : my $vfs = $vfa->fetch_all_somatic_by_Slice_Source($slice, $source);
  Description: Retrieves a list of variation features representing somatic mutations on the given Slice for a given source.
  Returntype : listref of Bio::EnsEMBL::VariationFeatures
  Exceptions : throw on incorrect argument
  Caller     : Bio::EnsEMBL::Slice
  Status     : Stable
=cut

sub fetch_all_somatic_by_Slice_Source {
  my ($self, $slice, $source) = @_;

  if (!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Bio::EnsEMBL::Slice arg expected');
  }

  if (!ref($source) || !$source->isa('Bio::EnsEMBL::Variation::Source')) {
    throw('Bio::EnsEMBL::Variation::Source arg expected');
  }

  if (!$source->dbID()) {
    warning("Source does not have dbID, cannot retrieve s variants");
    return [];
  }

  my $constraint = 'vf.source_id = ' . $source->dbID;

  return $self->fetch_all_somatic_by_Slice_constraint($slice, $constraint);
}

=head2 fetch_all_by_Variation

  Arg [1]    : Bio::EnsEMBL:Variation::Variation $var
  Example    : my @vfs = @{$vfa->fetch_all_by_Variation($var)};
  Description: Retrieves all variation features for a given variation.  Most
               variations should only hit the genome once and only a return
               a single variation feature.
  Returntype : reference to list Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut


sub fetch_all_by_Variation {
  my $self = shift;
  my $var  = shift;

  if(!ref($var) || !$var->isa('Bio::EnsEMBL::Variation::Variation')) {
    throw('Bio::EnsEMBL::Variation::Variation arg expected');
  }

  if(!defined($var->dbID())) {
    throw("Variation arg must have defined dbID");
  }

  return $self->generic_fetch("vf.variation_id = ".$var->dbID());
}

=head2 fetch_all_genotyped_by_Slice

  Arg [1]    : Bio::EnsEMBL:Variation::Slice $slice
  Example    : my @vfs = @{$vfa->fetch_all_genotyped_by_Slice($slice)};
  Description: Retrieves all variation features that have been gentoyped for a given slice.
               Most variations should only hit the genome once and only a return
               a single variation feature.
  Returntype : reference to list Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_genotyped_by_Slice{
    my $self = shift;
    my $slice = shift;

    my $constraint = "vf.flags & 1 AND vf.somatic = 0";
    #call the method fetch_all_by_Slice_constraint with the genotyped constraint
    return $self->fetch_all_by_Slice_constraint($slice,$constraint);
}

sub _internal_fetch_all_with_phenotype_by_Slice{

	my $self = shift;
	my $slice = shift;
	my $v_source = shift;
	my $p_source = shift;
	my $phenotype = shift;
	my $constraint = shift;

	if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
		throw('Bio::EnsEMBL::Slice arg expected');
	}
	
  my $extra_sql = '';
  my $extra_sql_in = '';
  my $extra_table = '';
  
  if(defined $v_source) {
    $extra_sql .= qq{ AND s.name = '$v_source' };
  }
    
  if(defined $p_source) {
    $extra_sql_in .= qq{ AND pf.source_id = ps.source_id AND ps.name = '$p_source' };
    $extra_table .= qq{, source ps};
  }
    
  if(defined $phenotype) {
    if($phenotype =~ /^[0-9]+$/) {
      $extra_sql_in .= qq{ AND pf.phenotype_id = $phenotype };
    }
    else {
      $extra_sql_in .= qq{ AND pf.phenotype_id = p.phenotype_id 
                           AND (p.name = '$phenotype' OR p.description LIKE '%$phenotype%') 
                         };
      $extra_table .= qq{, phenotype p};
    }
  }
  
  if ($constraint) {
    $extra_sql .= qq{ AND $constraint };
  }
    
  # Add the constraint for failed variations
  $extra_sql .= " AND vf.display = 1 " unless $self->db->include_failed_variations();
    
  my $cols = join ",", $self->_columns();
    
  my $sth = $self->prepare(qq{
    SELECT DISTINCT $cols
    FROM (variation_feature vf, phenotype_feature pf,
        source s $extra_table) # need to link twice to source
    WHERE vf.seq_region_id = pf.seq_region_id
    AND vf.seq_region_start = pf.seq_region_start
    AND vf.seq_region_end = pf.seq_region_end
    AND vf.variation_name = pf.object_id
    $extra_sql_in
    AND vf.source_id = s.source_id 
		$extra_sql
    AND vf.seq_region_id = ?
    AND vf.seq_region_end >= ?
    AND vf.seq_region_start <= ?
    ORDER BY vf.variation_feature_id
  });
  
  $sth->execute($slice->get_seq_region_id, $slice->start, $slice->end);
    
  return $self->_objs_from_sth($sth, undef, $slice);
}

=head2 fetch_all_with_maf_by_Slice

  Arg [1]    : Bio::EnsEMBL:Variation::Slice $slice
  Arg [2]    : float $threshold_frequency [optional]
  Arg [3]    : int $greater_than_threshold [optional]
  Example    : my @vfs = @{$vfa->fetch_all_with_maf_by_Slice($slice,$maximum_frequency)};
  Description: Retrieves all germline variation features associated with a minor
               allele frequency (MAF) for a given slice.
               The optional $threshold_frequency argument can be used to
               retrieve only variation features with a MAF lesser or equal than
               $threshold_frequency will be returned by default.
               To retrieve the MAF greater than the $threshold_frequency argument, the optional
               flag $greater_than_threshold must me populated (with the value 1 for instance).
  Returntype : reference to list Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_with_maf_by_Slice {
  my $self = shift;
  my ($slice, $freq, $greater_than_freq) = @_;

  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Bio::EnsEMBL::Slice arg expected');
  }

  my $maf_col = 'vf.minor_allele_freq';

  # Add the MAF constraint
  my $constraint = " $maf_col IS NOT NULL ";

  # Add the MAF threshold constraint
  if(defined $freq) {
   # adjust frequency if given a percentage
   $freq /= 100 if $freq > 1;
   if ($freq > 0.5) {
      throw("The maximum value of the minor allele frequency parameter should to be lesser than 0.5 (the value provided is $freq).");
   }
   my $operator = ($greater_than_freq) ? '>' : '<=';
   $constraint .= " AND $maf_col $operator $freq ";
  }

  # Add the constraint for failed variations
  $constraint .= " AND vf.display = 1 " unless $self->db->include_failed_variations();

  # Get the VariationFeatures by calling fetch_all_by_Slice_constraint
  my $vfs = $self->fetch_all_by_Slice_constraint($slice,$constraint);

  return $vfs;
}

=head2 fetch_all_with_phenotype_by_Slice

  Arg [1]    : Bio::EnsEMBL:Variation::Slice $slice
  Arg [2]    : $variation_feature_source [optional]
  Arg [3]    : $phenotype_source [optional]
  Arg [4]    : $phenotype_name [optional]
  Example    : my @vfs = @{$vfa->fetch_all_with_phenotype_by_Slice($slice)};
  Description: Retrieves all germline variation features associated with phenotypes for
               a given slice.
               The optional $variation_feature_source argument can be used to
               retrieve only variation features from a paricular source.
               The optional $phenotype source argument can be used to
               retrieve only variation features with phenotypes provided by
               a particular source.
               The optional $phenotype_name argument can
               be used to retrieve only variation features associated with
               that phenotype - this can also be a phenotype's dbID.
  Returntype : reference to list Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_with_phenotype_by_Slice {
    my $self = shift;
    my ($slice, $v_source, $p_source, $phenotype) = @_;
    my $constraint = 'vf.somatic = 0';
    return $self->_internal_fetch_all_with_phenotype_by_Slice($slice, $v_source, $p_source, $phenotype, $constraint);
}

=head2 fetch_all_somatic_with_phenotype_by_Slice

  Arg [1]    : Bio::EnsEMBL:Variation::Slice $slice
  Arg [2]    : $variation_feature_source [optional]
  Arg [3]    : $phenotype_source [optional]
  Arg [4]    : $phenotype_name [optional]
  Example    : my @vfs = @{$vfa->fetch_all_somatic_with_phenotype_by_Slice($slice)};
  Description: Retrieves all somatic variation features associated with phenotypes for
               a given slice.
               (see fetch_all_with_phenotype_by_Slice documentation for description of
               the other parameters)
  Returntype : reference to list Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_somatic_with_phenotype_by_Slice {
    my $self = shift;
    my ($slice, $v_source, $p_source, $phenotype) = @_;
    my $constraint = 'vf.somatic = 1';
    return $self->_internal_fetch_all_with_phenotype_by_Slice($slice, $v_source, $p_source, $phenotype, $constraint);
}

=head2 fetch_all_by_Slice_VariationSet

  Arg [1]    : Bio::EnsEMBL:Variation::Slice $slice
  Arg [2]    : Bio::EnsEMBL:Variation::VariationSet $set
  Example    : my @vfs =
@{$vfa->fetch_all_by_Slice_VariationSet($slice, $set)};
  Description: Retrieves all variation features in a slice that belong to a
			   given variation set.
  Returntype : reference to list Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Slice_VariationSet{

  my $self = shift;
  my $slice = shift;
  my $set = shift;
	
  #$self->{_get_variations} = 1;

  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Bio::EnsEMBL::Slice arg expected');
  }
  if(!ref($set) || !$set->isa('Bio::EnsEMBL::Variation::VariationSet')) {
    throw('Bio::EnsEMBL::Variation::VariationSet arg expected');
  }
  
  # fix for failed sets
  my $failed = $self->db->include_failed_variations;
  $self->db->include_failed_variations(1) if $failed == 0 && $set->name =~ /fail/;
  
  #Get the bitvalue for this set and its subsets
  my $bitvalue = $set->_get_bitvalue();
  
  # Add a constraint to only return VariationFeatures having the primary keys of the supplied VariationSet or its subsets in the variation_set_id column
  my $constraint = " vf.variation_set_id & $bitvalue ";
  
  #Get the VariationFeatures by calling fetch_all_by_Slice_constraint
  my $vfs = $self->fetch_all_by_Slice_constraint($slice,$constraint);
  
  # restore failed fetch flag
  $self->db->include_failed_variations($failed);

  $self->{_get_variations} = 0;

  return $vfs;
}


=head2 fetch_all_by_Slice_Population

  Arg [1]	 : Bio::EnsEMBL::Slice
  Arg [2]    : Bio::EnsEMBL::Variation::Population
  Arg [3]	 : $minimum_frequency (optional)
  Example    : $pop = $pop_adaptor->fetch_by_dbID(659);
			  $slice = $slice_adaptor->fetch_by_region("chromosome", 1, 1, 1e6);
              @vfs = @{$vf_adaptor->fetch_all_by_Slice_Population($pop,$slice)};
  Description: Retrieves all variation features in a slice which are stored for
			   a specified population. If $minimum_frequency is supplied, only
			   variations with a minor allele frequency (MAF) greater than
			   $minimum_frequency will be returned.
  Returntype : listref of Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_Slice_Population {
  my $self = shift;
 
  my $slice = shift;
  
  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Bio::EnsEMBL::Slice arg expected');
  }
  
  my $pop = shift;
  
  if(!ref($pop) || !$pop->isa('Bio::EnsEMBL::Variation::Population')) {
	throw('Bio::EnsEMBL::Variation::Population arg expected');
  }
  
  # default to 5% frequency
  my $freq = shift;
  my $extra_sql = '';
  
  if(defined $freq) {
	
	# adjust frequency if given a percentage
	$freq /= 100 if $freq > 1;
	$extra_sql = qq{ AND (IF(a.frequency > 0.5, 1-a.frequency, a.frequency) > $freq) }
  }
  
  # Add the constraint for failed variations
  $extra_sql .= " AND vf.display = 1 " unless $self->db->include_failed_variations();

  my $cols = join ",", $self->_columns();
  
  my $sth = $self->prepare(qq{
	SELECT DISTINCT $cols
	FROM (variation_feature vf, source s, allele a)
	WHERE vf.source_id = s.source_id
	AND vf.variation_id = a.variation_id
	AND a.population_id = ?
	$extra_sql
	AND vf.seq_region_id = ?
	AND vf.seq_region_end >= ?
	AND vf.seq_region_start <= ?
    ORDER BY a.variation_id
  });
  
  $sth->execute($pop->dbID, $slice->get_seq_region_id, $slice->start, $slice->end);

##### ALTERNATIVE USING SUB-QUERY - WORKS FOR MULTIALLELIC BUT TAKES FOREVER!!!
#  my $cols = join ",", $self->_columns();
#  
#  my $sth = $self->prepare(qq{
#	SELECT $cols
#	FROM variation_feature vf, source s, allele a_outer
#	WHERE
#	  vf.source_id = s.source_id AND
#	  vf.variation_id = a_outer.variation_id AND
#	  a_outer.population_id = ? AND
#	  a_outer.frequency > ? AND
#	  a_outer.frequency < 0.5 AND
#	  vf.seq_region_id = ? AND
#	  vf.seq_region_end >= ? AND
#	  vf.seq_region_start <= ? AND
#	  NOT EXISTS (
#		SELECT
#		  *
#		FROM
#		  allele a_inner
#		WHERE
#		  a_inner.variation_id = a_outer.variation_id AND
#		  a_inner.population_id = a_outer.population_id AND
#		  a_inner.frequency < a_outer.frequency
#	  );
#  });
  
  #$sth->execute($pop->dbID, $freq, $slice->get_seq_region_id, $slice->start, $slice->end);
  
  return $self->_objs_from_sth($sth, undef, $slice);
}

sub _internal_fetch_all_with_phenotype {
    
  my ($self, $v_source, $p_source, $phenotype, $constraint) = @_;
    
  my $extra_sql = '';
  my $extra_table = '';
    
  if(defined $v_source) {
    $extra_sql .= qq{ AND s.name = '$v_source' };
  }
    
  if(defined $p_source) {
    $extra_sql .= qq{ AND pf.source_id = ps.source_id AND ps.name = '$p_source' };
    $extra_table .= qq{, source ps};
  }
    
  if(defined $phenotype) {
    if($phenotype =~ /^[0-9]+$/) {
      $extra_sql .= qq{ AND pf.phenotype_id = $phenotype };
    }
    else {
      $extra_sql .= qq{ AND pf.phenotype_id = p.phenotype_id 
					              AND (p.name = '$phenotype' OR p.description LIKE '%$phenotype%')
                      };
			$extra_table .= qq{, phenotype p};
    }
  }
    
  if ($constraint) {
    $extra_sql .= qq{ AND $constraint };
  }
    
  # Add the constraint for failed variations
  $extra_sql .= " AND vf.display = 1 " unless $self->db->include_failed_variations();
    
  my $cols = join ",", $self->_columns();
    
  my $sth = $self->prepare(qq{
        SELECT DISTINCT $cols
        FROM (variation_feature vf, phenotype_feature pf,
        source s $extra_table) # need to link twice to source
        WHERE vf.source_id = s.source_id
        AND vf.seq_region_id = pf.seq_region_id
				AND vf.seq_region_start = pf.seq_region_start
				AND vf.seq_region_end = pf.seq_region_end
				AND vf.variation_name = pf.object_id
        $extra_sql
        ORDER BY vf.variation_feature_id
  });
    
  $sth->execute;
    
  return $self->_objs_from_sth($sth);
}

=head2 fetch_all_with_phenotype

  Arg [1]    : $variation_feature_source [optional]
  Arg [2]    : $phenotype_source [optional]
  Arg [3]    : $phenotype_name [optional]
  Example    : my @vfs = @{$vfa->fetch_all_with_phenotype('EGA', undef, 123)};
  Description: Retrieves all germline variation features associated with the given phenotype
  Returntype : reference to list Bio::EnsEMBL::Variation::VariationFeature
  Caller     : webcode
  Status     : Experimental

=cut

sub fetch_all_with_phenotype {
    
    my ($self, $v_source, $p_source, $phenotype, $constraint) = @_;
    
    my $somatic_constraint = 'vf.somatic = 0';
    
    if ($constraint) {
        $constraint .= " AND $somatic_constraint";
    }
    else {
        $constraint = $somatic_constraint;
    }
    
    return $self->_internal_fetch_all_with_phenotype($v_source, $p_source, $phenotype, $constraint);
}

=head2 fetch_all_somatic_with_phenotype

  Arg [1]    : $variation_feature_source [optional]
  Arg [2]    : $phenotype_source [optional]
  Arg [3]    : $phenotype_name [optional]
  Example    : my @vfs = @{$vfa->fetch_all_somatic_with_phenotype('COSMIC', undef, 807)};
  Description: Retrieves all somatic variation features associated with the given phenotype
  Returntype : reference to list Bio::EnsEMBL::Variation::VariationFeature
  Caller     : webcode
  Status     : Experimental

=cut

sub fetch_all_somatic_with_phenotype {
    
    my ($self, $v_source, $p_source, $phenotype, $constraint) = @_;
    
    my $somatic_constraint = 'vf.somatic = 1';
    
    if ($constraint) {
        $constraint .= " AND $somatic_constraint";
    }
    else {
        $constraint = $somatic_constraint;
    }
    
    return $self->_internal_fetch_all_with_phenotype($v_source, $p_source, $phenotype, $constraint);
}

=head2 fetch_all_by_Slice_SO_terms

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : listref of SO terms
  Description: Fetch all germline VariationFeatures on the given slice with
               consequences with given SO terms
  Returntype : listref of Bio::EnsEMBL::Variation::VariationFeatures
  Status     : At risk

=cut

sub fetch_all_by_Slice_SO_terms {
  my ($self, $slice, $terms, $without_children, $included_so) = @_;

  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Bio::EnsEMBL::Slice arg expected');
  }
  
  if(!defined($terms) || scalar @$terms == 0) {
	return $self->fetch_all_by_Slice($slice);
  }

  my $constraint = $self->_get_consequence_constraint($terms, $without_children, $included_so);
  if (!$constraint) {
    return [];
  }
  
  my $vfs = $self->fetch_all_by_Slice_constraint($slice,$constraint);

  return $vfs;
}

=head2 fetch_all_somatic_by_Slice_SO_terms

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : listref of SO terms
  Description: Fetch all somatic VariationFeatures on the given slice with
               consequences with given SO terms
  Returntype : listref of Bio::EnsEMBL::Variation::VariationFeatures
  Status     : At risk

=cut

sub fetch_all_somatic_by_Slice_SO_terms {
  my ($self, $slice, $terms, $without_children, $included_so) = @_;

  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Bio::EnsEMBL::Slice arg expected');
  }
  
  if(!defined($terms) || scalar @$terms == 0) {
	return $self->fetch_all_somatic_by_Slice($slice);
  }

  my $constraint = $self->_get_consequence_constraint($terms, $without_children, $included_so);
  if (!$constraint) {
    return [];
  }
  
  my $vfs = $self->fetch_all_somatic_by_Slice_constraint($slice,$constraint);

  return $vfs;
}

=head2 fetch_all_by_Slice_VariationSet_SO_terms

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : Bio::EnsEMBL::Variation::VariationSet
  Arg [3]    : Optional listref of SO terms
  Description: Fetch all VariationFeatures in the given set, on the given
               slice, optionally with consequences with given SO terms 
  Returntype : listref of Bio::EnsEMBL::Variation::VariationFeatures
  Status     : At risk

=cut

sub fetch_all_by_Slice_VariationSet_SO_terms{

my ($self, $slice, $set, $terms, $without_children, $included_so) = @_;

  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Bio::EnsEMBL::Slice arg expected');
  }

  if( !ref($set) || !$set->isa('Bio::EnsEMBL::Variation::VariationSet')) {
    throw('Bio::EnsEMBL::Variation::VariationSet arg expected');
  }

  if(!defined($terms) || scalar @$terms == 0 ) {
    return $self->fetch_all_by_Slice_VariationSet($slice, $set);
  }

  my $constraint = $self->_get_consequence_constraint($terms, $without_children, $included_so);

  # Get the bitvalue for this set and its subsets
  my $bitvalue = $set->_get_bitvalue();  
  $constraint .= " and vf.variation_set_id & $bitvalue ";
  
  my $vfs = $self->fetch_all_by_Slice_constraint($slice,$constraint);

  return $vfs;
}


# call to method in BaseAdaptor
sub _get_consequence_constraint {
	my $self = shift;
	return $self->SUPER::_get_consequence_constraint('variation_feature', @_);
}

=head2 fetch_Iterator

  Arg [1]    : int $cache_size (optional)
  Example    : $var_iterator = $varf_adaptor->fetch_Iterator;
  Description: returns an iterator over all variation features in the database
  Returntype : Bio::EnsEMBL::Utils::Iterator
  Exceptions : none
  Caller     : general
  Status     : Experimental

=cut

sub fetch_Iterator {
    my ($self, $cache_size) = @_;
    return $self->_generic_fetch_Iterator($cache_size);
}


sub _generic_fetch_Iterator {

    my ($self, $cache_size, $constraint) = @_;

    # prepare and execute a query to fetch all dbIDs
    my $sth = $self->prepare(qq{
        SELECT      variation_feature_id
        FROM        variation_feature
    });
    $sth->execute;

    my $varf_id;
    $sth->bind_columns(\$varf_id);

    # we probably can't fit all of these into memory at once though,
    # so create an iterator that fetches $cache_size dbIDs from the
    # statement handle at a time and then fetches these objects,
    # storing them in a cache. We then return variation feature
    # objects from this cache one by one, before filling it again if
    # necessary

    $cache_size ||= $DEFAULT_ITERATOR_CACHE_SIZE;

    my @cache;
    my $items_to_fetch = 1;

    return Bio::EnsEMBL::Utils::Iterator->new(sub{

        if (@cache == 0 && $items_to_fetch) {
            # Our cache is empty, and there are still items to fetch, so
            # fetch the next chunk of dbIDs and create objects from them

            my @dbIDs;
            my $item_count = 0;

            while( $sth->fetch ) {
                push @dbIDs, $varf_id;

                if (++$item_count == $cache_size) {
                    # we have fetched a cache's worth of dbIDs, so
                    # flag that there are still items to fetch and
                    # last out of the loop
                    $items_to_fetch = 1;
                    last;
                }

                # If this is the last row, this flag will be 0 outside
                # the loop
                $items_to_fetch = 0;
            }

            $sth->finish unless $items_to_fetch;

            # This looks tricky...
            @cache = @{ $self->fetch_all_by_dbID_list(\@dbIDs) } if @dbIDs;
        }

        return shift @cache;
    });
}



# method used by import VCF script
sub _fetch_all_by_coords {
    my ($self, $seq_region_id, $start, $end, $somatic) = @_;
	
	$somatic ||= 0;
    
    return $self->generic_fetch(qq{
        vf.seq_region_id = $seq_region_id AND
        vf.seq_region_start = $start AND
        vf.seq_region_end = $end AND
		vf.somatic = $somatic
    });
}

# method used by superclass to construct SQL
sub _tables {
	my $self = shift;
	
  my $dbc = $self->dbc;
  my @tables;
  if ($dbc->driver eq 'SQLite') {
	  @tables = (
		  [ 'variation_feature', 'vf'],
		  [ 'source', 's']
	  );
  } else {
    @tables = (
		  [ 'variation_feature', 'vf', 'IGNORE INDEX(consequence_type_idx)'],
		  [ 'source', 's']
	  );
  }

  return @tables;
}


sub _default_where_clause {
  my $self = shift;
  
  my $clause = 'vf.source_id = s.source_id';

  $clause .= ' AND vf.display = 1 ' unless $self->db->include_failed_variations();

  return $clause;
}

sub _columns {
  return qw( vf.variation_feature_id vf.seq_region_id vf.seq_region_start
             vf.seq_region_end vf.seq_region_strand vf.variation_id
             vf.allele_string vf.variation_name vf.map_weight vf.source_id vf.somatic 
             vf.consequence_types vf.class_attrib_id
             vf.minor_allele vf.minor_allele_freq vf.minor_allele_count vf.alignment_quality 
             vf.evidence_attribs vf.clinical_significance vf.display s.name);
}

sub _objs_from_sth {
    my ($self, $sth, $mapper, $dest_slice) = @_;
	
    #warn $sth->sql;

    # 
    # This code is ugly because an attempt has been made to remove as many
    # function calls as possible for speed purposes.  Thus many caches and
    # a fair bit of gymnastics is used.
    #

    my $sa = $self->db()->dnadb()->get_SliceAdaptor();

    my $aa = $self->db->get_AttributeAdaptor;

    my @features;
    my %slice_hash;
    my %sr_name_hash;
    my %sr_cs_hash;

    my ($variation_feature_id, $seq_region_id, $seq_region_start,
      $seq_region_end, $seq_region_strand, $variation_id,
      $allele_string, $variation_name, $map_weight, $source_id,
      $is_somatic, $consequence_types, $class_attrib_id,
      $minor_allele, $minor_allele_freq, $minor_allele_count, $last_vf_id,
      $alignment_quality,$evidence_attribs,$clin_sig,$display, $source_name );

    $sth->bind_columns(\$variation_feature_id, \$seq_region_id,
                     \$seq_region_start, \$seq_region_end, \$seq_region_strand,
                     \$variation_id, \$allele_string, \$variation_name,
                     \$map_weight, \$source_id, \$is_somatic,
                     \$consequence_types, \$class_attrib_id,
                     \$minor_allele, \$minor_allele_freq, \$minor_allele_count,
                     \$alignment_quality, \$evidence_attribs, \$clin_sig, \$display, \$source_name);

    my $asm_cs;
    my $cmp_cs;
    my $asm_cs_vers;
    my $asm_cs_name;
    my $cmp_cs_vers;
    my $cmp_cs_name;
    
    if($mapper) {
        $asm_cs = $mapper->assembled_CoordSystem();
        $cmp_cs = $mapper->component_CoordSystem();
        $asm_cs_name = $asm_cs->name();
        $asm_cs_vers = $asm_cs->version();
        $cmp_cs_name = $cmp_cs->name();
        $cmp_cs_vers = $cmp_cs->version();
    }

    my $dest_slice_start;
    my $dest_slice_end;
    my $dest_slice_strand;
    my $dest_slice_length;
    
    if($dest_slice) {
        $dest_slice_start  = $dest_slice->start();
        $dest_slice_end    = $dest_slice->end();
        $dest_slice_strand = $dest_slice->strand();
        $dest_slice_length = $dest_slice->length();
    }

    my $finished = 0;
    
    my $iterator = Bio::EnsEMBL::Utils::Iterator->new(sub{    
        
        return undef if $finished;

        FEATURE: while( $sth->fetch ) {
        
            # Skip if we are getting multiple rows because of the left join to failed variation
            next if (defined($last_vf_id) && $last_vf_id == $variation_feature_id);
            $last_vf_id = $variation_feature_id;
    
            #get the slice object
            my $slice = $slice_hash{"ID:".$seq_region_id};
            if(!$slice) {
                $slice = $sa->fetch_by_seq_region_id($seq_region_id);
                $slice_hash{"ID:".$seq_region_id} = $slice;
                $sr_name_hash{$seq_region_id} = $slice->seq_region_name();
                $sr_cs_hash{$seq_region_id} = $slice->coord_system();
            }
            
            # remap the feature coordinates to another coord system
            # if a mapper was provided
            
            if($mapper) {
                my $sr_name = $sr_name_hash{$seq_region_id};
                my $sr_cs   = $sr_cs_hash{$seq_region_id};
        
                ($sr_name,$seq_region_start,$seq_region_end,$seq_region_strand) =
                    $mapper->fastmap($sr_name, $seq_region_start, $seq_region_end,
                                        $seq_region_strand, $sr_cs);
        
                #skip features that map to gaps or coord system boundaries
                next FEATURE if(!defined($sr_name));
        
                #get a slice in the coord system we just mapped to
                if($asm_cs == $sr_cs || ($cmp_cs != $sr_cs && $asm_cs->equals($sr_cs))) {
                    $slice = $slice_hash{"NAME:$sr_name:$cmp_cs_name:$cmp_cs_vers"} ||=
                    $sa->fetch_by_region($cmp_cs_name, $sr_name,undef, undef, undef,
                                        $cmp_cs_vers);
                } else {
                    $slice = $slice_hash{"NAME:$sr_name:$asm_cs_name:$asm_cs_vers"} ||=
                    $sa->fetch_by_region($asm_cs_name, $sr_name, undef, undef, undef,
                                        $asm_cs_vers);
                }
            }
        
            #
            # If a destination slice was provided convert the coords
            # If the dest_slice starts at 1 and is foward strand, nothing needs doing
            #
            if($dest_slice) {
                if($dest_slice_start != 1 || $dest_slice_strand != 1) {
                    if($dest_slice_strand == 1) {
                        $seq_region_start = $seq_region_start - $dest_slice_start + 1;
                        $seq_region_end   = $seq_region_end   - $dest_slice_start + 1;
                    } else {
                        my $tmp_seq_region_start = $seq_region_start;
                        $seq_region_start = $dest_slice_end - $seq_region_end + 1;
                        $seq_region_end   = $dest_slice_end - $tmp_seq_region_start + 1;
                        #$seq_region_strand *= -1;
                    }
        
                    #throw away features off the end of the requested slice
                    if($seq_region_end < 1 || $seq_region_start > $dest_slice_length) {
                        next FEATURE;
                    }
                }
                $slice = $dest_slice;
            }

            my @evidence;
            if (defined($evidence_attribs)) {
              my $aa  = $self->db->get_AttributeAdaptor;
              my @attrib_ids = split(/,/,$evidence_attribs);
              foreach my $attrib_id (@attrib_ids){
                my $evidence_value = $aa->attrib_value_for_id($attrib_id);
                push @evidence, $evidence_value ;
              }
            }

            my @clin_sig;
            if (defined($clin_sig)) {
              @clin_sig = split(/,/,$clin_sig );
            }

            #my $overlap_consequences = $self->_variation_feature_consequences_for_set_number($consequence_types);
            
            my $overlap_consequences = [ map { $OVERLAP_CONSEQUENCES{$_} } split /,/, ($consequence_types || 'sequence_variant') ];

            $seq_region_start = $seq_region_end if $seq_region_start == $seq_region_end + 2;  ## patch problem fix
            # consequence_types
            return $self->_create_feature_fast('Bio::EnsEMBL::Variation::VariationFeature',
            #push @features, Bio::EnsEMBL::Variation::VariationFeature->new_fast(
            #if use new_fast, then do not need "-" infront of key, i.e 'start' => $seq_region_start,
        
                {'start'    => $seq_region_start,
                'end'      => $seq_region_end,
                'strand'   => $seq_region_strand,
                'slice'    => $slice,
                'allele_string' => $allele_string,
                'variation_name' => $variation_name,
                'adaptor'  => $self,
                'dbID'     => $variation_feature_id,
                'map_weight' => $map_weight,
                '_source_id'   => $source_id,
                '_source_name' => $source_name,
                'is_somatic' => $is_somatic,
                'overlap_consequences' => $overlap_consequences,
                '_variation_id' => $variation_id,
                'class_SO_term' => $aa->attrib_value_for_id($class_attrib_id),
                'minor_allele' => $minor_allele,
                'minor_allele_frequency' => $minor_allele_freq,
                'minor_allele_count' => $minor_allele_count,
                'flank_match'  => $alignment_quality,
                'evidence'     => \@evidence,
                'clinical_significance' => \@clin_sig,
                'display'   => $display
                }
            );
        }
        
        unless ($finished) {
            $sth->finish;
            $finished = 1;
        }
        
        return undef;
    });
    
    if ($self->{_iterator}) {
        return $iterator;
    }
    else {
        if ($self->{_get_variations}) {
            my $vfs = $iterator->to_arrayref;
            my @v_ids = map { $_->{_variation_id} } @$vfs;
            my $vs = $self->db->get_VariationAdaptor->fetch_all_by_dbID_list(\@v_ids);
            my %vs_by_id = map { $_->dbID => $_ } @$vs;
            #warn "Got variations";
            map { $_->variation( $vs_by_id{ $_->{_variation_id} }) } @$vfs;
            return $vfs;
        }
        if ($self->{_get_transcript_variations}) {
            my $vfs = $iterator->to_arrayref;
            return $vfs unless @$vfs;
            #warn "getting transcript variations";
            my $tvs = $self->db->get_TranscriptVariationAdaptor->fetch_all_by_VariationFeatures($vfs);
            for my $tv (@$tvs) {
                $tv->variation_feature->add_TranscriptVariation($tv);
                #$tv->variation_feature->{transcript_variations}->{$tv->transcript_stable_id} = $tv;
            }
            return $vfs;
        }
        else {
            my $vfs = $iterator->to_arrayref;
            #warn "Got ".scalar(@$vfs). "VFs";
            return $vfs;
        }
    }   
}


=head2 list_dbIDs

  Arg [1]    : none
  Example    : @feature_ids = @{$simple_feature_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all simple features in 
               the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub list_dbIDs {
  my $self = shift;
  return $self->_list_dbIDs('variation_feature');
}


=head2 get_all_synonym_sources

    Args[1]     : Bio::EnsEMBL::Variation::VariationFeature vf
    Example     : my @sources = @{$vf_adaptor->get_all_synonym_sources($vf)};
    Description : returns a list of all the sources for synonyms of this
                  VariationFeature
    ReturnType  : reference to list of strings
    Exceptions  : none
    Caller      : general
    Status      : At Risk
                : Variation database is under development.

=cut

sub get_all_synonym_sources{
    my $self = shift;
    my $vf = shift;
    my %sources;
    my @sources;

    if(!ref($vf) || !$vf->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
   throw("Bio::EnsEMBL::Variation::VariationFeature argument expected");
    }
    
    if (!defined($vf->{'_variation_id'}) && !defined($vf->{'variation'})){
  warning("Not possible to get synonym sources for the VariationFeature: you need to attach a Variation first");
  return \@sources;
    }
    #get the variation_id
    my $variation_id;
    if (defined ($vf->{'_variation_id'})){
  $variation_id = $vf->{'_variation_id'};
    }
    else{
  $variation_id = $vf->variation->dbID();
    }
    #and go to the varyation_synonym table to get the extra sources
    my $source_name;
    my $sth = $self->prepare(qq{SELECT s.name 
        FROM variation_synonym vs, source s 
        WHERE s.source_id = vs.source_id
              AND   vs.variation_id = ?
          });
    $sth->bind_param(1,$variation_id,SQL_INTEGER);
    $sth->execute();
    $sth->bind_columns(\$source_name);
    while ($sth->fetch){
  $sources{$source_name}++;
    }
    @sources = keys(%sources); 
    return \@sources;
}

=head2 new_fake

  Arg [1]    : string $species
  Example    :
  $vfa = Bio::EnsEMBL::Variation::VariationFeatureAdaptor->new_fake('human');
  Description: Creates a VariationFeatureAdaptor with no underlying database
         attached. Should be used only when getting consequence types for
         species with no variation database available.
  Returntype : Bio::EnsEMBL::Variation::VariationFeatureAdaptor
  Exceptions : throw if no species given
  Caller     : called from webcode for species where no variation database present
  Status     : Stable

=cut

sub new_fake {
  my $class = shift;
  my $species = shift;
  
  throw("No species defined") unless defined $species;
  
  my $self = bless {}, $class;
  
  $self->{'species'} = $species;
  
  return $self;
}

sub _parse_hgvs_genomic_position {

  my $description  = shift;
    
  my ($start, $start_offset, $end, $end_offset) = $description=~ m/^([\-\*]?\d+)((?:[\+\-]\d+)?)(?:_([\-\*]?\d+)((?:[\+\-]\d+)?))?/;
  ## end information needed even if same as start
  unless ($end){$end = $start;}  
  unless ($end_offset){$end_offset = $start_offset;} 

  return ($start, $end);
}


sub _parse_hgvs_transcript_position {
  ### Work out genomic coordinates from hgvs coding or non coding annotation

  ### Non-exonic notation - 
  ### P+n  => n bases after listed exonic base
  ### P-n  => n bases before listed exonic base
  ### -P+n => n bases after listed 5'UTR base
  ### *P-n => n bases before listed 3'UTR base

  my $description = shift;
  my $transcript  = shift;
            
  my ($start,$start_offset, $end, $end_offset) = $description =~ m/^([\-\*]?\d+)((?:[\+\-]\d+)?)(?:_([\-\*]?\d+)((?:[\+\-]\d+)?))?/;

 my $is_exonic = 1;
 $is_exonic = 0 if ($start_offset || $end_offset || substr($start,0,1) eq '*' || (defined $end && substr($end,0,1) eq '*') || $start < 0);

  ## go back or forward into intron
  my $start_direction = 0;
  my $end_direction   = 0;

  if($start_offset ){
  ### extract + or - for intronic positions in coding nomenclature
    if (substr($start_offset,0,1) eq '+' || substr($start_offset,0,1) eq '-'){
      $start_direction  = substr($start_offset,0,1);  
      $start_offset     = substr($start_offset,1) ;
      $start            = $start;
    }
  }
  else{  ### exonic
    $start_offset = 0 ;
  }
    
  if($end_offset ){
  ###  this is needed for long intronic events eg. ENST00000299272.5:c.98-354_98-351dupGAAA
    if (substr($end_offset,0,1) eq '+' || substr($end_offset,0,1) eq '-'){
       $end_direction   = substr($end_offset,0,1); 
       $end_offset      = substr($end_offset,1) ;
       $end             = $end
     }
  }
  else{
    $end_offset   =  $start_offset;    
  }
  ### add missing values if single-location variant - needed for refseq check later
  unless (defined $end) {
    $end          = $start;
    $end_direction= $start_direction  ;
  }
    
    ### Variant in the 3' UTR =>  convert the coordinates by setting them to be the stop codon position + the UTR offset
  if (substr($start,0,1) eq '*'){
     $start = ($transcript->cdna_coding_end() - $transcript->cdna_coding_start() + 1) + int(substr($start,1)) ;
  }
  if (substr($end,0,1) eq '*'){
     $end   = ($transcript->cdna_coding_end() - $transcript->cdna_coding_start() + 1) + int(substr($end,1));
  }
    
  ### Variant in the 5' UTR =>  convert the coordinates by setting them to be the start codon position(0) - the UTR offset
  if (substr($start,0,1) eq '-'){
     $start =  0 - int(substr($start,1)) ;
  }
  if (substr($end,0,1) eq '-'){
     $end   =  0 - int(substr($end,1));
  }
        
  # Get the TranscriptMapper to convert to genomic coords
  my $tr_mapper = $transcript->get_TranscriptMapper();    
    
  if($DEBUG ==1){print "About to convert to genomic $start $end, ccs:". $transcript->cdna_coding_start() ."\n";}
    #The mapper can only convert cDNA coordinates, but we have CDS (relative to the start codon), so we need to convert them
    my ($cds_start, $cds_end) ;
    if( defined $transcript->cdna_coding_start()){
         ($cds_start, $cds_end)  = (($start + $transcript->cdna_coding_start() - ($start > 0)),($end + $transcript->cdna_coding_start() - ($end > 0)));
    }
    else{
  #### non coding transcript
        ($cds_start, $cds_end)  = ($start, $end);
    }
    # Convert the cDNA coordinates to genomic coordinates.
    my @coords = $tr_mapper->cdna2genomic($cds_start,$cds_end);
    if($DEBUG ==1){
      print "In parser: cdna2genomic coords: ". $coords[0]->start() . "-". $coords[0]->end() . " and strand ". $coords[0]->strand()." from $cds_start,$cds_end\n";}
    
    #Throw an error if we didn't get an unambiguous coordinate back
    throw ("Unable to map the cDNA coordinates $start\-$end to genomic coordinates for Transcript " .$transcript->stable_id()) if (scalar(@coords) != 1 || !$coords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate'));
    
    my  $strand = $coords[0]->strand();    
    
    ### overwrite exonic location with genomic coordinates
    $start = $coords[0]->start(); 
    $end   = $coords[0]->end();

    #### intronic variants are described as after or before the nearest exon 
    #### - add this offset to genomic start & end positions
   if( $start_direction ne 0 ||  $end_direction ne 0 ){
      if($strand  == 1){

        if($start_direction eq "+"){ $start = $start + $start_offset; }
        if($end_direction   eq "+"){ $end   = $end   + $end_offset;   }
      
        if($start_direction eq "-"){ $start = $start - $start_offset; }
        if($end_direction   eq "-"){ $end   = $end   - $end_offset;   }
     }
      elsif($strand  == -1 ){

        ## start & end swapped for rev strand features => apply offsets
        if($end_direction   eq "+"){ $start = $start - $end_offset;}
        if($start_direction eq "+"){ $end   = $end   - $start_offset;  }

        if($end_direction   eq "-"){ $start = $start + $end_offset;}
        if($start_direction eq "-"){ $end   = $end   + $start_offset;  }
 
     }
   }

  return ($start, $end, $strand, $is_exonic);
}

=head2 fetch_by_hgvs_notation

    Arg[1]      : String $hgvs
    Example     : my $hgvs = 'LRG_8t1:c.3891A>T';
                  $vf = $vf_adaptor->fetch_by_hgvs_notation($hgvs);
    Description : Parses an HGVS notation and tries to create a VariationFeature object
    based on the notation. The object will have a Variation and Alleles attached.
    ReturnType  : Bio::EnsEMBL::Variation::VariationFeature, undef on failure
    Exceptions  : thrown on error
    Caller      : general
    Status      : Stable

=cut
   
sub fetch_by_hgvs_notation {
  my $self = shift;

  my ($hgvs, $user_slice_adaptor, $user_transcript_adaptor, $multiple_ok, $replace_ref) = rearrange([qw(
    HGVS
    SLICE_ADAPTOR
    TRANSCRIPT_ADAPTOR
    MULTIPLE_OK
    REPLACE_REF
  )], @_);

  if($DEBUG ==1){print "\nStarting fetch_by_hgvs_notation for $hgvs\n";}

 
  ########################### Check & split input ###########################

  #Split the HGVS notation into the reference, notation type and variation description
  my ($reference,$type,$description) = $hgvs =~ m/^([^\:]+)\:.*?([cgmnrp]?)\.?(.*?[\*\-0-9]+.*)$/i;

  #If any of the fields are unknown, return undef
  throw ("Could not parse the HGVS notation $hgvs") 
      unless (defined($reference) && defined($type) && defined($description));

      
  my $extra;
  if($description =~ m/\(.+\)/) {        
    ($description, $extra) = $description=~ /(.+?)(\(.+\))/;        
    throw ("Could not parse the HGVS notation $hgvs - can't interpret \'$extra\'") unless $extra eq '(p.=)';
  }
  
  # strip version number from reference
  if ($reference=~ /^ENS|^LRG_\d+/){
    $reference =~ s/\.\d+//g;
   }

  $description =~ s/\s+//;

  #######################  extract genomic coordinates and reference seq allele  #######################
  my ($start, $end, $strand, $ref_allele, $alt_allele);

  #Get a slice adaptor to enable check of supplied reference allele
  my $slice_adaptor = $user_slice_adaptor || $self->db()->dnadb()->get_SliceAdaptor();
  my $slice ;

  if($type =~ m/c|n/i) {   

    #Get the Transcript object to convert coordinates
    my $transcript_adaptor = $user_transcript_adaptor || $self->db()->dnadb()->get_TranscriptAdaptor();
    my $transcript = $transcript_adaptor->fetch_by_stable_id($reference);

    my @transcripts;
    
    #try and fetch via gene
    if(!defined($transcript)) {
      push @transcripts, @{$self->_get_gene_transcripts($transcript_adaptor, $reference, $multiple_ok)};
    }
    else {
      push @transcripts, $transcript;
    }
    
    throw ("Could not get a Transcript object for '$reference'") unless @transcripts;

    my @results;
    my %errors;

    foreach my $transcript(@transcripts) {
      my ($is_exonic, $result);

      eval {        
        ($start, $end, $strand, $is_exonic) = _parse_hgvs_transcript_position($description, $transcript) ; 

        $slice = $slice_adaptor->fetch_by_region($transcript->coord_system_name(),$transcript->seq_region_name());   
        ($ref_allele, $alt_allele) = get_hgvs_alleles($hgvs);

        $result = $self->_hgvs_from_components(
          $hgvs, $reference, $type, $description,
          $slice, $ref_allele, $alt_allele, $start, $end, $strand, $replace_ref
        )
      };

      if($@) {
        $errors{$@} = 1;
      }
      else {
        push @results, $result;
      }
    }

    throw(join("\n", keys %errors)) unless @results;

    return $multiple_ok ? \@results : $results[0];
  }

  elsif($type =~ m/g|m/i) {
    ($start, $end) = _parse_hgvs_genomic_position($description) ;  
    
    throw ("Region requested must be smaller than 5kb") unless ($end - $start < 5000);
        
    ## grab reference allele; second call after "||" allows for LRG regions to be fetched
    $slice = $slice_adaptor->fetch_by_region('chromosome', $reference ) || $slice_adaptor->fetch_by_region(undef, $reference);    
    $strand =1; ## strand should be genome strand for HGVS genomic notation
    ($ref_allele, $alt_allele) = get_hgvs_alleles($hgvs);

    my $result = $self->_hgvs_from_components(
      $hgvs, $reference, $type, $description,
      $slice, $ref_allele, $alt_allele, $start, $end, $strand, $replace_ref
    );
    return $multiple_ok ? [$result] : $result;
  }
         
  elsif($type =~ m/p/i) {
  
    #Get the Transcript object to convert coordinates
    my $transcript_adaptor = $user_transcript_adaptor || $self->db()->dnadb()->get_TranscriptAdaptor();
    my $transcript = $transcript_adaptor->fetch_by_translation_stable_id($reference);
    
    my @transcripts;

    #try and fetch via gene
    if(!defined($transcript)) {
      push @transcripts, @{$self->_get_gene_transcripts($transcript_adaptor, $reference, $multiple_ok)};
    }
    else {
      push @transcripts, $transcript;
    }
    
    throw ("Could not get a Transcript object for '$reference'") unless @transcripts;

    my @results;
    my %errors;

    foreach my $transcript(@transcripts) {
      
      my $possible_prot;
      eval {
        $possible_prot = _parse_hgvs_protein_position($description, $reference, $transcript);
        throw("Could not uniquely determine nucleotide change from $hgvs") if scalar @$possible_prot > 1 && !$multiple_ok;
        $slice = $slice_adaptor->fetch_by_region($transcript->coord_system_name(), $transcript->seq_region_name());
      };
      if($@) {
        $errors{$@} = 1;
        next;
      }

      foreach my $poss(@$possible_prot) {

        my $result;

        eval {
          $result = $self->_hgvs_from_components(
            $hgvs, $reference, $type, $description, $slice,
            @$poss, $replace_ref
          )
        };

        if($@) {
          $errors{$@} = 1;
        }
        else {
          push @results, $result;
        }
      }
    }

    throw(join("\n", keys %errors)) unless @results;

    return $multiple_ok ? \@results : $results[0];
  }

  else {
    throw ("HGVS type $type not recognised for $hgvs");
  }
}

sub _get_gene_transcripts {
  my ($self, $transcript_adaptor, $reference, $multiple_ok) = @_;

  my $gene_adaptor = $transcript_adaptor->db->get_GeneAdaptor();
  my ($gene) = grep {($_->external_name || '') eq $reference} @{$gene_adaptor->fetch_all_by_external_name($reference)};
  
  my @transcripts;

  if($gene) {
    if($multiple_ok) {
      push @transcripts, @{$gene->get_all_Transcripts};
    }
    else {
      push @transcripts, $self->_pick_likely_transcript($gene->get_all_Transcripts);
      print STDERR "Picked transcript ".$transcripts[0]->stable_id."\n" if $DEBUG==1;
    }
  }

  return \@transcripts;
}


=head2 fetch_all_possible_by_hgvs_notation

    Arg[1]      : String $hgvs
    Example     : my $hgvs = 'ENSP00000284967.6:p.Glu2Asp';
                  $vfs = $vf_adaptor->fetch_all_possible_by_hgvs_notation($hgvs);
    Description : Parses an HGVS notation and tries to create VariationFeature objects.
                  Some protein (p.) notations may reverse translate to multiple
                  possible genomic changes; this method returns all possible
                  variants representing the shortest path from reference to alternate
                  allele.
    ReturnType  : arrayref of Bio::EnsEMBL::Variation::VariationFeature
    Exceptions  : thrown on error
    Caller      : general
    Status      : Stable

=cut
   
sub fetch_all_possible_by_hgvs_notation {
  my $self = shift;

  my ($hgvs, $user_slice_adaptor, $user_transcript_adaptor, $replace_ref) = rearrange([qw(
    HGVS
    SLICE_ADAPTOR
    TRANSCRIPT_ADAPTOR
    REPLACE_REF
  )], @_);

  return $self->fetch_by_hgvs_notation(
    -hgvs               => $hgvs,
    -slice_adaptor      => $user_slice_adaptor,
    -transcript_adaptor => $user_transcript_adaptor,
    -replace_ref        => $replace_ref,
    -multiple_ok        => 1
  );
}


sub _hgvs_from_components {
  my ($self, $hgvs, $reference, $type, $description, $slice, $ref_allele, $alt_allele, $start, $end, $strand, $replace_ref) = @_;

  #######################  check alleles  #######################
  
  #Get the reference allele based on the coordinates - need to supply lowest coordinate first to slice->subseq()
  my $refseq_allele;
  if($start > $end){ $refseq_allele = $slice->subseq($end,   $start, $strand);}
  else{              $refseq_allele = $slice->subseq($start, $end,  $strand);}

  # take reference allele from genomic reference & coordinates if not supplied in HGVS string for a deletion
  $ref_allele = $refseq_allele  if( $type =~ m/g|c|n|m/ && $description =~ m/del/ );

  # If the reference allele was omitted, set it to undef
  $ref_allele = undef unless (defined($ref_allele) && length($ref_allele));      
  
  # take alternate allele from genomic reference & coordinates if not supplied in HGVS string for a duplication
  if($description =~ /dup/){ 
     ## special case: handle as insertion for ensembl object purposes 
     $start = $end ;
     if($strand  == 1){ $start++; }
     else{             $end--;   }

    $ref_allele = "-" ;
    $alt_allele = $refseq_allele;
  }
  elsif ($description =~ m/ins/i && $description =~ m/del/i) {
    ## check for malformed indels
    ($ref_allele, $alt_allele,  $start, $end) = @{trim_sequences($ref_allele, $alt_allele,  $start, $end, 1)};
  }
  elsif ($description =~ m/ins/i && $description !~ m/del/i) {
     # insertion: the start & end positions are inverted by convention
      if($end > $start){ ($start, $end  ) = ( $end , $start); }   
  }
 
  else {    
    if($DEBUG==1){print "Reference allele: $refseq_allele expected allele: $ref_allele\n";}
    
    if($replace_ref && defined($ref_allele)) {
      $ref_allele = $refseq_allele;
    }
    # If the reference from the sequence does not correspond to the reference given in the HGVS notation, throw an exception 
    elsif (defined($ref_allele) && $ref_allele ne $refseq_allele){        
      throw ("Reference allele extracted from $reference:$start-$end ($refseq_allele) does not match reference allele given by HGVS notation $hgvs ($ref_allele)");
    }
  }
  if (defined($ref_allele) && $ref_allele eq $alt_allele){         
    throw ("Reference allele extracted from $reference:$start-$end ($refseq_allele) matches alt allele given by HGVS notation $hgvs ($alt_allele)");
  }

  
  ####################### Create objects #######################

  #Create Allele objects
  my @allele_objs;
  if($self->db) {
    foreach my $allele ($ref_allele,$alt_allele) {
       push(@allele_objs,Bio::EnsEMBL::Variation::Allele->new('-adaptor' => $self, '-allele' => $allele));
    }
  }

  #Create a variation object. Use the HGVS string as its name
  my $variation = Bio::EnsEMBL::Variation::Variation->new(
    '-adaptor' => $self->db ? $self->db()->get_VariationAdaptor() : undef,
    '-name'    => $hgvs,
    '-source'  => 'Parsed from HGVS notation',
    '-alleles' => \@allele_objs
  );

  #Create a variation feature object
  my $variation_feature = Bio::EnsEMBL::Variation::VariationFeature->new(
    '-adaptor'       => $self,
    '-start'         => $start,
    '-end'           => $end,
    '-strand'        => $strand,
    '-slice'         => $slice,
    '-map_weight'    => 1,
    '-variation'     => $variation,
    '-allele_string' => "$ref_allele/$alt_allele"
  );
    
  if($DEBUG==1){print "Created object $hgvs allele_string: $ref_allele/$alt_allele, start:$start, end:$end\n";}

  return $variation_feature;
}

## when a user gives a HGVS notation on a gene
## we want to pick the "most likely" transcript
sub _pick_likely_transcript {
  my $self = shift;
  my $transcripts = shift;
  
  return $transcripts->[0] if scalar @$transcripts == 1;
  
  my @tr_info;

  foreach my $tr(@$transcripts) {

    # create a hash of info for this transcript that will be used to rank it
    my $info = {
      canonical => 1,
      ccds => 1,
      length => 0,
      biotype => 1,
      tsl => 100,
      appris => 100,
      tr => $tr
    };
     
    # 0 is "best"
    $info->{canonical} = $tr->is_canonical ? 0 : 1;
    $info->{biotype} = $tr->biotype eq 'protein_coding' ? 0 : 1;
    $info->{ccds} = (grep {$_->database eq 'CCDS'} @{$tr->get_all_DBEntries}) ? 0 : 1;
    
    # "invert" length so longer is best
    $info->{length} = 0 - $tr->length();
    
    # lower TSL is best
    if(my ($tsl) = @{$tr->get_all_Attributes('TSL')}) {
      if($tsl->value =~ m/tsl(\d+)/) {
        $info->{tsl} = $1 if $1;
      }
    }

    # same for APPRIS
    if(my ($appris) = @{$tr->get_all_Attributes('APPRIS')}) {
      if($appris->value =~ m/(.+?)(\d+)/) {
        if(my $num = $2) {
          $num += 10 unless $1 eq 'principal';
          $info->{appris} = $num if $num;
        }
      }
    }
    
    push @tr_info, $info;
  }
  
  my @order = qw(appris tsl canonical biotype ccds rank length);
  my $picked;
  
  # go through each category in order
  foreach my $cat(@order) {
    
    # sort on that category
    @tr_info = sort {$a->{$cat} <=> $b->{$cat}} @tr_info;
    
    # take the first (will have the lowest value of $cat)
    $picked = shift @tr_info;
    my @tmp = ($picked);
    
    # now add to @tmp those vfoas that have the same value of $cat as $picked
    push @tmp, shift @tr_info while @tr_info && $tr_info[0]->{$cat} eq $picked->{$cat};
    
    # if there was only one, return
    return $picked->{tr} if scalar @tmp == 1;
    
    # otherwise shrink the array to just those that had the lowest
    # this gives fewer to sort on the next round
    @tr_info = @tmp;
  }
  
  # probably shouldn't get here, but if we do, return the first
  return $tr_info[0]->{tr};
}

## Extract enough information to make a variation_feature from HGVS protein nomenclature
## Only attempts substitutions
##    - assumes protein change results from minimum number of nucleotide changes
##    - returns VF information only if one minimal solution found
sub _parse_hgvs_protein_position{

  my ($description, $reference, $transcript ) = @_;

  ## only supporting the parsing of hgvs substitutions [eg. Met213Ile]
  my ($from, $pos, $to) = $description =~ /^(\w+?)(\d+)(\w+?|\*)$/; 

  throw("Could not parse HGVS protein notation " . $reference . ":p.". $description ) unless $from and $pos and $to;

  # convert three letter AA to single
  $from = $Bio::SeqUtils::ONECODE{$from} || $from;
  $to   = $Bio::SeqUtils::ONECODE{$to} || $to;

  # get genomic position 
  my $tr_mapper = $transcript->get_TranscriptMapper(); 

  my @coords = $tr_mapper->pep2genomic($pos, $pos); 

  throw ("Unable to map the peptide coordinate $pos to genomic coordinates for protein $reference") if (scalar(@coords) != 1 || !$coords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate')); 

  my $strand = $coords[0]->strand(); 
  my $start  = $strand > 0 ? $coords[0]->start() : $coords[0]->end(); 
  my $end    = $strand > 0 ? $coords[0]->start() : $coords[0]->end(); 

  ## find reference sequence 
  my $slice = $transcript->slice();

  ## make a small slice for sequence look-up
  my $from_slice = Bio::EnsEMBL::Slice->new(-coord_system => $slice->coord_system(),
                                            -start => $coords[0]->start(),
                                            -end => $coords[0]->start()+2,
                                            -strand => $slice->strand(),
                                            -seq_region_name => $slice->seq_region_name,
                                            -seq_region_length => 3,
                                            -adaptor => $slice->adaptor);

  my $from_codon_ref = $from_slice->seq();

  throw ("Unable to find the genomic reference sequence for protein $reference") unless defined $from_codon_ref; 

  ## correct for strand
  reverse_comp(\$from_codon_ref) if $strand <0;

  # get correct codon table 
  my $attrib = $transcript->slice->get_all_Attributes('codon_table')->[0];

  # default to the vertebrate codon table which is denoted as 1 
  my $codon_table = Bio::Tools::CodonTable->new( -id => ($attrib ? $attrib->value : 1)); 

  # check genomic codon is compatible with input HGVS
  my $check_prot   = $codon_table->translate($from_codon_ref);
  my @from_codons;
  if ($check_prot eq $from){
    push @from_codons, $from_codon_ref ;
  }
  else{
    # rev-translate input ref sequence if the genome sequence does not match
    @from_codons   = $codon_table->revtranslate($from);
  }

  # rev-translate alt sequence
  my @to_codons   = $codon_table->revtranslate($to); 

  # now iterate over all possible mutation paths 
  my %paths; 
  foreach my $from_codon (@from_codons) {
    foreach my $to_codon (@to_codons) { 

      my $key = $from_codon .'_'. $to_codon;
      for my $i(0..2) { 
       
        my ($a, $b) = (substr($from_codon, $i, 1), substr($to_codon, $i, 1)); 
        next if uc($a) eq uc($b); 
        push @{$paths{$key}}, $i.'_'.uc($a).'/'.uc($b); 
      } 

      # non consecutive paths 
      if(scalar @{$paths{$key}} == 2 and $paths{$key}->[0] =~ /^0/ and $paths{$key}->[1] =~ /^2/) { 
         splice(@{$paths{$key}}, 1, 0, '1_'.substr($from_codon, 1, 1).'/'.substr($from_codon, 1, 1)); 
      } 

      $paths{$key} = join ",", @{$paths{$key}}; 
    }
  } 

  # get shortest dist and best paths with that dist 
  my $shortest_dist = length((sort {length($a) <=> length($b)} values %paths)[0]); 
  my %best_paths = map {$_ => 1} grep {length($_) eq $shortest_dist} values %paths;

  throw("Could not determine nucleotide change from peptide change $from \-\> $to") unless scalar keys %best_paths;

  my @results;

  foreach my $best_path(keys %best_paths) {

    my ($ref_allele, $alt_allele) = ('', '');
    my ($this_start, $this_end) = ($start, $end);
    my @path = split(/\,/, $best_path);

    # coords
	  if($strand > 0) {
  		$this_start += (split /\_/, $path[0])[0]; 
  		$this_end   += (split /\_/, $path[-1])[0];
	  }
	  else {
    	$this_start -= (split /\_/, $path[0])[0];
    	$this_end   -= (split /\_/, $path[-1])[0];
	  }

    # alleles 
    $ref_allele .= (split /\_|\//, $path[$_])[1] for 0..$#path;
    $alt_allele .= (split /\_|\//, $path[$_])[2] for 0..$#path;  

    push @results, [$ref_allele, $alt_allele, $this_start, $this_end, $strand];
  }

  return \@results;
  # 
  #use Data::Dumper; 
  #$Data::Dumper::Maxdepth = 3; 
  #warn Dumper \@from_codons; 
  #warn Dumper \@to_codons; 
  #warn Dumper \%paths; 
  #warn Dumper \%best_paths; 
  #exit(0); 
}


=head2 fetch_by_dbID

  Arg [1]    : string $id
               The unique database identifier for the VariationFeature to be obtained
  Example    : $vf = $adaptor->fetch_by_dbID(1234));
  Description: Returns the feature created from the database defined by the
               the id $id.  ID may be either a variation_feature_id or a
               "location_identifier" chr:start:allele_string:source_name
  Returntype : Bio::EnsEMBL::VariationFeature or undef
  Exceptions : thrown if $id arg is not provided
               does not exist
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $id = shift;

  throw("ERROR: No dbID given\n") unless $id;

  if(looks_like_number($id)) {
    return $self->SUPER::fetch_by_dbID($id);
  }
  else {
    my $vfs = $self->fetch_all_by_location_identifier($id);
    return @$vfs ? $vfs->[0] : undef;
  }
}


=head2 fetch_all_by_location_identifier

  Arg [1]    : string $location_identifier
  Example    : $vf = $adaptor->fetch_by_dbID('1:230710048:A_G');
  Description: Fetches VariationFeatures by location identifier.
               Primarily used to fetch variants from VCFCollections
               as optional 4th component is source_name or
               VCFCollection name.

               location_identifier format chr:start:allele_string:source_name

               chr = chromosome name
               start = start coordinate of variant
               allele_string = "_"-separated allele string (optional)
               source_name = source name or name of VCFCollection (optional)
               
  Returntype : listref of Bio::EnsEMBL::VariationFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_location_identifier {
  my $self = shift;
  my $vl = shift;

  my $sa = $self->db->dnadb->get_SliceAdaptor();
  return [] unless $sa;

  my $vfs;
      
  my ($chr, $start, $alleles, $vcf_id) = split(':', $vl);
  my $slice = $sa->fetch_by_region(undef, $chr, $start, $start);

  if($slice) {
    $vfs = $self->fetch_all_by_Slice($slice);

    # filter based on $alleles and $vcf_id
    if($alleles) {
      $alleles =~ s/\_/\//g;
      @$vfs = grep {$_->allele_string eq $alleles} @$vfs if $alleles;
    }

    @$vfs = grep {$_->source_name eq $vcf_id} @$vfs if $vcf_id;
  }

  return $vfs;
}

1;
