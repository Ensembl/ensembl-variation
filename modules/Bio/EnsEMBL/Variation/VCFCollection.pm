=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

#
# Ensembl module for Bio::EnsEMBL::Variation::VCFCollection
#
#

=head1 NAME

Bio::EnsEMBL::DBSQL::VCFCollection

=head1 SYNOPSIS
  my $reg = 'Bio::EnsEMBL::Registry';

  # get VCF Collection Adaptor
  my $vca = $reg->get_adaptor('human', 'variation', 'vcfcollection');

  # iterate over collections
  foreach my $c(@{$vca->fetch_all})

    # get samples
    my $samples = $c->get_all_Samples;

    # get genotypes for a VariationFeature
    my $gts = $c->get_all_SampleGenotypeFeatures_by_VariationFeature($vf);
  }

=head1 DESCRIPTION

This module represents a collection of VCF files, e.g. one for each chromosome
from the 1000 Genomes set. Each file is represented by a
Bio::EnsEMBL::IO::Parser::VCF4Tabix object

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::VCFCollection;

use Cwd;
use Scalar::Util qw(weaken);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Scalar qw(check_ref assert_ref);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::VEP qw(parse_line);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(MAX_DISTANCE_FROM_TRANSCRIPT);

use Bio::EnsEMBL::IO::Parser::VCF4Tabix;
use Bio::EnsEMBL::Variation::SampleGenotypeFeature;
use Bio::EnsEMBL::Variation::Sample;
use Bio::EnsEMBL::Variation::Individual;
use Bio::EnsEMBL::Variation::Population;
use Bio::EnsEMBL::Variation::VCFVariationFeature;
use Bio::EnsEMBL::Variation::IntergenicVariation;

use base qw(Bio::EnsEMBL::Variation::BaseAnnotation);

my $MAX_OPEN_FILES = 2;

=head2 new

  Arg [-ID]:                     string - identifier for this collection
  Arg [-DESCRIPTION]:            string - description for this collection
  Arg [-TYPE]:                   string - "local" or "remote"
  Arg [-USE_AS_SOURCE]:          boolean
  Arg [-FILENAME_TEMPLATE]:      string
  Arg [-CHROMOSOMES]:            arrayref of strings
  Arg [-SAMPLE_PREFIX]:          string
  Arg [-POPULATION_PREFIX]:      string
  Arg [-SAMPLE_POPULATIONS]:     hashref - { 'sample1': ['pop1','pop2'] }
  Arg [-STRICT_NAME_MATCH]:      boolean
  Arg [-REF_FREQ_INDEX]:         int - index position of ref frequency in INFO field, if given
  Arg [-USE_SEQ_REGION_SYNONYMS]:boolean
  Arg [-ADAPTOR]:                Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor

  Example    : my $collection = Bio::EnsEMBL::Variation::VCFCollection->new(
                -id                 => 'test',
                -type               => 'local',
                -filename_template  => '/path/to/vcfs/test_###CHR###.vcf.gz',
                -chromosomes        => [1, 2, 3],
                -sample_populations => {
                    'sample1': ['pop1','pop2'],
                    'sample2': ['pop3']
                 }
               );

  Description: Constructor.  Instantiates a new VCFCollection object.
  Returntype : Bio::EnsEMBL::Variation::DBSQL::VCFCollection
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my (
    $sample_prefix,
    $individual_prefix,
    $pop_prefix,
    $sample_pops,
    $populations,
    $ref_freq_index,
  ) = rearrange(
    [qw(
      SAMPLE_PREFIX
      INDIVIDUAL_PREFIX
      POPULATION_PREFIX
      SAMPLE_POPULATIONS
      POPULATIONS
      REF_FREQ_INDEX
    )],
    @_
  ); 
 
 my $self = $class->SUPER::new(@_);

  $self->sample_prefix($sample_prefix);
  $self->population_prefix($pop_prefix);
  $self->{individual_prefix} = $individual_prefix;
  $self->{populations} = $populations;
  $self->{_raw_populations} = $sample_pops;
  $self->ref_freq_index($ref_freq_index);
  $self->{_use_db} = 1; 
  
  return $self;
}

=head2 ref_freq_index

  Arg [1]    : int $ref_freq_index (optional)
               The new value to set the attribute to
  Example    : my $ref_freq_index = $collection->ref_freq_index()
  Description: Getter/Setter for the ref_freq_index parameter. Some VCF
               files contain the reference (REF) allele frequency in those
               listed in AF-type INFO fields. Some (ESP) give the reference
               frequency last, others (dbSNP) give it as the first. Set this
               to a defined value to extract the REF frequency from the
               given indexed position in the list given in the INFO field
               value. -1 may be used if the REF is the last.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub ref_freq_index {
  my $self = shift;
  $self->{ref_freq_index} = shift if @_;
  return $self->{ref_freq_index};
}

=head2 list_chromosomes

  Example    : my $chrs = $collection->list_chromosomes()
  Description: Get list of chromosome names covered by this collection.
  Returntype : arrayref of strings
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub list_chromosomes {
  my $self = shift;

  if(!defined($self->{chromosomes})) {

    # if the collection only contains one file we can use tabix to get the chr list
    # test this by passing a dummy value to get the filename
    my $dummy_chr = $$.'_dummy_VCFCollection_chr';

    my $fn = $self->_get_vcf_filename_by_chr($dummy_chr);

    if($fn !~ /$dummy_chr/) {
      $self->{chromosomes} = $self->_vcf_parser_obj($fn)->{tabix_file}->seqnames;
    }
    else {
      $self->{chromosomes} = [];
    }
  }
  return $self->{chromosomes};
}

=head2 use_db

  Arg [1]    : int $use_db (optional)
               The new value to set the use_db attribute to.
  Example    : my $use_db = $collection->use_db()
  Description: Getter/Setter for the use_db attribute of this collection.
               If set to 1 (default), the API will attempt to retrieve
               Sample and Population objects from the database, using
               sample_prefix and population_prefix as appropriate.
               If set to 0, the API will create "fake" Sample and
               Population objects. Fake objects will also be created if
               the DB fetch fails for a sample.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub use_db {
  my $self = shift;
  $self->{_use_db} = shift if @_;
  return $self->{_use_db};
}

=head2 sample_prefix

  Arg [1]    : string $sample_prefix (optional)
               The new value to set the sample_prefix attribute to
  Example    : my $sample_prefix = $collection->sample_prefix()
  Description: Getter/Setter for the sample_prefix of this collection.
               This property can be useful to align sample names from
               the VCF header with entries in the sample database table.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub sample_prefix {
  my $self = shift;
  $self->{sample_prefix} = shift if @_;
  return $self->{sample_prefix} || $self->{individual_prefix} || '';
}


=head2 population_prefix

  Arg [1]    : string $population_prefix (optional)
               The new value to set the population_prefix attribute to
  Example    : my $population_prefix = $collection->population_prefix()
  Description: Getter/Setter for the population_prefix of this collection.
               This property can be useful to align population names from
               the config file with entries in the population database table.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub population_prefix {
  my $self = shift;
  $self->{population_prefix} = shift if @_;
  return $self->{population_prefix} || '';
}


=head2 get_all_Samples

  Example    : my $samples = $collection->get_all_Samples()
  Description: Get all sample objects that will have genotypes in
               this collection.
  Returntype : arrayref of Bio::EnsEMBL::Variation::Sample
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Samples {
  my $self = shift;
  if(!defined($self->{samples})) {

    my @samples;
    
    # we should only need to get samples from one chromosome's VCF
    my $chr = $self->list_chromosomes ? $self->list_chromosomes->[0] : '';
    my $vcf = $self->_get_vcf_by_chr($chr);
    
    throw("ERROR: VCF file ".$self->_get_vcf_filename_by_chr($chr)." not found\n") unless $vcf;
    
    my $prefix = $self->sample_prefix;
    my $sample_adpt = $self->use_db ? $self->adaptor->db->get_SampleAdaptor() : undef;
    my $individual_adpt = $self->use_db ? $self->adaptor->db->get_IndividualAdaptor() : undef;
    
    my $vcf_sample_names = $vcf->get_samples();
    
    # do a fetch_all_by_name_list
    my %sample_objs;
    %sample_objs = map {$_->name() => $_} @{$sample_adpt->fetch_all_by_name_list([map {$prefix.$_} @$vcf_sample_names])} if $self->use_db;
   
    # Get the list of sample synonyms
    my %synonyms;
    for my $sample_name (keys %sample_objs) {
      my $sample = $sample_objs{$sample_name};
      for my $syn (@{$sample->get_all_synonyms}) {
        $synonyms{$syn} = $sample->name;
      }
    }

    # some may not be in DB
    foreach my $vcf_sample_name (@$vcf_sample_names) {
      # Use the main sample name to retrieve its metadata
      my $sample_name = $synonyms{ $vcf_sample_name } || $vcf_sample_name;
      
      # either use the DB one or create one
      my $sample = $sample_objs{ $prefix.$sample_name } ||
        Bio::EnsEMBL::Variation::Sample->new_fast({
          name            => $prefix.$sample_name,
          adaptor         => $sample_adpt,
          display         => 'UNDISPLAYABLE',
          dbID            => --($self->{_sample_id}),
          individual      => Bio::EnsEMBL::Variation::Individual->new_fast({
            name     => $prefix.$sample_name,
            adaptor  => $individual_adpt,
            type_individual => 'outbred',
            dbID     => --($self->{_sample_id}),
          }),
        });
      # store the raw name to easily match to data returned from other methods
      $sample->{_raw_name} = $vcf_sample_name;
      push @samples, $sample;
    }

    $self->{samples} = \@samples;
  }
  
  return $self->{samples};
}


=head2 get_all_Populations

  Example    : my $pops = $collection->get_all_Populations()
  Description: Get all population objects for the samples in this
               collection.
  Returntype : arrayref of Bio::EnsEMBL::Variation::Population
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Populations {
  my $self = shift;

  if(!exists($self->{populations}) || !defined $self->{populations}->[0] ) {
    my $hash = $self->_get_Population_Sample_hash;

    if( $self->use_db) {
      my $pa = $self->adaptor->db->get_PopulationAdaptor;
      $self->{populations} = $pa->fetch_all_by_dbID_list([keys %$hash]);
    }
  }
  
  return $self->{populations};
}


=head2 has_Population
  Arg[1]     : string $pop OR Bio::EnsEMBL::Variation::Population $pop
  Example    : my $has_pop = $collection->has_Population($pop)
  Description: Returns true if this collection contains this population
  Returntype : boolean
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub has_Population {
  my $self = shift;
  my $pop = shift;
  
  my $name = ref($pop) eq '' ? $pop : $pop->name;
  
  return grep {$name eq $_} @{$self->_get_all_population_names};
}


=head2 get_all_VariationFeatures_by_Slice

  Arg[1]     : Bio::EnsEMBL::Slice $slice
  Arg[2]     : (optional) bool $dont_fetch_overlaps
  Example    : my $vfs = $collection->get_all_VariationFeatures_by_Slice($slice)
  Description: Get all VariationFeatures (actually VCFVariationFeatures) for this
               slice. By default feature overlap objects are also generated for
               any overlapping Transcripts, RegulatoryFeatures and MotifFeatures.
               Set $dont_fetch_overlaps to a true value to disable this;
               VariationFeatures in this scenario will have their consequence type
               set to "intergenic".
  Returntype : arrayref of Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : none
  Caller     : Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor
  Status     : Stable

=cut

sub get_all_VariationFeatures_by_Slice {
  my $self = shift;
  my $slice = shift;
  my $dont_fetch_vf_overlaps = shift;
  # Variants from a VCF file are always on the forward strand. If a slice on the reverse strand is passed to the sub we
  # create the same slice but on the forward strand using fetch_by_region and attach it to the VCFVariationFeature object.
  # Attaching the slice on the reverse strand to the VCFVariationFeature object would cause problems when
  # calculating the consequence for an overlapping transcript which is also on the reverse strand. The variant alleles would not be reverse complemented
  # correctly because the code assumes based on the negative slice that the alleles are already reverse complemented. This can lead to wrongly calculated
  # consequences.
  if ($slice->strand == -1) {
    my $sa = $self->adaptor->db->dnadb->get_SliceAdaptor();
    $slice = $sa->fetch_by_region($slice->coord_system_name, $slice->seq_region_name, $slice->start, $slice->end);
  }
  return [] unless $self->_seek_by_Slice($slice);
  
  my $vcf = $self->_current();
  my $vfa = $self->use_db ? $self->adaptor->db->get_VariationFeatureAdaptor : Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor->new_fake($self->species);
  
  my @vfs;
  
  while($vcf->{record} && $vcf->get_start <= $slice->end) {
    
    my $copy = $vcf->get_frozen_copy();
    
    foreach my $parsed_vf(@{parse_line({format => 'vcf'}, join("\t", @{$vcf->{record}}))}) {
      next unless $parsed_vf->isa('Bio::EnsEMBL::Variation::VariationFeature');
      delete $parsed_vf->{_line};
      
      my $vcf_vf = Bio::EnsEMBL::Variation::VCFVariationFeature->new_from_VariationFeature(
        -variation_feature => $parsed_vf,
        -collection => $self,
        -vcf_record => $copy,
        -slice => $slice,
        -adaptor => $vfa
      );
      
      push @vfs, $vcf_vf;
    }
    
    $vcf->next();
  }

  if($dont_fetch_vf_overlaps || !$self->use_db) {
    foreach my $vf(@vfs) {
      $vf->{intergenic_variation} = Bio::EnsEMBL::Variation::IntergenicVariation->new(
        -variation_feature  => $vf,
        -no_ref_check       => 1,
      );
      weaken($vf->{intergenic_variation}->{base_variation_feature});

      $vf->_finish_annotation();
    }
  }

  else {

    ## we need to up-front fetch overlapping transcripts, regfeats and motiffeatures
    ## this prevents the API loading them per-variant later at great cost in speed
    ## not an easy way to do this generically in a loop, so done type-by-type
    my $db = $self->adaptor->db;

    # transcripts
    my @transcripts =
      map {$_->transfer($slice)}
      @{$slice->expand(MAX_DISTANCE_FROM_TRANSCRIPT, MAX_DISTANCE_FROM_TRANSCRIPT)->get_all_Transcripts(1)};

    $db->get_TranscriptVariationAdaptor->fetch_all_by_VariationFeatures(
      \@vfs,
      \@transcripts
    ) if @transcripts;

    # funcgen types
    foreach my $type(qw(RegulatoryFeature MotifFeature)) {
      if(
        my $fg_adaptor = Bio::EnsEMBL::DBSQL::MergedAdaptor->new(
          -species  => $db->species, 
          -type     => $type,
        )
      ) {
        my $get_adaptor_method = 'get_'.$type.'VariationAdaptor';
        my @features = map {$_->transfer($slice)} @{$fg_adaptor->fetch_all_by_Slice($slice)};

        $db->$get_adaptor_method->fetch_all_by_VariationFeatures(
          \@vfs,
          \@features,
        ) if @features;
      }
    }

    # this "fills in" the hashes that store e.g. TranscriptVariations
    # so that the API doesn't go and try to fill them again later
    $_->_finish_annotation() for @vfs;
  }
  
  return \@vfs;
}


=head2 get_all_SampleGenotypeFeatures_by_VariationFeature

  Arg[1]     : Bio::EnsEMBL::Variation::VariationFeature $vf
  Arg[2]     : (optional) Bio::EnsEMBL::Variation::Population OR
               Bio::EnsEMBL::Variation::Sample
  Example    : my $gts = $collection->get_all_SampleGenotypeFeatures_by_VariationFeature($vf)
  Description: Get all SampleGenotypeFeatures for a given
               VariationFeature object
  Returntype : arrayref of Bio::EnsEMBL::Variation::SampleGenotypeFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_SampleGenotypeFeatures_by_VariationFeature {
  my $self = shift;
  my $vf = shift;
  my $sample = shift;
  my $vcf = shift;
  
  assert_ref($vf, 'Bio::EnsEMBL::Variation::VariationFeature');
  
  if(!$vcf) {
    # seek to record for VariationFeature
    return [] unless $self->_seek_by_VariationFeature($vf); 
    $vcf = $self->_current();
  }
  else {
    $self->_current($vcf);
  }
  
  my $samples = $self->_limit_Samples($self->get_all_Samples, $sample);
  my @sample_names = map {$_->{_raw_name}} @$samples;
  
  return [] unless scalar @$samples;

  return $self->_create_SampleGenotypeFeatures(
    $samples,
    $vcf->get_samples_genotypes(\@sample_names),
    $vf
  );
}


=head2 get_all_Alleles_by_VariationFeature

  Arg[1]     : Bio::EnsEMBL::Variation::VariationFeature $vf
  Arg[2]     : (optional) Bio::EnsEMBL::Variation::Population $pop
  Arg[3]     : (optional) Bio::EnsEMBL::IO::Parser::VCF4Tabix $vcf
  Example    : my $alleles = $collection->get_all_Alleles_by_VariationFeature($vf)
  Description: Get all Alleles for a given VariationFeature object
  Returntype : arrayref of Bio::EnsEMBL::Variation::Allele
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Alleles_by_VariationFeature {
  my $self = shift;
  my $vf = shift;
  my $given_pop = shift;
  my $vcf = shift;

  assert_ref($vf, 'Bio::EnsEMBL::Variation::VariationFeature');
  assert_ref($given_pop,'Bio::EnsEMBL::Variation::Population') if (defined($given_pop));

  # if given $vcf, we don't want to call next() on it
  # as it is a frozen copy created in VCFVF object creation
  my $next_ok = 1;
  
  if(!$vcf) {
    # seek to record for VariationFeature
    return [] unless $self->_seek_by_VariationFeature($vf);    
    $vcf = $self->_current();
  }
  else {
    $self->_current($vcf);
    $next_ok = 0;
  }

  # get data from VF
  my $vf_ref = $vf->ref_allele_string;
  my $vf_alts = $vf->alt_alleles;

  # check initial VCF match
  my $matched = $self->_get_matched_alleles_VF_VCF($vf, $vcf);
  return [] unless @$matched;

  # get generic pop data and collection-specific param
  my @pops = @{$self->get_all_Populations};
  @pops = grep {$_->name eq $given_pop->name || ($_->{_raw_name} || '') eq $given_pop->name} @pops if $given_pop;
  my $ref_freq_index = $self->ref_freq_index();

  # log any alleles present in the VCF not present in the VF alleles
  my $missing_alleles = {};
  my $freqs = {};
  my $counts = {};
  my $total_acs = {};
  my $total_afs = {};
  my $ans = {};

  # loop over VCF as there may be more than one VCF line containing freq data for this VF
  # e.g. VF is G/A/C
  # VCF line 1 is G/A, line 2 is G/C
  while(1) {

    # get data from VCF entry
    my $info = $vcf->get_info;
    my $vcf_alts = $vcf->get_alternatives;
    my %allele_map = map {$_->{b_index} => $_->{a_allele}} @$matched;
     POP : {
       foreach my $pop(@pops) {

        # this allows for an empty population name to be looked up as e.g. AF
        my $raw_name = exists($pop->{_raw_name}) ? $pop->{_raw_name} : $pop->name;
        my $suffix = $raw_name ? '_'.$raw_name : '';
        my $pop_id = $pop->dbID;
        # Due to remapping errors of gnomAD data there are cases
        # where there are more than one VCF line with the same alt.
        # This check handles multiple lines with the same alt
        # and will only consider the first VCF line and ignore the rest.
        foreach my $allele (values %allele_map) {
          if (exists $freqs->{$pop_id}->{$allele}) {
            last POP;
          }
        }


        no warnings 'uninitialized';
        my ($ac, $an, $af) = (
          $info->{'AC'.$suffix} // $info->{$pop->{_ac}},
          $info->{'AN'.$suffix} // $info->{$pop->{_an}},
          $info->{'AF'.$suffix} // $info->{$pop->{_af}},
        );
        if (defined($ac) && $ac ne '.' && !defined($an)) {
          my @split = split(',', $ac);
          $an += $_ for @split;
        }
        # log AN for later use
        $ans->{$pop_id} = $an if $an;

        # Have AC and AN (AN must be defined and non-zero)
        # Used first - instead of using AF first
        if (defined($ac) && $ac ne '.' && $an && $an ne '.') {
          $total_acs->{$pop_id} ||= 0;
          my @split = split(',', $ac);
          # is ref AC included? This may have been set as ref_freq_index()
          # or we can auto-detect by comparing size of @split to @$vcf_alts
          if (defined($ref_freq_index) || scalar @split > scalar @$vcf_alts) {
            my $ref_allele_count = splice(@split, defined($ref_freq_index) ? $ref_freq_index : -1, 1);
            $freqs->{$pop_id}->{$vf_ref} = sprintf('%.4g', $ref_allele_count  / $an);
            $counts->{$pop_id}->{$vf_ref} = $ref_allele_count;
          }

          for my $i(0..$#split) {
            my $c = $split[$i];
            $total_acs->{$pop_id} += $c;
            if(my $allele = $allele_map{$i}) {
              $counts->{$pop_id}->{$allele} = $c;
              $freqs->{$pop_id}->{$allele} = sprintf('%.4g', $c / $an);
            }
            else {
              $missing_alleles->{$pop_id}->{$vcf_alts->[$i]} = 1;
            }
          }
        }
        # check for AF
        elsif (defined($af) && $af ne '.') {
          $total_afs->{$pop_id} ||= 0;
          my @split = split(',', $af);

          # is ref AC included? This may have been set as ref_freq_index()
          # or we can auto-detect by comparing size of @split to @$vcf_alts
          if (defined($ref_freq_index) || scalar @split > scalar @$vcf_alts) {
            $freqs->{$pop_id}->{$vf_ref} = splice(@split, defined($ref_freq_index) ? $ref_freq_index : -1, 1);
          }

          for my $i(0..$#split) {
            my $f = $split[$i];
            next if $f eq '.';
            $total_afs->{$pop_id} += $f;

            if(my $allele = $allele_map{$i}) {
              $freqs->{$pop_id}->{$allele} = $f;
              $counts->{$pop_id}->{$allele} = sprintf('%.0f', $f * $an) if $an;
            }
            else {
              $missing_alleles->{$pop_id}->{$vcf_alts->[$i]} = 1;
            }
          }
        }
      } # end foreach population
    }

    last unless $next_ok;

    # get next VCF entry
    $vcf->next();

    # check a record exists
    last unless $vcf->{record};

    # check matches current VF
    $matched = $self->_get_matched_alleles_VF_VCF($vf, $vcf);
    last unless @$matched;
  }

  # create the actual allele objects
  my @alleles;

  # we need these to create allele objects
  my $variation = $vf->variation;
  my $allele_adaptor = $self->use_db ? $self->adaptor->db->get_AlleleAdaptor : undef;
  foreach my $pop(@pops) {
    my $pop_id = $pop->dbID;

    # interpolate ref freq
    unless(exists($freqs->{$pop_id}->{$vf_ref})) {

      no warnings 'uninitialized';
      my ($ac, $an, $af) = (
        $total_acs->{$pop_id},
        $ans->{$pop_id},
        $total_afs->{$pop_id}
      );
      # have AC and AN, can work out freq and count
      if(defined($ac) && $an) {
        $counts->{$pop_id}->{$vf_ref} = $an - $ac;
        $freqs->{$pop_id}->{$vf_ref} = sprintf('%.4g', 1 - ($ac / $an));
      }

      # only have AF
      elsif(defined($af)) {
        $freqs->{$pop_id}->{$vf_ref} = 1 - $af;

        # if have AN can do counts too
        $counts->{$pop_id}->{$vf_ref} = sprintf('%.0f', $an * (1 - $af)) if $an;
      }
    }

    foreach my $a(keys %{$freqs->{$pop_id}}) {
      push @alleles, Bio::EnsEMBL::Variation::Allele->new_fast({
        allele     => $a,
        count      => $counts && $counts->{$pop_id} ? ($counts->{$pop_id}->{$a} || 0) : undef,
        frequency  => $freqs->{$pop_id}->{$a},
        population => $pop,
        variation  => $variation,
        adaptor    => $allele_adaptor,
        subsnp     => undef,
        _missing_alleles => $missing_alleles->{$pop_id} || {},
      })
    }
  }

  return \@alleles;
}

=head2 get_all_PopulationGenotypes_by_VariationFeature

  Arg[1]     : Bio::EnsEMBL::Variation::VariationFeature $vf
  Arg[2]     : (optional) Bio::EnsEMBL::Variation::Population $pop
  Arg[3]     : (optional) Bio::EnsEMBL::IO::Parser::VCF4Tabix $vcf
  Example    : my $population_genotypes = $collection->get_all_PopulationGenotypes_by_VariationFeature($vf)
  Description: Get all PopulationGenotypes for a given VariationFeature object
               PopulationGenotype counts and frequencies are only retrieved for ESP from the info field. 
  Returntype : arrayref of Bio::EnsEMBL::Variation::PopulationGenotype
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_PopulationGenotypes_by_VariationFeature {
  my $self = shift;
  my $vf = shift;
  my $given_pop = shift;
  my $vcf = shift;

  assert_ref($vf, 'Bio::EnsEMBL::Variation::VariationFeature');
  assert_ref($given_pop,'Bio::EnsEMBL::Variation::Population') if (defined($given_pop));

  # if given $vcf, we don't want to call next() on it
  # as it is a frozen copy created in VCFVF object creation
  my $next_ok = 1;

  if(!$vcf) {
    # seek to record for VariationFeature
    return [] unless $self->_seek_by_VariationFeature($vf);    
    $vcf = $self->_current();
  }
  else {
    $self->_current($vcf);
    $next_ok = 0;
  }
  # get data from VF
  my $vf_ref = $vf->ref_allele_string;
  my $vf_alts = $vf->alt_alleles;
  # check initial VCF match
  my $matched = $self->_get_matched_alleles_VF_VCF($vf, $vcf);
  return [] unless @$matched;

  my @pops = @{$self->get_all_Populations};
  if ($given_pop) {
    @pops = grep {$_->name eq $given_pop->name || ($_->{_raw_name} || '') eq $given_pop->name} @pops;
  }

  my $variation = $vf->variation;
  my $population_genotype_adaptor = $self->use_db ? $self->adaptor->db->get_PopulationGenotypeAdaptor : undef;

  my @population_genotypes = ();

  while(1) {
    my $info = $vcf->get_info;
    my $vcf_alts = $vcf->get_alternatives;
    my @alt_alleles = @$vcf_alts;
    my $ref = $vcf->get_reference;
    my %allele_map = map {$_->{b_index} => $_->{a_allele}} @$matched;
    my @sorted_gts = ();
    my $info_gts = $info->{'GTS'};
    if (defined $info_gts) {
      my @gts = split(',', $info_gts);
      foreach my $gt (@gts) {
        my @tmp_gt = ();
        if ($gt =~ /^[ACGT]+$/i) {
          my @gt_alleles = split('', $gt);
          foreach my $gt_allele (@gt_alleles) {
            push @tmp_gt, $gt_allele;
          }
        } elsif ($gt =~ /[A\d|R]+/) { 
          # Observed Genotypes. For INDELs, A1, A2, or An refers
          # to the N-th alternate allele while R refers to the reference allele.
          my @gt_variables  = split('', $gt);
          while (@gt_variables) {
            my $gt_variable = shift @gt_variables; 
            if ($gt_variable eq 'R') {
              push @tmp_gt, $vf_ref;
            } elsif ($gt_variable eq 'A') {
              my $idx = shift @gt_variables;
              if (my $allele = $allele_map{$idx - 1}) {
                push @tmp_gt, $allele;
              } else {
                my $allele = $vcf_alts->[$idx - 1];
                push @tmp_gt, $allele;
              }
            }  
            else {
              warn "Unexpected genotype variable. Expected values are R, A1, A2, An"
            }
          }
        } else {
          warn "Unexpected genotype variable. Expected values are R, A1, A2, An, or SNV represented by ACGT\n";
        } 
        my $sorted_gt = join('|', sort @tmp_gt);
        push @sorted_gts, $sorted_gt;
      }
      POP : {
        foreach my $pop (@pops) {
          my $raw_name = exists($pop->{_raw_name}) ? $pop->{_raw_name} : $pop->name;
          my $pop_id = $pop->dbID;
          my @gtc = split(',', $info->{$raw_name . '_GTC'});
          my $total_gtc = 0;
          $total_gtc += $_ for @gtc;
          my $idx = 0;
          foreach my $pop_gt (@sorted_gts) {
            if ($gtc[$idx]) {
              my @gt_alleles = sort {$a cmp $b} (split /\|/, $pop_gt);
              my $pop_gt_obj = Bio::EnsEMBL::Variation::PopulationGenotype->new_fast({
                  genotype   => \@gt_alleles,
                  population => $pop,
                  variation  => $variation,
                  frequency  => sprintf('%.4g', $gtc[$idx] / $total_gtc),
                  count      => $gtc[$idx],
                });
              push @population_genotypes, $pop_gt_obj;
            }
            $idx++;
          }
        }
      }
    }
    last unless $next_ok;

    # get next VCF entry
    $vcf->next();

    # check a record exists
    last unless $vcf->{record};

    # check matches current VF
    $matched = $self->_get_matched_alleles_VF_VCF($vf, $vcf);
    last unless @$matched;
  }
  return \@population_genotypes;
}


=head2 get_location_to_name_map

  Arg[1]     : Bio::EnsEMBL::Slice $slice
  Example    : my $location_to_name = $collection->get_location_to_name($slice);
  Description: Match variation name from a variation database to entry
               in the VCF file by location and alleles.
  Returntype : hashref
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_location_to_name_map {
  my $self = shift;
  my $slice = shift;

  return {} unless $self->_seek_by_Slice($slice);
  
  my $vcf = $self->_current();
  
  # get VariationFeatures if using DB
  my %vfs_by_pos;
  my $use_db = $self->use_db;
  my %location_to_name; 
  if($use_db || $self->use_as_source) {
    my $vfa = $self->adaptor->db->get_VariationFeatureAdaptor();
    foreach my $vf(@{$vfa->fetch_all_by_Slice($slice, 1)}) {
      push @{$vfs_by_pos{$vf->seq_region_start}}, $vf;
    }
    # reset seek
    $self->_seek_by_Slice($slice);
  }
  
  while($vcf->{record} && $vcf->get_start <= $slice->end) {
    my $start = $vcf->get_start;
    my $vf;
    # try to match this VCF record to VariationFeature at this position
    if(scalar keys %vfs_by_pos) {
      foreach my $tmp_vf(@{$vfs_by_pos{$start} || []}) {
        if (@{$self->_get_matched_alleles_VF_VCF($tmp_vf, $vcf)}) {
          $location_to_name{$vcf->get_raw_start} = $tmp_vf->variation_name;
          last;
        }
      }
    }
    $vcf->next();
  }

  return \%location_to_name;
}

=head2 get_all_SampleGenotypeFeatures_by_Slice

  Arg[1]     : Bio::EnsEMBL::Slice $slice
  Arg[2]     : (optional) Bio::EnsEMBL::Variation::Population OR
               Bio::EnsEMBL::Variation::Sample
  Arg[3]     : (optional) $non_ref_only
  Example    : my $gts = $collection->get_all_SampleGenotypeFeatures_by_Slice($slice)
  Description: Get all SampleGenotypeFeatures for a given
               genomic region represented by a slice. Set $non_ref_only to a true
               value to skip homozygous reference genotypes.
  Returntype : arrayref of Bio::EnsEMBL::Variation::SampleGenotypeFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_SampleGenotypeFeatures_by_Slice {
  my $self = shift;
  my $slice = shift;
  my $sample = shift;
  my $non_ref_only = shift;
  my $sample_cmp = shift;
  my $cmp = shift;

  return [] unless $self->_seek_by_Slice($slice);
  
  my $vcf = $self->_current();
  
  # get VariationFeatures if using DB
  my %vfs_by_pos;
  my $use_db = $self->use_db;
  
  if($use_db || $self->use_as_source) {
    my $vfa = $self->adaptor->db->get_VariationFeatureAdaptor();

    foreach my $vf(@{$vfa->fetch_all_by_Slice($slice, 1)}) {
      push @{$vfs_by_pos{$vf->seq_region_start}}, $vf;
    }

    # reset seek
    $self->_seek_by_Slice($slice);
  }
  
  my @genotypes;
  my $samples = $self->_limit_Samples($self->get_all_Samples, $sample);
  my @sample_names = map {$_->{_raw_name}} @$samples;
  
  return [] unless scalar @$samples;

  my $strict = $self->strict_name_match;
  
  while($vcf->{record} && $vcf->get_start <= $slice->end) {
    my $start = $vcf->get_start;
    
    my $vf;
    
    # try to match this VCF record to VariationFeature at this position
    if(scalar keys %vfs_by_pos) {

      if($strict) {
        foreach my $tmp_vf(@{$vfs_by_pos{$start} || []}) {
          $vf = $tmp_vf if(grep {$tmp_vf->variation_name eq $_ || $tmp_vf->variation_name eq 'ss'.$_} @{$vcf->get_IDs});
          last if $vf;
        }
      }
      else {
        foreach my $tmp_vf(@{$vfs_by_pos{$start} || []}) {
          $vf = $tmp_vf if @{$self->_get_matched_alleles_VF_VCF($tmp_vf, $vcf)};
          last if $vf;
        }
      }
    }
    
    # otherwise create a VariationFeature object
    else {
      $vf = parse_line({format => 'vcf'}, join("\t", @{$vcf->{record}}))->[0];
    }
    
    if($vf) {
      push @genotypes, @{$self->_create_SampleGenotypeFeatures($samples, $vcf->get_samples_genotypes(\@sample_names, $non_ref_only), $vf, $sample_cmp, $cmp)};
    }
    
    $vcf->next();
  }

  return \@genotypes;
}

sub _limit_Samples {
  my $self = shift;
  my $samples = shift;
  my $sample = shift;
  my $limited = $samples;
  
  # limit by sample or population?
  if(defined($sample)) {
    if($sample->isa('Bio::EnsEMBL::Variation::Sample')) {
      $limited = [grep {$_->name eq $sample->name} @$samples];
    }
    elsif($sample->isa('Bio::EnsEMBL::Variation::Population')) {
      my %limit = map {$_->name => 1} @{$sample->get_all_Samples};
      $limited = [grep {defined($limit{$_->name})} @$samples];
    }
    else {
      throw("Argument $sample is not a Bio::EnsEMBL::Variation::Sample or ::Popluation");
    }
  }
  
  return $limited;
}

sub _create_SampleGenotypeFeatures {
  my $self = shift;
  my $samples = shift;
  my $raw_gts = shift;
  my $vf = shift;
  my $sample_cmp = shift;
  my $cmp = shift;
  
  # $vf could be a StructuralVariationFeature if it is generated by parsing the VCF
  return [] unless check_ref($vf, 'Bio::EnsEMBL::Variation::VariationFeature');
  
  my @genotypes;
  
  $self->{_gta} ||= $self->use_db ? $self->adaptor->db->get_SampleGenotypeFeatureAdaptor : undef;

  my $vcf = $self->_current();

  my $matched = $self->_get_matched_alleles_VF_VCF($vf, $vcf);
  my %allele_map = map {$_->{b_allele} => $_->{a_allele}} @$matched;
  $allele_map{$vcf->get_reference} = $vf->ref_allele_string;
  my $missing_alleles = {};

  my $sample2genotype = {}; 
  SAMPLE: foreach my $sample(@$samples) {
    next unless defined($raw_gts->{$sample->{_raw_name}});
    
    my $raw_gt = $raw_gts->{$sample->{_raw_name}};
    my $phased = ($raw_gt =~ /\|/ ? 1 : 0);
    
    my @vcf_bits = split(/\||\/|\\/, $raw_gt);

    # map to VF
    my @vf_bits;
    my $missing = 0;
    foreach my $bit(@vcf_bits) {
      if(my $mapped_bit = $allele_map{$bit}) {
        push @vf_bits, $mapped_bit;
      }
      else {
        $missing_alleles->{$sample->dbID}->{$bit}++;
        $missing = 1;
      }
    }

    next SAMPLE unless @vf_bits && !$missing;

    if (!$cmp) {
      my $sgf = $self->_create_SampleGenotypeFeature($sample, \@vf_bits, $phased, $vf);
      push @genotypes, $sgf;
    } else { 
      $sample2genotype->{$sample->{name}}->{genotype} = [sort @vf_bits];
      $sample2genotype->{$sample->{name}}->{phased} = $phased;
      $sample2genotype->{$sample->{name}}->{sample} = $sample;
    }
  } 

  if (!$cmp) {
    # log the missing alleles on the first genotype object
    # the adaptor object upstream will then find this
    $genotypes[0]->{_missing_alleles} = $missing_alleles if scalar keys %$missing_alleles && @genotypes;
    return \@genotypes;
  }
  my $sample_cmp_name = $sample_cmp->name;
  my $genotype_cmp = $sample2genotype->{$sample_cmp_name}->{genotype};
  return [] if (!$genotype_cmp);
  my $phased = $sample2genotype->{$sample_cmp_name}->{phased};

  if ($cmp eq 'differences') {
    my $differences = {};
    my $reference = $vcf->get_reference;

    foreach my $sample_name (keys %$sample2genotype) {
      next if $sample_name eq $sample_cmp_name;
      my $genotype = $sample2genotype->{$sample_name}->{genotype};
      next if ($self->_same_as_reference($genotype, $reference));
      next if ($self->_compare_genotypes($genotype, $genotype_cmp));
      $differences->{$sample_name} = join('|', @$genotype);
    }
    my $sgf = $self->_create_SampleGenotypeFeature($sample_cmp, $genotype_cmp, $phased, $vf);
    $sgf->differences($differences);
    return [$sgf];
  }

  if ($cmp eq 'unique') {
    foreach my $sample_name (keys %$sample2genotype) {
      next if $sample_name eq $sample_cmp_name;
      my $genotype = $sample2genotype->{$sample_name}->{genotype};
      if ($self->_compare_genotypes($genotype, $genotype_cmp)) {
        return [];
      }
    }
    my $sgf = $self->_create_SampleGenotypeFeature($sample_cmp, $genotype_cmp, $phased, $vf);
    return [$sgf];
  }

}

sub _create_SampleGenotypeFeature {
  my $self = shift;
  my $sample = shift;
  my $bits = shift;
  my $phased = shift;
  my $vf = shift;
  return Bio::EnsEMBL::Variation::SampleGenotypeFeature->new_fast({
      _variation_id     => $vf->{_variation_id},
      _sample_id        => $sample->dbID,
      variation_feature => $vf,
      sample            => $sample,
      genotype          => $bits,
      phased            => $phased,
      adaptor           => $self->{_gta},
      start             => $vf->start,
      end               => $vf->end,
      strand            => $vf->seq_region_strand,
      slice             => $vf->slice,
  });
}

sub _same_as_reference {
  my $self = shift;
  my $genotype = shift;
  my $reference = shift;
  my $same_as_reference = 1;
  foreach my $allele (@$genotype) {
    if ($allele ne $reference) {
      $same_as_reference = 0;
    }
  }
  return $same_as_reference;
}

sub _compare_genotypes {
  my $self = shift;
  my $genotype1 = shift;
  my $genotype2 = shift;
  return 0 if (scalar(@$genotype1) != scalar(@$genotype2));

  foreach my $i (0 .. $#{$genotype1}) {
    return 0 if($genotype1->[$i] ne $genotype2->[$i]);
  }
  return 1;
}

sub _get_vcf_by_chr {
  my $self = shift;
  my $chr = shift;
  
  if(!exists($self->{files}) || !exists($self->{files}->{$chr})) {
    my $obj;
    
    # check we have this chromosome
    if(my $chrs = $self->list_chromosomes) {
      return unless grep {$chr eq $_} @$chrs;
    }
    
    my $file = $self->_get_vcf_filename_by_chr($chr);
    
    if($self->type eq 'local' && !-e $file) {
      $self->{files}->{chr} = undef;
    }
    else {

      # close first opened handle to prevent going over max allowed connections
      while(scalar @{$self->{_open_files} || []} >= $MAX_OPEN_FILES) {
        my $close = shift @{$self->{_open_files}};
        $self->{files}->{$close}->close();
        delete $self->{files}->{$close};
      }
    
      my $vcf_obj = $self->_vcf_parser_obj($file);
      if(defined $vcf_obj) {
        $self->{files}->{$chr} = $vcf_obj;
        push @{$self->{_open_files}}, $chr;
      }
    }
  }
  
  return $self->{files}->{$chr};
}

sub _get_vcf_filename_by_chr {
  my ($self, $chr) = @_;
  my $file = $self->filename_template;
  $file =~ s/\#\#\#CHR\#\#\#/$chr/;
  return $file;
}

sub _vcf_parser_obj {
  my ($self, $file) = @_;  

  # change dir
  my $cwd = cwd();
  chdir($self->tmpdir);

  # open obect (remote indexes get downloaded)
  my $obj = undef;
  eval { $obj = Bio::EnsEMBL::IO::Parser::VCF4Tabix->open($file); }; warn "$file: $@" if $@;
  # change back
  chdir($cwd);

  return $obj;
}

sub _seek {
  my $self = shift;
  my ($c, $s, $e) = @_;
  
  my $vcf = $self->_get_vcf_by_chr($c);
  return unless $vcf;
  
  # set current to the correct VCF
  $self->_current($vcf);
  
  # now seek
  $s = 0 if ($s < 0);
  $vcf->seek($c, $s, $e);

  if($self->use_seq_region_synonyms && !defined($vcf->{iterator}->{_tabix_iter}) && (my @synonyms = @{$self->_get_synonyms_by_chr($c)})) {
    while(!defined($vcf->{iterator}->{_tabix_iter}) && @synonyms) {
      $vcf->seek(shift @synonyms, $s, $e);
    }
  }
  
  return $vcf->{iterator}->{_tabix_iter} ? $vcf : undef;
}

sub _seek_by_Slice {
  my $self = shift;
  my $slice = shift;
  
  my $vcf = $self->_seek($slice->seq_region_name, $slice->start - 1, $slice->end);
  return unless $vcf;

  $vcf->next();
  
  return defined($vcf->{record});
}

sub _seek_by_VariationFeature {
  my $self = shift;
  my $vf = shift;
  
  my ($seek_start, $seek_end);
  if ($vf->var_class() eq 'SNP') {
    ## only need single location
    $seek_start = $seek_end = $vf->seq_region_start;
  } else {
    $seek_start = $vf->seq_region_start - 2;
    $seek_end =  $vf->seq_region_end + 2;
  }
  my $vcf = $self->_seek($vf->seq_region_name, $seek_start, $seek_end);

  return unless $vcf;
  
  # compare IDs
  my $count = 0;
  my $strict = $self->strict_name_match;
  my @names = ($vf->variation_name);
  push @names, @{$vf->variation->get_all_synonyms} if ($vf->variation);
  
  RECORD: while($count++ < 10 && $vcf->next()) {
    last if !$vcf->{record};
    
    # if it has an ID, we can use that
    foreach my $id(@{$vcf->get_IDs}) {
      last RECORD if grep {($_ eq $id && $id ne '.') || $_ eq 'ss'.$id} @names;
    }
    
    # otherwise compare coords
    last if !$strict && @{$self->_get_matched_alleles_VF_VCF($vf, $vcf)};
    
    # if we've gone too far, quit out
    return 0 if $vcf->get_start > $vf->seq_region_end + 1;
  }
  
  return defined($vcf->{record});
}

sub _get_matched_alleles_VF_VCF {
  my ($self, $vf, $vcf) = @_;

  return get_matched_variant_alleles(
    {
      ref    => $vf->ref_allele_string,
      alts   => $vf->alt_alleles,
      pos    => $vf->seq_region_start || $vf->start,
      strand => $vf->strand
    },
    {
      ref  => $vcf->get_reference,
      alts => $vcf->get_alternatives,
      pos  => $vcf->get_raw_start,
    }
  );
}

sub _get_Population_Sample_hash {
  my $self = shift;
  
  if(!exists($self->{_population_hash})) {
    my $hash = {};
    
    # populations defined in config?
    if(defined($self->{_raw_populations})) {
      my $pops;
      my $pa;
      my $prefix = $self->population_prefix;
      
      foreach my $sample(@{$self->get_all_Samples}) {
        foreach my $pop(@{$self->{_raw_populations}->{$sample->name} || $self->{_raw_populations}->{$sample->{_raw_name}} || []}) {
          
          # try and fetch from DB
          if(!defined($pops->{$pop})) {
            if($self->use_db) {
              $pa ||= $self->adaptor->db->get_PopulationAdaptor();
              $pops->{$pop} = $pa->fetch_by_name($prefix.$pop);
            }
          }
          
          $pops->{$pop} ||= Bio::EnsEMBL::Variation::Population->new_fast({
            name => $prefix.$pop,
            dbID => --($self->{_population_id}),
            _raw_name => $pop,
          });
          
          $hash->{$pops->{$pop}->dbID}->{$sample->dbID} = 1;
          push @{$sample->{populations}}, $pops->{$pop};
        }
        
        $self->{populations} = [values %$pops];
      }
      
      $self->{_population_hash} = $hash;
    }
    
    # otherwise we'll have to fetch from the samples
    else {
      my $samples = $self->get_all_Samples();
      
      my @dbIDs = grep {defined($_)} map {$_->dbID || undef} @$samples;
      
      my $pa = $self->adaptor->db->get_PopulationAdaptor();
      $hash = $pa->_get_sample_population_hash(\@dbIDs);
    }

    $self->{_population_hash} = $hash;
  }
  
  return $self->{_population_hash};
}

sub _add_Populations_to_Samples {
  my $self = shift;
  
  # get hash
  my $hash = $self->_get_Population_Sample_hash();
  
  # get all population objects
  my $pa = $self->adaptor->db->get_PopulationAdaptor();
  my %pop_objs_by_dbID = map {$_->dbID => $_} @{$pa->fetch_all_by_dbID_list([keys %$hash]) || []};
  
  # get sample objects
  my %sam_objs_by_dbID = map {$_->dbID => $_} @{$self->get_all_Samples()};
  
  foreach my $pop_id(keys %$hash) {
    if(my $pop_obj = $pop_objs_by_dbID{$pop_id}) {
      foreach my $sam_id(keys %{$hash->{$pop_id}}) {
        if(my $sam_obj = $sam_objs_by_dbID{$sam_id}) {
          push @{$sam_obj->{populations}}, $pop_obj;
        }
      }
    }
  }
}

sub _create_Population {
  my $self = shift;
  my $name = shift;

  my $prefix = $self->population_prefix;
  
  my ($existing) = grep {($_->{_raw_name} && $_->{_raw_name} eq $name) || $_->{name} eq $name} @{$self->{populations} || []};
  return $existing if $existing;
  
  my $pop = Bio::EnsEMBL::Variation::Population->new_fast({
    name => $prefix.$name,
    dbID => --($self->{_population_id}),
    _raw_name => $name,
  });
  
  push @{$self->{populations}}, $pop;
  
  return $pop;
}

sub _get_all_population_names {
  my $self = shift;
  
  if(!exists($self->{_population_names})) {
    
    my @names;
    
    if(defined($self->{_raw_populations})) {
      my $prefix = $self->population_prefix;
      @names =
        map {$prefix.$_}
        map {@{$self->{_raw_populations}->{$_}}}
        keys %{$self->{_raw_populations}};
    }
    
    else {
      @names = map {$_->name} @{$self->get_all_Populations};
    }
    
    $self->{_population_names} = \@names;
  }
  return $self->{_population_names};
}

=head2 vcf_collection_close

  Example    : $collection->vcf_collection_close()
  Description: Close the filehandle of the VCF collection
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub vcf_collection_close {
  my $self = shift;
  my $vcf  = $self->{current};

  $vcf->close() if ($vcf);
}

1;
