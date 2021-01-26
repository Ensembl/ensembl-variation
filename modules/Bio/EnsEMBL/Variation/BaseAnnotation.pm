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

=head1 NAME

Bio::EnsEMBL::Variation::BaseAnnotation

=head1 SYNOPSIS

Abstract class - should not be instantiated.  Implementation of
abstract methods must be performed by subclasses.

=head1 DESCRIPTION

This class provides methods which are shared by file based annotation modules.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::BaseAnnotation;

use Cwd;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

our %TYPES = (
  'remote' => 1,
  'local'  => 1,
);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my (
    $id,
    $description,
    $type,
    $use_as_source,
    $filename_template,
    $chromosomes,
    $assembly,
    $source,
    $strict,
    $created,
    $updated,
    $is_remapped,
    $adaptor,
    $use_seq_region_synonyms,
    $tmpdir,
  ) = rearrange(
    [qw(
      ID
      DESCRIPTION
      TYPE
      USE_AS_SOURCE
      FILENAME_TEMPLATE
      CHROMOSOMES
      ASSEMBLY
      SOURCE
      STRICT_NAME_MATCH
      CREATED
      UPDATED
      IS_REMAPPED
      ADAPTOR
      USE_SEQ_REGION_SYNONYMS
      TMPDIR
    )],
    @_
  ); 
  
  throw("ERROR: No id defined for collection") unless $id;
  throw("ERROR: Collection type $type invalid") unless $type && defined($TYPES{$type});
 
  if( defined  $source && !$source->isa('Bio::EnsEMBL::Variation::Source')) {
    throw("Bio::EnsEMBL::Variation::Source argument expected");
  }

  my %collection = (
    adaptor => $adaptor,
    id => $id,
    description => $description,
    type => $type,
    use_as_source => $use_as_source,
    chromosomes => $chromosomes,
    filename_template => $filename_template,
    assembly  => $assembly,
    source => $source,
    strict_name_match => defined($strict) ? $strict : 0,
    created => $created,
    updated => $updated,
    is_remapped => $is_remapped,
    use_seq_region_synonyms => $use_seq_region_synonyms,
    tmpdir => $tmpdir || cwd(),
  );
  
  bless(\%collection, $class);
  
  return \%collection;
}


=head2 adaptor

  Arg [1]    : Bio::EnsEMBL::Variation::DBSQL::BaseAnnotationAdaptor $adaptor (optional)
               Set the adaptor for this BaseAnnotation
  Example    : my $adaptor = $collection->adaptor()
  Description: Getter/Setter for the adaptor.
  Returntype : Bio::EnsEMBL::Variation::DBSQL::BaseAnnotationAdaptor
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub adaptor {
  my $self = shift;
  $self->{adaptor} = shift if @_;
  return $self->{adaptor};
}


=head2 id

  Arg [1]    : string $id (optional)
               The new value to set the ID attribute to
  Example    : my $id = $collection->id()
  Description: Getter/Setter for the ID of this collection
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub id {
  my $self = shift;
  $self->{id} = shift if @_;
  return $self->{id};
}


=head2 description

  Arg [1]    : string $description (optional)
               The new value to set the description attribute to
  Example    : my $description = $collection->description()
  Description: Getter/Setter for the description of this collection
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub description {
  my $self = shift;
  $self->{description} = shift if @_;
  return $self->{description};
}


=head2 type

  Arg [1]    : string $type (optional)
               The new value to set the type attribute to
  Example    : my $type = $collection->type()
  Description: Getter/Setter for the type of this collection
               ('local' or 'remote').
  Returntype : string
  Exceptions : invalid type
  Caller     : general
  Status     : Stable

=cut

sub type {
  my $self = shift;
  
  if(@_) {
    my $type = shift;
    throw("ERROR: Collection type $type invalid") unless $type && defined($TYPES{$type});
    $self->{type} = shift if @_;
  }
  
  return $self->{type};
}


=head2 use_as_source

  Arg [1]    : bool $use_as_source (optional)
               The new value to set the use_as_source attribute to
  Example    : my $use_as_source = $collection->use_as_source()
  Description: Getter/Setter for the use_as_source attribute of this
               collection. Indicates to web code if we should treat
               the variants in this collection as a source/track.
  Returntype : boolean
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub use_as_source {
  my $self = shift;
  $self->{use_as_source} = shift if @_;
  return $self->{use_as_source};
}


=head2 strict_name_match

  Arg [1]    : bool $strict (optional)
               The new value to set the attribute to
  Example    : my $strict = $collection->strict_name_match()
  Description: Getter/Setter for the strict_name_match parameter. If set
               to a true value, genotypes are returned against a
               VariationFeature only if one of the names in the VCF ID field
               matches $vf->variation_name
  Returntype : boolean
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub strict_name_match {
  my $self = shift;
  $self->{strict_name_match} = shift if @_;
  return $self->{strict_name_match};
}

=head2 filename_template

  Arg [1]    : string $filename_template (optional)
               The new value to set the filename_template attribute to
  Example    : my $filename_template = $collection->filename_template()
  Description: Getter/Setter for the filename template of this collection.
               The wildcard string '###CHR###' can be used in this template
               and will be replaced with the chromosome name when reading,
               allowing a collection to consist of e.g. one VCF file per
               chromosome.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub filename_template {
  my $self = shift;
  $self->{filename_template} = shift if @_;
  return $self->{filename_template};
}


=head2 tmpdir

  Arg [1]    : string $tmpdir (optional)
               The new value to set the tmpdir attribute to
  Example    : my $tmpdir = $collection->tmpdir()
  Description: Getter/Setter for the temporary directory path used when
               downloading indexes for remote tabix files.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub tmpdir {
  my $self = shift;
  $self->{tmpdir} = shift if @_;
  return $self->{tmpdir};
}


=head2 use_seq_region_synonyms

  Arg [1]    : int $use_seq_region_synonyms (optional)
               The new value to set the use_seq_region_synonyms attribute to
  Example    : my $use_seq_region_synonyms = $collection->use_seq_region_synonyms()
  Description: Getter/Setter for the parameter that tells the API to look
               up seq_region synonyms in file based annotation based queries (e.g. 
               from VCF)
  Returntype : bool
  Caller     : general
  Status     : Stable

=cut

sub use_seq_region_synonyms {
  my $self = shift;
  $self->{use_seq_region_synonyms} = shift if @_;
  return $self->{use_seq_region_synonyms};
}

=head2 created

  arg [1]    : string $created (optional)
               the new value to set the created attribute to
  example    : my $created = $collection->created()
  description: getter/setter for created attribute of this collection
  returntype : string
  exceptions : none
  caller     : general
  status     : stable

=cut

## used for GA4GH - milliseconds from the epoch 
## could store by file (chrom) rather than collection?
sub created {
  my $self = shift;
  $self->{created} = shift if @_;
  return $self->{created};
}

=head2 updated

  arg [1]    : string $updated (optional)
               the new value to set the updated attribute to
  example    : my $updated = $collection->updated()
  description: getter/setter for updated attribute of this collection
  returntype : string
  exceptions : none
  caller     : general
  status     : stable

=cut

## used for GA4GH - milliseconds from the epoch
sub updated {
  my $self = shift;
  $self->{updated} = shift if @_;
  return $self->{updated};
}

=head2 is_remapped

  arg [1]    : boolean $is_remapped (optional)
               the new value to set the is_remapped attribute to
  example    : my $is_remapped = $collection->is_remapped()
  description: getter/setter for is_remapped attribute of this collection
  returntype : boolean
  exceptions : none
  caller     : general
  status     : stable

=cut

## info values cannot all be trusted for lifted over positions
sub is_remapped {
  my $self = shift;
  $self->{is_remapped} = shift if @_;
  return $self->{is_remapped};
}

=head2 source

  Arg [1]    : Bio::EnsEMBL::Variation::Source $src (optional)
               The new value to set the source attribute to
  Example    : $source = $collection->source()
  Description: Getter/Setter for the source object attribute
  Returntype : Bio::EnsEMBL::Variation::Source
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub source {
  my $self = shift;
  
  if(@_) {
    if(!ref($_[0]) || !$_[0]->isa('Bio::EnsEMBL::Variation::Source')) {
      throw("Bio::EnsEMBL::Variation::Source argument expected");
    }
    $self->{'source'} = shift;
  }
  
  return $self->{'source'};
}

=head2 source_name

  Arg [1]    : string $source_name (optional)
               The new value to set the source_name to
  Example    : my $source_name = $collection->source_name()
  Description: Getter/Setter for the source_name of this collection
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub source_name {
  my $self = shift;
  my $source = $self->{source};
  return unless defined $source;
  
  $source->name(@_) if(@_);
  return $source->name;
}

=head2 source_url

  Arg [1]    : string $source_url (optional)
               The new value to set the source_url to
  Example    : my $source_url = $collection->source_url()
  Description: Getter/Setter for the source_url of this collection
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub source_url {
  my $self = shift;
  my $source = $self->{source};
  return unless defined $source;

  $source->url(@_) if(@_);
  return $source->url;
}

=head2 assembly

  arg [1]    : string $assembly (optional)
               the new value to set the assembly attribute to
  example    : my $assembly = $collection->assembly()
  description: getter/setter for the assembly of this collection
  returntype : string
  exceptions : none
  caller     : general
  status     : stable

=cut

sub assembly {
  my $self = shift;
  $self->{assembly} = shift if @_;
  return $self->{assembly};
}

sub _get_synonyms_by_chr {
  my $self = shift;
  my $chr = shift;

  my $cache = $self->{_seq_region_synonyms} ||= {};

  if(!exists($cache->{$chr})) {
    my @synonyms = ();

    if(my $db = $self->adaptor->db) {
      if(my $sa = $db->dnadb->get_SliceAdaptor()) {
        if(my $s = $sa->fetch_by_region(undef, $chr)) {
          @synonyms = map {$_->name} @{$s->get_all_synonyms};
        }
      }
    }
    if (!$self->use_db) {
      push @synonyms, 'chr'.$chr;
    }
    $cache->{$chr} = \@synonyms;
  }

  return $cache->{$chr};
}

sub _current {
  my $self = shift;
  $self->{current} = shift if @_;
  return $self->{current};
}

1;
