=head1 LICENSE

 Copyright (c) 1999-2011 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

#
# Ensembl module Bio::EnsEMBL::Variation::DBSQL::VariationGroupFeatureAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::VariationGroupFeatureAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $vgfa = $reg->get_adaptor("human","variation","variationgroupfeature");
  $vga = $reg->get_adaptor("human","variation","variationgroup");
  $sa = $reg->get_adaptor("human","core","slice");

  # Get a VariationFeature by its internal identifier
  $vgf = $vgfa->fetch_by_dbID(145);

  # get all VariationFeatures in a region
  $slice = $sa->fetch_by_region('chromosome', 'X', 1e6, 2e6);
  foreach $vgf (@{$vgfa->fetch_all_by_Slice($slice)}) {
    print $vgf->start(), '-', $vgf->end(), ' ', $vgf->allele_string(), "\n";
  }


  # fetch all genome hits for a particular variation group
  $vg = $vga->fetch_by_name('PERLEGEN:B000002');

  foreach $vgf (@{$vgfa->fetch_all_by_VariationGroup($vg)}) {
    print $vgf->seq_region_name(), $vgf->seq_region_start(), '-',
          $vgf->seq_region_end(),"\n";
  }

=head1 DESCRIPTION

This adaptor provides database connectivity for VariationGroupFeature objects.
Genomic locations of VariationGroups can be obtained from the variation
database using this adaptor.  See the base class BaseFeatureAdaptor for more
information.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::VariationGroupFeatureAdaptor;

use Bio::EnsEMBL::Variation::VariationGroupFeature;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;


our @ISA = ('Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor');



=head2 fetch_all_by_VariationGroup

  Arg [1]    : Bio::EnsEMBL::Variation::Variation $var
  Example    : my @vfs = @{$vfa->fetch_all_by_VariationGroup($var)};
  Description: Retrieves all variation features for a given variation.  Most
               variations should only hit the genome once and only a return
               a single variation feature.
  Returntype : reference to list Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : At Risk

=cut


sub fetch_all_by_VariationGroup {
  my $self = shift;
  my $vg  = shift;

  if(!ref($vg) || !$vg->isa('Bio::EnsEMBL::Variation::VariationGroup')) {
    throw('Bio::EnsEMBL::Variation::VariationGroup arg expected');
  }

  if(!defined($vg->dbID())) {
    throw("VariationGroup arg must have defined dbID");
  }

  return $self->generic_fetch("vgf.variation_group_id = ".$vg->dbID());
}



# method used by superclass to construct SQL
sub _tables { return ['variation_group_feature', 'vgf']; }


sub _columns {
  return qw( vgf.variation_group_feature_id vgf.seq_region_id
             vgf.seq_region_start vgf.seq_region_end vgf.seq_region_strand
             vgf.variation_group_id vgf.variation_group_name );
}



sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;

  #
  # This code is ugly because an attempt has been made to remove as many
  # function calls as possible for speed purposes.  Thus many caches and
  # a fair bit of gymnastics is used.
  #
  my $sa = $self->db()->dnadb()->get_SliceAdaptor();

  my @features;
  my %slice_hash;
  my %sr_name_hash;
  my %sr_cs_hash;

  my ($variation_group_feature_id, $seq_region_id, $seq_region_start,
      $seq_region_end, $seq_region_strand, $variation_group_id,
      $variation_group_name );

  $sth->bind_columns(\$variation_group_feature_id, \$seq_region_id,
                     \$seq_region_start, \$seq_region_end, \$seq_region_strand,
                     \$variation_group_id, \$variation_group_name);

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

  FEATURE: while($sth->fetch()) {
    #get the slice object
    my $slice = $slice_hash{"ID:".$seq_region_id};

    if(!$slice) {
      $slice = $sa->fetch_by_seq_region_id($seq_region_id);
      $slice_hash{"ID:".$seq_region_id} = $slice;
      $sr_name_hash{$seq_region_id} = $slice->seq_region_name();
      $sr_cs_hash{$seq_region_id} = $slice->coord_system();
    }

    #
    # remap the feature coordinates to another coord system
    # if a mapper was provided
    #
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
          $seq_region_strand *= -1;
        }

        #throw away features off the end of the requested slice
        if($seq_region_end < 1 || $seq_region_start > $dest_slice_length) {
          next FEATURE;
        }
      }
      $slice = $dest_slice;
    }

    push @features, Bio::EnsEMBL::Variation::VariationGroupFeature->new
      (-start    => $seq_region_start,
       -end      => $seq_region_end,
       -strand   => $seq_region_strand,
       -slice    => $slice,
       -variation_group_name => $variation_group_name,
       -variation_group_id => $variation_group_id,
       -adaptor  => $self,
       -dbID     => $variation_group_feature_id);
  }

  return \@features;
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



1;
