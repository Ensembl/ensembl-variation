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
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::VariationGroupAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::VariationGroupAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $vga = $reg->get_adaptor("human","variation","variationgroup");
  
  # retrieve a variation group by its name
  $vg = $vga->fetch_by_name('PERLEGEN:B000009');

  # retrieve a variation group by its internal identifier
  $vg = $vga->fetch_by_dbID(63211);

  # retrieve all variation groups which a variation is a part of
  @vgs = @{$vga->fetch_all_by_Variation($var)};


=head1 DESCRIPTION

This adaptor provides database connectivity for VariationGroup objects.
VariationGroups may be retrieved from the Ensembl variation database by
several means using this module.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::VariationGroupAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Variation::VariationGroup;

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');




=head2 fetch_by_dbID

  Arg [1]    : int $dbID
  Example    : $vg = $vg_adaptor->fetch_by_dbID(5526);
  Description: Retrieves a VariationGroup object via its internal identifier.
               If no such variation group exists undef is returned.
  Returntype : Bio::EnsEMBL::Variation::VariationGroup
  Exceptions : throw if dbID arg is not defined
  Caller     : general
  Status     : At Risk

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  throw('dbID argument expected') if(!defined($dbID));

  # left join allows variation groups without any variations to be fetched

  my $sth = $self->prepare
    (q{SELECT vg.variation_group_id, vg.name, s.name, vg.type, vgv.variation_id
       FROM   (variation_group vg, source s)
       LEFT JOIN variation_group_variation vgv ON
                 vgv.variation_group_id = vg.variation_group_id
       WHERE  vg.source_id = s.source_id
       AND    vg.variation_group_id = ?});
  $sth->bind_param(1,$dbID,SQL_INTEGER);
  $sth->execute();

  my $result = $self->_objs_from_sth($sth);
  $sth->finish();

  return undef if(!@$result);

  return $result->[0];
}



=head2 fetch_by_name

  Arg [1]    : string $name
  Example    : $vg = $vga->fetch_by_name('PERLEGEN:B000009');
  Description: Retrieves a variation group by its name
  Returntype : Bio::EnsEMBL::Variation::VariationGroup
  Exceptions : throw if name argument is not provided
  Caller     : general
  Status     : At Risk

=cut

sub fetch_by_name {
  my $self = shift;
  my $name = shift;

  throw('name argument expected') if(!defined($name));

  # left join allows variation groups without any variations to be fetched

  my $sth = $self->prepare
    (q{SELECT vg.variation_group_id, vg.name, s.name, vg.type, vgv.variation_id
       FROM   (variation_group vg, source s)
       LEFT JOIN variation_group_variation vgv
       ON    vgv.variation_group_id = vg.variation_group_id
       WHERE  vg.source_id = s.source_id
       AND    vg.name = ?});
  $sth->bind_param(1,$name,SQL_VARCHAR);
  $sth->execute();

  my $result = $self->_objs_from_sth($sth);
  $sth->finish();

  return undef if(!@$result);

  return $result->[0];
}



=head2 fetch_all_by_Variation

  Arg [1]    : Bio::EnsEMBL::Variation::Variation
  Example    : my $vgs = $vga->fetch_all_by_Variation($var);
  Description: Retrieves all variation groups which a particular variation
               is present in.
  Returntype : reference to list of Bio::EnsEMBL::Variation::VariationGroups
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_Variation {
  my $self = shift;
  my $var  = shift;

  if(!ref($var) || !$var->isa('Bio::EnsEMBL::Variation::Variation')) {
    throw('Bio::EnsEMBL::Variation::Variation argument expected');
  }

  my $sth = $self->prepare
    (q{SELECT vg.variation_group_id, vg.name, s.name, vg.type, vgv.variation_id
       FROM   variation_group vg, source s, variation_group_variation vgv
       WHERE  vg.source_id = s.source_id
       AND    vgv.variation_group_id = vg.variation_group_id
       AND    vgv.variation_id = ?
       ORDER BY vg.variation_group_id});
  $sth->bind_param(1,$var->dbID,SQL_INTEGER);
  $sth->execute();

  my $result = $self->_objs_from_sth($sth);
  $sth->finish();

  return $result;
}




sub _objs_from_sth {
  my $self = shift;
  my $sth  = shift;

  my ($vg_id, $name, $source, $type, $var_id);
  $sth->bind_columns(\$vg_id, \$name, \$source, \$type, \$var_id);

  my %seen_vars;

  my @results;

  my ($cur_vg, $cur_vg_id);

  # construct all variation groups without their associated variations

  while($sth->fetch()) {
    if(!defined($cur_vg) || $cur_vg_id != $vg_id) {
      $cur_vg = Bio::EnsEMBL::Variation::VariationGroup->new
        (-dbID => $vg_id,
         -adaptor => $self,
         -name    => $name,
         -type    => $type,
         -source  => $source);
      $cur_vg_id = $vg_id;
      push @results, $cur_vg;
    }

    if(defined($var_id)) {
      $seen_vars{$var_id} ||= [];
      push @{$seen_vars{$var_id}}, $cur_vg;
    }
  }

  # fetch all of the variations at once and add them to the variation groups

  my @var_ids = keys %seen_vars;
  my $va = $self->db->get_VariationAdaptor();
  my $vars = $va->fetch_all_by_dbID_list(\@var_ids);

  foreach my $var (@$vars) {
    foreach my $vg (@{$seen_vars{$var->dbID()}}) {
      $vg->add_Variation($var);
    }
  }

  return \@results;
}




1;
