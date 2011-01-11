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
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::AlleleGroupAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::AlleleGroupAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $aga = $reg->get_adaptor("human","variation","allelegroup");
  
  # retrieve an allele group by its name
  $ag = $aga->fetch_by_name('ABDR-10');

  # retrieve a variation group by its internal identifier
  $ag = $aga->fetch_by_dbID(123);

  # retrieve all allele groups which are part of the same variation group
  @ags = @{$aga->fetch_all_by_VariationGroup($vg)};


=head1 DESCRIPTION

This adaptor provides database connectivity for AlleleGroup objects.
AlleleGroups may be retrieved from the Ensembl variation database by
several means using this module.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::AlleleGroupAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Variation::AlleleGroup;

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

  # left join allows allele groups without any alleles to be fetched

  my $sth = $self->prepare
    (q{SELECT ag.allele_group_id, ag.variation_group_id, ag.sample_id,
              ag.name, s.name, ag.frequency, aga.allele, aga.variation_id
       FROM   (allele_group ag, source s)
       LEFT JOIN allele_group_allele aga
       ON     aga.allele_group_id = ag.allele_group_id
       WHERE  ag.source_id = s.source_id
       AND    ag.allele_group_id = ?});
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

  # left join allows allele groups without any alleles to be fetched

  my $sth = $self->prepare
    (q{SELECT ag.allele_group_id, ag.variation_group_id, ag.sample_id,
              ag.name, s.name, ag.frequency, aga.allele, aga.variation_id
       FROM   (allele_group ag, source s)
       LEFT JOIN allele_group_allele aga
       ON     aga.allele_group_id = ag.allele_group_id
       WHERE  ag.source_id = s.source_id
       AND    ag.name = ?});
  $sth->bind_param(1,$name,SQL_VARCHAR);
  $sth->execute();

  my $result = $self->_objs_from_sth($sth);
  $sth->finish();

  return undef if(!@$result);

  return $result->[0];
}




=head2 fetch_all_by_VariationGroup

  Arg [1]    : Bio::EnsEMBL::Variation::VariationGroup $vg
  Example    : @a_grps = @{$a_grp_adp->fetch_all_by_VariationGroup($var_grp)};
  Description: Retrieves all allele groups which are part of a VariationGroup
  Returntype : Bio::EnsEMBL::Variation::VariationGroup
  Exceptions : throw if VariationGroup argument is not provided
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_VariationGroup {
  my $self = shift;
  my $vg   = shift;

  if(!ref($vg) || !$vg->isa('Bio::EnsEMBL::Variation::VariationGroup')) {
    throw('Variation group argument expected');
  }

  # Add the constraint for failed variations
  my $constraint = " AND " . $self->db->exclude_failed_variations_condition();
    
  # left join allows allele groups without any alleles to be fetched

  my $sth = $self->prepare
    (qq{SELECT ag.allele_group_id, ag.variation_group_id, ag.sample_id,
              ag.name, s.name, ag.frequency, aga.allele, aga.variation_id
       FROM   (allele_group ag, source s)
       LEFT JOIN allele_group_allele aga
       ON     aga.allele_group_id = ag.allele_group_id
       LEFT JOIN failed_variation fv ON (fv.variation_id = aga.variation_id)
       WHERE  ag.source_id = s.source_id
       AND    ag.variation_group_id = ?
       $constraint
       ORDER BY ag.allele_group_id});
  $sth->bind_param(1,$vg->dbID,SQL_INTEGER);
  $sth->execute();

  my $result = $self->_objs_from_sth($sth);
  $sth->finish();

  return $result;
}



sub _objs_from_sth {
  my $self = shift;
  my $sth  = shift;

  my ($ag_id,  $vg_id, $sample_id, $name, $source, $freq, $allele, $var_id);
  $sth->bind_columns(\$ag_id, \$vg_id, \$sample_id, \$name, \$source,
                     \$freq, \$allele, \$var_id);

  my %pop_cache;
  my %var_cache;
  my %vg_cache;

  my $pop_adaptor = $self->db->get_PopulationAdaptor();
  my $vg_adaptor  = $self->db->get_VariationGroupAdaptor();

  my @results;

  my ($cur_ag, $cur_ag_id);

  # construct all variation groups without their associated variations

  while($sth->fetch()) {
    if(!defined($cur_ag) || $cur_ag_id != $ag_id) {

      # obtain the population for this allele group
      my $pop;
      if(defined($sample_id)) {
        $pop = $pop_cache{$sample_id} ||= $pop_adaptor->fetch_by_dbID($sample_id);
      }

      # obtain the variation group for this allele group
      my $vg;
      if(defined($vg_id)) {
        $vg = $vg_cache{$vg_id};
        if(!$vg) {
          $vg = $vg_adaptor->fetch_by_dbID($vg_id);

          # cache variations already obtained by variation group to
          # save going to DB again
          foreach my $v (@{$vg->get_all_Variations()}) {
            $var_cache{$v->dbID()} = $v;
          }
        }
      }

      # construct this allele group
      $cur_ag = Bio::EnsEMBL::Variation::AlleleGroup->new
        (-dbID        => $ag_id,
         -adaptor     => $self,
         -population  => $pop,
         -name        => $name,
         -source      => $source,
         -variation_group => $vg,
         -frequency   => $freq);

      $cur_ag_id = $ag_id;

      push @results, $cur_ag;
    }

    # add associated alleles and variations to this allele group
    if($var_id) {
      # all variations should be already obtained by variation group
      my $var = $var_cache{$var_id};
      if($var) {
        $cur_ag->add_Variation($var, $allele);
      } else {
        warning("AlleleGroup's VariationGroup does not contain variation " .
                "$var_id referenced by AlleleGroup.");
      }
    }
  }

  return \@results;
}




1;
