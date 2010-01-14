#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::VariationSetAdaptor
#
# Copyright (c) 2010 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::VariationSetAdaptor

=head1 SYNOPSIS

  $db = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(...);

  $vsa = $db->get_VariationSetAdaptor();

  # retrieve a variation set by its name
  $vs = $vsa->fetch_by_name('Phenotype-associated variations');

  # retrieve a variation set by its internal identifier
  $vs = $vsa->fetch_by_dbID(12);

  # retrieve all variation sets which a variation is a part of
  @vs = @{$vsa->fetch_all_by_Variation($var)};


=head1 DESCRIPTION

This adaptor provides database connectivity for VariationSet objects.
VariationSets may be retrieved from the Ensembl variation database by
several means using this module.

=head1 AUTHOR - Pontus Larsson

=head1 CONTACT

Post questions to the Ensembl development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::VariationSetAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Variation::VariationSet;

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');


=head2 fetch_by_dbID

  Arg [1]    : int $dbID
  Example    : $vs = $vs_adaptor->fetch_by_dbID(2);
  Description: Retrieves a VariationSet object via its internal identifier.
               If no such variation set exists undef is returned.
  Returntype : Bio::EnsEMBL::Variation::VariationSet
  Exceptions : throw if dbID arg is not defined
  Caller     : general
  Status     : At Risk

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  throw('dbID argument expected') unless (defined($dbID));

  my $cols = join(',',$self->_columns());
  
  my $stmt = qq{
    SELECT
      $cols
    FROM
      variation_set vs
    WHERE
      vs.variation_set_id = ?
  };
  my $sth = $self->prepare($stmt);
  $sth->bind_param(1,$dbID,SQL_INTEGER);
  $sth->execute();

  my $result = $self->_objs_from_sth($sth);
  $sth->finish();

  return undef if(!@$result);

  return $result->[0];
}



=head2 fetch_by_name

  Arg [1]    : string $name
  Example    : $vg = $vga->fetch_by_name('Phenotype-associated variations');
  Description: Retrieves a variation set by its name. Name lookup is done with
              MySQL's LIKE function, allowing for partial matches.
  Returntype : Bio::EnsEMBL::Variation::VariationSet
  Exceptions : throw if name argument is not provided
  Caller     : general
  Status     : At Risk

=cut

sub fetch_by_name {
  my $self = shift;
  my $name = shift;

  throw('name argument expected') unless (defined($name));

  my $cols = join(',',$self->_columns());
  
  my $stmt = qq{
    SELECT
      $cols
    FROM
      variation_set vs
    WHERE
      vs.name LIKE ?
  };

  my $sth = $self->prepare($stmt);
  $sth->bind_param(1,'%' . $name . '%',SQL_VARCHAR);
  $sth->execute();

  my $result = $self->_objs_from_sth($sth);
  $sth->finish();

  return undef if(!@$result);

  return $result->[0];
}



=head2 fetch_all_by_Variation

  Arg [1]    : Bio::EnsEMBL::Variation::Variation
  Example    : my $vgs = $vga->fetch_all_by_Variation($var);
  Description: Retrieves all variation sets which a particular variation
               is present in.
  Returntype : reference to list of Bio::EnsEMBL::Variation::VariationSets
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

  my $cols = join(',',$self->_columns());
  
  my $stmt = qq{
    SELECT
      $cols
    FROM
      variation_set vs,
      variation_set_variation vsv
    WHERE
      vs.variation_set_id = vsv.variation_set_id AND
      vsv.variation_id = ?
  };

  my $sth = $self->prepare($stmt);
  $sth->bind_param(1,$var->dbID,SQL_INTEGER);
  $sth->execute();

  my $result = $self->_objs_from_sth($sth);
  $sth->finish();

  return $result;
}


sub _columns {
  return qw( vs.variation_set_id vs.name vs.description );
}
sub _tables {
  return ( ['variation_set','vs'] );
}
sub _default_where_clause {
  return '1';
}

sub _objs_from_sth {
  my $self = shift;
  my $sth  = shift;

  my ($vs_id, $name, $description);
  $sth->bind_columns(\$vs_id, \$name, \$description);

  my @results;
  my ($cur_vs, $cur_vs_id);

# Construct all variation sets

  while($sth->fetch()) {
    if (!defined($cur_vs) || $vs_id != $cur_vs_id) {
      $cur_vs = Bio::EnsEMBL::Variation::VariationSet->new
        (
          -dbID => $vs_id,
          -adaptor => $self,
          -name    => $name,
          -description    => $description
        );
      $cur_vs_id = $vs_id;
      push(@results,$cur_vs);
    }
  }

  return \@results;
}


1;
