use strict;
use warnings;

package dbSNP::DBManager;

use Bio::EnsEMBL::Registry;

# A class that keeps the DBAdaptors
sub new {
  my $class = shift;
  my $registryfile = shift;
  my $species = shift;
  
  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_all( $registryfile );
  
  my $ref;
  $ref->{'registryfile'} = $registryfile;
  $ref->{'registry'} = $registry;
  $ref->{'species'} = $species;
  $ref->{'bulk_insert_buffer_size'} = (512 * 1024 * 1024);
  
  return bless($ref,$class);
}

sub dbSNP_shared {
  my $self = shift;
  my $shared = shift;
  
  if (defined($shared)) {
    $self->{'dbSNP_shared'} = $shared;
  }
  
  return $self->{'dbSNP_shared'};
}
sub dbSNP {
  my $self = shift;
  
  return $self->get_dbAdaptor('dbsnp');
}
sub dbVar {
  my $self = shift;
  
  return $self->get_dbAdaptor('variation');
}
sub dbCore {
  my $self = shift;
  
  return $self->get_dbAdaptor('core');
}

sub get_dbAdaptor {
  my $self = shift;
  my $type = shift;

  if (!defined($self->{$type})) {
    my $dba = $self->{'registry'}->get_DBAdaptor($self->{'species'},$type) or die ("Could not get DBadaptor to $type database");
    #$dba->dbc->disconnect_when_inactive(1);
    $dba->dbc->{mysql_auto_reconnect} = 1;
    
    # If we opened a new connection to the variation db, set the insert buffer size
    if ($type eq 'variation') {
      # Set some variables on the MySQL server that can speed up table read/write/loads
      my $stmt = qq{
        SET SESSION
          bulk_insert_buffer_size=$self->{'bulk_insert_buffer_size'}
      };
      $dba->dbc->do($stmt);
    }
    $self->{$type} = $dba;
  }
  
  return $self->{$type};
}

sub registryfile {
  my $self = shift;
  
  return $self->{'registryfile'};
}

sub species {
  my $self = shift;
  
  return $self->{'species'};
}

1;
