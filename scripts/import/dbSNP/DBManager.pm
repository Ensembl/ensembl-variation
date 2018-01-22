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
  <helpdesk.org>.

=cut

use strict;
use warnings;

package dbSNP::DBManager;

use Bio::EnsEMBL::Registry;

# A class that keeps the DBAdaptors
sub new {
  my $class = shift;
  my $registryfile = shift;
  my $species = shift;
  my $schema = shift; ##for postgresSQL
  
  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_all( $registryfile );
  
  my $ref;
  $ref->{'registryfile'} = $registryfile;
  $ref->{'registry'} = $registry;
  $ref->{'species'} = $species;
  $ref->{'bulk_insert_buffer_size'} = (512 * 1024 * 1024);
  $ref->{'schema_name'}  = $schema if defined $schema;    

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
  
  my $adaptor =  $self->get_dbAdaptor('dbsnp');

  if(defined $self->{'schema_name'} && $self->{'schema_name'} =~/\w+/){
      print "setting search path for postgreSQL: " .$self->{'schema_name'} . ", " . $self->{'dbSNP_shared'} ."\n";
      my $sth = "SET search_path TO $self->{'schema_name'},$self->{'dbSNP_shared'},public";
      $adaptor->dbc()->do($sth);
  }
  return $adaptor;


}
sub dbVar {
  my $self = shift;
  
  return $self->get_dbAdaptor('variation');
}
sub dbCore {
  my $self = shift;
  
  return $self->get_dbAdaptor('core');
}
sub dbInt {
  my $self = shift;

  return $self->get_dbAdaptor('intvar', 'multi');
}
sub get_dbAdaptor {
  my $self = shift;
  my $type = shift;
  my $species = shift || $self->{'species'};

  if (!defined($self->{$type})) {
    my $dba = $self->{'registry'}->get_DBAdaptor($species,$type) or die ("Could not get DBadaptor to $type database");
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
sub registry {
  my $self = shift;
  
  return $self->{'registry'};
}

sub species {
  my $self = shift;
  
  return $self->{'species'};
}

## for postgreSQL
sub schema {
  my $self = shift;

  return $self->{'schema'};
}


1;
