=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
package Bio::EnsEMBL::Variation::Pipeline::DumpVEP::InitDump;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

use DBI;

my $DEBUG = 0;

sub param_defaults {
  return {
    'refseq' => 0,
    'merged' => 0,
    'convert' => 0,
    'include_pattern' => '',
    'exclude_pattern' => '',
  };
}

sub fetch_input {
  my $self = shift;

  my $servers = $self->required_param('dump_servers');
  
  my @species;
  
  foreach my $server(@$servers) {
    push @species, @{$self->get_species_list($server)};
  }
  
  $self->param('species_list', \@species);
  $self->param('refseq_list', [grep {$_->{species_refseq}} @species]);
  $self->param('variation_list', [grep {$_->{variation}} @species]);
  return;
}

sub write_output {
  my $self = shift;

  $self->dataflow_output_id($self->param('species_list'), 2);
  $self->dataflow_output_id($self->param('refseq_list'), 3) if $self->param('merged');
  $self->dataflow_output_id($self->param('variation_list'), 4) if $self->param('convert');
  return;
}

sub get_species_list {
  my $self = shift;
  my $server = shift;

  my $connection_string = sprintf(
    "DBI:mysql(RaiseError=>1):host=%s;port=%s",
    $server->{host},
    $server->{port}
  );
  
  # connect to DB
  my $dbc = DBI->connect(
      $connection_string, $server->{user}, $server->{pass}
  );
  
  my $version = $self->required_param('ensembl_release');

  my $sth = $dbc->prepare(qq{
    SHOW DATABASES LIKE '%\_core\_$version%'
  });
  $sth->execute();
  
  my $db;
  $sth->bind_columns(\$db);
  
  my @dbs;
  push @dbs, $db while $sth->fetch;
  $sth->finish;
  
  # refseq?
  if($self->param('refseq')) {
    $sth = $dbc->prepare(qq{
      SHOW DATABASES LIKE '%\_otherfeatures\_$version%'
    });
    $sth->execute();
    $sth->bind_columns(\$db);
    
    push @dbs, $db while $sth->fetch;
    $sth->finish;
  }

  # remove master and coreexpression
  @dbs = grep {$_ !~ /master|express/} @dbs;

  # filter on pattern if given
  my $pattern = $self->param('include_pattern');
  my $exclude = $self->param('exclude_pattern');
  @dbs = grep {$_ =~ /$pattern/i} @dbs if $pattern;
  @dbs = grep {$_ !~ /$exclude/i} @dbs if $exclude;

  my @species;

  foreach my $current_db_name (@dbs) {
    
    # special case otherfeatures
    # check it has refseq transcripts
    if($current_db_name =~ /otherfeatures/) {
    
      # get assembly name
      $sth = $dbc->prepare("SELECT version FROM ".$current_db_name.".coord_system ORDER BY rank LIMIT 1;");
      $sth->execute();
      my $assembly = $sth->fetchall_arrayref()->[0]->[0];
      die("ERROR: Could not get assembly name from meta table for $current_db_name\n") unless $assembly;
      
      $sth = $dbc->prepare(qq{
        SELECT COUNT(*)
        FROM $current_db_name\.transcript
        WHERE stable_id LIKE 'NM%'
        OR source = 'refseq'
      });
      $sth->execute;
      
      my $count;
      $sth->bind_columns(\$count);
      $sth->fetch;
      $sth->finish();
      
      if($count) {

        my $var_db_name = $current_db_name;
        $var_db_name =~ s/otherfeatures/variation/;
        my $has_var_db = $dbc->do("SHOW DATABASES LIKE '$var_db_name';");
        
        $current_db_name =~ s/^([a-z]+\_[a-z,1-9]+)(\_[a-z]+)?(.+)/$1$2/;
        $current_db_name =~ s/\_otherfeatures$//;
        
        # copy server details
        my %species_hash = %$server;
      
        $species_hash{species} = $current_db_name;
        $species_hash{assembly} = $assembly;
        $species_hash{species_refseq} = 1;
        $species_hash{variation} = $has_var_db;

        push @species, \%species_hash;
      }
    }
    
    else {
      # get assembly and species names
      $sth = $dbc->prepare("select species_id, meta_value from ".$current_db_name.".meta where meta_key = 'species.production_name';");
      $sth->execute();
      
      my ($species_id, $value, $species_ids);
      $sth->bind_columns(\$species_id, \$value);
      
      $species_ids->{$species_id} = $value while $sth->fetch();
      $sth->finish();
      
      my $count = 0;
      
      # do we have a variation DB?
      my $var_db_name = $current_db_name;
      $var_db_name =~ s/core/variation/;
      my $has_var_db = $dbc->do("SHOW DATABASES LIKE '$var_db_name';");
      
      foreach $species_id(keys %$species_ids) {
        $sth = $dbc->prepare("SELECT version FROM ".$current_db_name.".coord_system WHERE species_id = ".$species_id." ORDER BY rank LIMIT 1;");
        $sth->execute();
        my $assembly;
        $sth->bind_columns(\$assembly);
        $sth->execute();
        $sth->fetch();
        
        next unless $assembly;
        
        # copy server details
        my %species_hash = %$server;
        
        $species_hash{species} = $species_ids->{$species_id};
        $species_hash{assembly} = $assembly;
        $species_hash{variation} = $has_var_db;
        $sth->finish();

        push @species, \%species_hash;
        $count++;
      }
      
      die("ERROR: Problem getting species and assembly names from $current_db_name; check coord_system table\n") unless $count;
    }
  }
  
  return \@species;
}

1;

