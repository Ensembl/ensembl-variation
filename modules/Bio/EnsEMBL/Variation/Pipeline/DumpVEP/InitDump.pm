=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

  # make some lists
  my (@var, @normal, @highmem, @refseq, @refseq_highmem);

  foreach my $species(@species) {
    push @var, $species if $species->{variation};

    if($species->{species} eq 'homo_sapiens') {
      if($species->{species_refseq}) {
        push @refseq_highmem, $species;
      }
      else {
        push @highmem, $species;
      }
    }
    else {
      if($species->{species_refseq}) {
        push @refseq, $species;
      }
      else {
        push @normal, $species;
      }
    }
  }
  
  $self->param('normal_list', \@normal);
  $self->param('normal_highmem_list', \@highmem);
  $self->param('refseq_list', \@refseq);
  $self->param('refseq_highmem_list', \@refseq_highmem);
  $self->param('var_list', \@var);
  return;
}

sub write_output {
  my $self = shift;

  # 1 = distribute dumps (not set here)

  # 2 = normal
  # 3 = highmem
  # 4 = refseq
  # 5 = refseq highmem

  # 6 = all (finish_dump)
  # 7 = all refseqs (merge)
  # 8 = var (convert)

  $self->dataflow_output_id($self->param('normal_list'), 2);
  $self->dataflow_output_id($self->param('normal_highmem_list'), 3);
  $self->dataflow_output_id($self->param('refseq_list'), 4);
  $self->dataflow_output_id($self->param('refseq_highmem_list'), 5);

  $self->dataflow_output_id(
    [
      @{$self->param('normal_list')},
      @{$self->param('normal_highmem_list')},
      @{$self->param('refseq_list')},
      @{$self->param('refseq_highmem_list')}
    ],
    6
  ) if $self->param('merged');
  
  $self->dataflow_output_id(
    [
      @{$self->param('refseq_list')},
      @{$self->param('refseq_highmem_list')}
    ],
    7
  ) if $self->param('merged');

  $self->dataflow_output_id($self->param('var_list'), 8) if $self->param('convert');
  
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
  my $pattern = exists($server->{include_pattern}) ? $server->{include_pattern} : $self->param('include_pattern');
  my $exclude = exists($server->{exclude_pattern}) ? $server->{exclude_pattern} : $self->param('exclude_pattern');
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
        
        # do we have a variation DB?
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
        $species_hash{variation} = $has_var_db == 1 ? 1 : 0;
        
        # do we have SIFT or PolyPhen?
        if($has_var_db == 1) {
          $sth = $dbc->prepare("SELECT meta_key, meta_value FROM $var_db_name\.meta WHERE meta_key in ('sift_version','polyphen_version')");
          $sth->execute();
          
          my ($key, $val);
          $sth->bind_columns(\$key, \$val);
          
          while($sth->fetch) {
            $key =~ s/\_version//;
            $species_hash{$key} = 'b';
          }
          $sth->finish();
        }
        
        # do we have a regulation DB?
        my $reg_db_name = $var_db_name;
        $reg_db_name =~ s/variation/funcgen/;
        my $has_reg_db = $dbc->do("SHOW DATABASES LIKE '$reg_db_name';");
  
        my $has_reg_build = $dbc->do(qq{
          SELECT meta_value
          FROM $reg_db_name\.meta
          WHERE meta_key = 'regbuild.version'
        }) if $has_reg_db == 1;
        $species_hash{regulatory} = 1 if $has_reg_build && $has_reg_build == 1;

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
      
      # do we have a regulation DB?
      my $reg_db_name = $var_db_name;
      $reg_db_name =~ s/variation/funcgen/;
      my $has_reg_db = $dbc->do("SHOW DATABASES LIKE '$reg_db_name';");

      my $has_reg_build = $dbc->do(qq{
        SELECT meta_value
        FROM $reg_db_name\.meta
        WHERE meta_key = 'regbuild.version'
      }) if $has_reg_db == 1;
      
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
        $species_hash{variation} = $has_var_db == 1 ? 1 : 0;
        $species_hash{regulatory} = 1 if $has_reg_build && $has_reg_build == 1;
        $sth->finish();
        
        # do we have SIFT or PolyPhen?
        if($has_var_db == 1) {
          $sth = $dbc->prepare("SELECT meta_key, meta_value FROM $var_db_name\.meta WHERE meta_key in ('sift_version','polyphen_version') AND (species_id = $species_id OR species_id IS NULL)");
          $sth->execute();
          
          my ($key, $val);
          $sth->bind_columns(\$key, \$val);
          
          while($sth->fetch) {
            $key =~ s/\_version//;
            $species_hash{$key} = 'b';
          }
          $sth->finish();
        }

        push @species, \%species_hash;
        $count++;
      }
      
      die("ERROR: Problem getting species and assembly names from $current_db_name; check coord_system table\n") unless $count;
    }
  }
  
  return \@species;
}

1;

