=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
package Bio::EnsEMBL::Variation::Pipeline::DumpVEP::DumpVEP;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::DumpVEP::BaseVEP);

sub param_defaults {
  return {
    'species_refseq' => 0,
    'species_flags'  => {},
    'overwrite'      => 0,
  };
}


sub run {
  my $self = shift;

  # basic params
  my $eg      = $self->param('eg');
  my $debug   = $self->param('debug');
  my $species = $self->required_param('species');
  my $refseq  = $self->required_param('species_refseq') ? '--refseq' : '';
  my $dir     = $self->required_param('pipeline_dir');

  my ($assembly, $version);
  my ($host, $port, $user, $pass);

  if($eg){
     my $meta_container = Bio::EnsEMBL::Registry->get_adaptor($species,'core','MetaContainer');

     $assembly = $meta_container->single_value_by_key('assembly.default');
     $version  = $meta_container->schema_version();
     $host     = $meta_container->dbc->host();
     $port     = $meta_container->dbc->port();
     $user     = $meta_container->dbc->username();
     $pass     = $meta_container->dbc->password() ? '--pass '.$meta_container->dbc->password() : '';
     
     $self->param('assembly', $assembly);
     $self->param('ensembl_release', $version);
  }
  else {
     $assembly = $self->required_param('assembly');
     $version  = $self->required_param('ensembl_release');
     $host = $self->required_param('host');
     $port = $self->required_param('port');
     $user = $self->required_param('user');
     $pass = $self->required_param('pass') ? '--pass '.$self->required_param('pass') : '';
  }

  # species-specific
  my $species_flags = $self->param('species_flags');
  
  my $species_flags_cmd = $refseq.' ';
  if(my $flags = $species_flags->{$species}) {
    
    # assembly-specific
    if(my $as = $flags->{assembly_specific}) {
      delete $flags->{assembly_specific};
      
      if($as->{$assembly}) {
        
        foreach my $key(keys %{$as->{$assembly}}) {
          my $v = $as->{$assembly}->{$key};
          
          if(ref($v) eq 'ARRAY') {
            $species_flags_cmd .= sprintf(' --%s %s ', $key, $_ eq '1' ? '' : $_) for @{$v};
          }
          
          else {
            $species_flags_cmd .= sprintf(' --%s %s ', $key, $v eq '1' ? '' : $_);
          }
        }
      }
    }
    
    $species_flags_cmd .= join(' ', map {$flags->{$_} eq '1' ? '--'.$_ : '--'.$_.' '.$flags->{$_}} keys %$flags);
  }

  # cmnd line
  my $perl    = $self->required_param('perl_command');
  my $vep_dir = $self->required_param('ensembl_cvs_root_dir').'/ensembl-tools/scripts/variant_effect_predictor';
  my $vep     = $self->required_param('vep_command');
  
  # construct command
  my $cmd = sprintf(
    '%s %s/variant_effect_predictor.pl %s --host %s --port %i --user %s %s --species %s --assembly %s --db_version %s --dir %s %s',
    $perl,
    $vep_dir,
    $vep,
    
    $host,
    $port,
    $user,
    $pass,
    
    $species,
    $assembly,
    $version,
    $dir,
    $species_flags_cmd
  );
  
  my $finished = 0;
  
  if($debug) {
    print STDERR "$cmd\n";
  }
  else {
    open CMD, "$cmd 2>&1 |" or die "ERROR: Failed to run command $cmd";
    my @buffer;
    while(<CMD>) {
      $finished = 1 if /Finished/;
      push @buffer, $_;
      shift @buffer if scalar @buffer > 5;
    }
    close CMD;
  
    die "ERROR: Encountered an error running VEP\n".join("", @buffer)."\n" unless $finished;
  
    # healthcheck resultant cache
    my $script_dir = $self->required_param('ensembl_cvs_root_dir').'/ensembl-variation/scripts/misc';
    
    # we use Test::Harness as the test script itself will run many 1000's of tests
    # Test::Harness gives a nice summary output instead of splurging everything to STDOUT/ERR
    # to run this under Test::Harness we need to set ENV variables
    $ENV{HC_VEP_HOST}     = $host;
    $ENV{HC_VEP_PORT}     = $port;
    $ENV{HC_VEP_USER}     = $user;
    $ENV{HC_VEP_SPECIES}  = $species.($refseq ? '_refseq' : '');
    $ENV{HC_VEP_VERSION}  = $version;
    $ENV{HC_VEP_DIR}      = $dir;
    $ENV{HC_VEP_NO_FASTA} = 1;
    $ENV{HC_VEP_MAX_VARS} = 100;
    $ENV{HC_VEP_RANDOM}   = $self->param('hc_random');
    
    $cmd = sprintf(
      '%s -MTest::Harness -e"runtests(@ARGV)" %s/healthcheck_vep_caches.pl',
      $perl,
      $script_dir
    );
    
    $finished = 0;
  
    open CMD, "$cmd 2>&1 |" or die "ERROR: Failed to run command $cmd";
  
    my $pipedir = $self->required_param('pipeline_dir');

    open REPORT, ">>", "$pipedir/$species\_$version\_$assembly\_QC_report.txt" or die "Failed to open $pipedir/$species\_$version\_$assembly\_QC_report.txt : $!\n";
    while(<CMD>) { print REPORT; }
    close CMD;
    close REPORT;
  }
  
  $self->tar($self->param('species_refseq') ? 'refseq' : '');
  
  return;
}


1;
