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
package Bio::EnsEMBL::Variation::Pipeline::DumpVEP::BaseVEP;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

sub tar {
  my $self = shift;
  my $type = shift;
  my $mod  = shift;

  my $eg       = $self->param('eg'); 
  my $debug    = $self->param('debug');
  my $species  = $self->required_param('species');
  my $assembly = $self->required_param('assembly');
  my $version  = $self->required_param('ensembl_release');
  my $dir      = $self->required_param('pipeline_dir');
 
  if($eg){
     $version = $self->param_required('eg_version');      

     my $meta_container = Bio::EnsEMBL::Registry->get_adaptor($species,'core','MetaContainer');
   
     if($meta_container->is_multispecies()==1){
        my $collection_db=$1 if($meta_container->dbc->dbname()=~/(.+)\_core/);
        $dir = $dir."/".$collection_db;
     }
  } 
 
  $species .= $type ? '_'.$type : '';
  $mod ||= '';
  
  my $tar_file = sprintf(
    '%s/%s_vep_%i_%s%s.tar.gz',
    $dir,
    $species,
    $version,
    $assembly,
    $mod
  );
  
  # check if tar exists
  if(!$self->param('overwrite') && -e $tar_file) {
    print STDERR "Existing dump file found for $species, skipping (use --overwrite to overwrite)\n";
    return;
  }
  
  # check dir exists
  my $root_dir = $dir;
  my $sub_dir  = $species."/".$version."_".$assembly;
  
  die("ERROR: VEP dump directory $root_dir/$sub_dir not found") unless -e $root_dir.'/'.$sub_dir;
  
  my $command = "tar -cz -C $root_dir -f $tar_file $sub_dir";
  
  if($debug) {
    print STDERR "$command\n";
  }
  else {
    my $output = `$command`;
    die "ERROR: Failed to create tar file $tar_file\n$output\n" if $output;
  }
  
  return;
}

sub run_cmd {
  my $self = shift;
  my $cmd = shift;
  if (my $return_value = system($cmd)) {
    $return_value >>= 8;
    die "system($cmd) failed: $return_value";
  }
}

sub get_vep_params {
  my $self = shift;

  my $params = {};

  # basic params
  $params->{eg}      = $self->param('eg');
  $params->{debug}   = $self->param('debug');
  $params->{species} = $self->required_param('species');
  $params->{refseq}  = $self->required_param('species_refseq') ? '--refseq' : '';
  $params->{dir}     = $self->required_param('pipeline_dir');

  $params->{is_multispecies} = 0;

  if($params->{eg}){
     my $meta_container = Bio::EnsEMBL::Registry->get_adaptor($params->{species},'core','MetaContainer');

     if($meta_container->is_multispecies()==1){
        my $collection_db=$1 if($meta_container->dbc->dbname()=~/(.+)\_core/);
        $params->{dir} .= "/".$collection_db;
        make_path($params->{dir});

        $params->{is_multispecies} = 1;
     }

     $params->{assembly}   = $meta_container->single_value_by_key('assembly.default');
     $params->{version}    = $meta_container->schema_version();
     $params->{eg_version} = $self->param('eg_version');
     $params->{host}       = $meta_container->dbc->host();
     $params->{port}       = $meta_container->dbc->port();
     $params->{user}       = $meta_container->dbc->username();
     $params->{pass}       = $meta_container->dbc->password() ? '--pass '.$meta_container->dbc->password() : '';

     $meta_container->dbc()->disconnect_if_idle();
     
     $self->param('assembly', $params->{assembly});
     $self->param('ensembl_release', $params->{version});
  }
  else {
     $params->{assembly} = $self->required_param('assembly');
     $params->{version}  = $self->required_param('ensembl_release');
     $params->{host}     = $self->required_param('host');
     $params->{port}     = $self->required_param('port');
     $params->{user}     = $self->required_param('user');
     $params->{pass}     = $self->required_param('pass') ? '--pass '.$self->required_param('pass') : '';
  }

  # species-specific
  my $species_flags = $self->param('species_flags');
  
  # copy in sift, polyphen, regulatory
  $species_flags->{$params->{species}}->{$_} = $self->param($_) for grep {$self->param($_)} qw(sift polyphen regulatory);
  
  my $species_flags_cmd = $params->{refseq}.' ';
  if(my $flags = $species_flags->{$params->{species}}) {
    
    # assembly-specific
    if(my $as = $flags->{assembly_specific}) {
      delete $flags->{assembly_specific};

      my $assembly = $params->{assembly};
      
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

  $params->{species_flags_cmd} = $species_flags_cmd;

  # cmnd line
  $params->{perl}    = $self->required_param('perl_command');
  $params->{vep_dir} = $self->required_param('ensembl_cvs_root_dir').'/ensembl-tools/scripts/variant_effect_predictor';
  $params->{vep}     = $self->required_param('vep_command');

  return $params;
}

sub healthcheck_cache {
  my $self = shift;
  my $params = shift;

  # healthcheck resultant cache
  my $script_dir = $self->required_param('ensembl_cvs_root_dir').'/ensembl-variation/scripts/misc';
  
  # we use Test::Harness as the test script itself will run many 1000's of tests
  # Test::Harness gives a nice summary output instead of splurging everything to STDOUT/ERR
  # to run this under Test::Harness we need to set ENV variables
  $ENV{HC_VEP_HOST}     = $params->{host};
  $ENV{HC_VEP_PORT}     = $params->{port};
  $ENV{HC_VEP_USER}     = $params->{user};
  $ENV{HC_VEP_SPECIES}  = $params->{species}.($params->{refseq} ? '_refseq' : '');
  $ENV{HC_VEP_VERSION}  = $params->{version};
  $ENV{HC_VEP_DIR}      = $params->{dir};
  $ENV{HC_VEP_NO_FASTA} = 1;
  $ENV{HC_VEP_MAX_VARS} = 100;
  $ENV{HC_VEP_RANDOM}   = $self->param('hc_random');
  
  my $cmd = sprintf(
    '%s -MTest::Harness -e"runtests(@ARGV)" %s/healthcheck_vep_caches.pl',
    $params->{perl},
    $script_dir
  );
  
  my $finished = 0;

  open CMD, "$cmd 2>&1 |" or die "ERROR: Failed to run command $cmd";

  my $pipedir = $self->required_param('pipeline_dir');
  my $failed = 0;

  my $hc_file = sprintf("%s/%s_%s_%s_QC_report.txt", $pipedir, map {$params->{$_}} qw(species version assembly));

  open REPORT, ">>", $hc_file or die "Failed to open $hc_file : $!\n";
  while(<CMD>) { print REPORT; $failed = 1 if /fail/i; }
  close CMD;
  close REPORT;

  die("ERROR: Cache healthcheck failed, see $hc_file for details\n") if $failed;
}
1;
