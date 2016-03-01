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
package Bio::EnsEMBL::Variation::Pipeline::DumpVEP::DumpVEP;

use strict;
use warnings;
use File::Path qw(make_path);

use base qw(Bio::EnsEMBL::Variation::Pipeline::DumpVEP::BaseVEP);

sub param_defaults {
  return {
    'species_refseq' => 0,
    'species_flags'  => {},
    'overwrite'      => 0,
    'sift'           => 0,
    'polyphen'       => 0,
    'regulatory'     => 0,
    'eg'             => 0,
  };
}


sub run {
  my $self = shift;

  my $params = $self->get_vep_params();
  
  # construct command
  my $cmd = sprintf(
    '%s %s/variant_effect_predictor.pl %s --host %s --port %i --user %s %s --species %s --assembly %s --db_version %s --dir %s %s --cache_version %s --is_multispecies %s',
    $params->{perl},
    $params->{vep_dir},
    $params->{vep},
    
    $params->{host},
    $params->{port},
    $params->{user},
    $params->{pass},
    
    $params->{species},
    $params->{assembly},
    $params->{version},
    $params->{dir},
    $params->{species_flags_cmd},

    $params->{eg_version} || $params->{version},
    $params->{is_multispecies}
  );
 
  my $finished = 0;
 
  if($params->{debug}) {
    print STDERR "$cmd\n";
  }
  else {
    open CMD, "$cmd 2>&1 |" or die "ERROR: Failed to run command $cmd";
    my @buffer;
    while(<CMD>) {
      $finished = 1 if /Finished/;
      push @buffer, $_;
      shift @buffer if scalar @buffer > 20;
    }
    close CMD;
  
    die "ERROR: Encountered an error running VEP\n".join("", @buffer)."\n" unless $finished;
  
    # healthcheck resultant cache
    $self->healthcheck_cache($params);
  
    $self->tar($self->param('species_refseq') ? 'refseq' : '');
  }
  
  return;
}


1;
