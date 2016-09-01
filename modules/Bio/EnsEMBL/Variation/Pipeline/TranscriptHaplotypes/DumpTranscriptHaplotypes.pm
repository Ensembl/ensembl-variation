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

package Bio::EnsEMBL::Variation::Pipeline::TranscriptHaplotypes::DumpTranscriptHaplotypes;

use Digest::MD5 qw(md5_hex);
use JSON;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

my $DEBUG = 0;

sub run {
  my $self = shift;

  delete($self->{_write_dir}) if $self->{_write_dir};

  my $tr_id = $self->required_param('transcript_stable_id');

  my $var_dba = $self->get_species_adaptor('variation');
  my $core_dba = $self->get_species_adaptor('core');

  my $ta = $core_dba->get_TranscriptAdaptor;
  my $tr = $ta->fetch_by_stable_id($tr_id) or die "ERROR: Failed to fetch transcript $tr_id\n";

  if(my $json_config = $self->param('json_config')) {
    $var_dba->vcf_config_file($json_config);
  }

  my $filter;
  if(my $freq = $self->param('filter_frequency')) {
    $filter = {frequency => {frequency => $freq}};
  }
  
  my $tha = $var_dba->get_TranscriptHaplotypeAdaptor;
  my $thc = $tha->get_TranscriptHaplotypeContainer_by_Transcript($tr, $filter);
  $thc->total_expected_frequency_delta;
  $thc->get_all_total_expected_population_frequency_deltas;

  $self->dump_ph_freqs($thc);
  $self->dump_ch_freqs($thc);
  $self->dump_pd_freqs($thc);
  $self->dump_cd_freqs($thc);

  $self->dump_diff_freqs($thc);

  $self->dump_to_json($thc);

  return;
}

sub dump_to_json {
  my $self = shift;
  my $thc = shift;

  my $dir = $self->get_write_dir();

  my $tr_id = $self->required_param('transcript_stable_id');

  my $file = sprintf('%s/%s.haplotypes.json.gz', $dir, $tr_id);
  open OUT, "| gzip -9 -c > $file" or die "ERROR: Could not open file $file for writing\n";

  my $json = JSON->new;
  print OUT $json->allow_blessed->convert_blessed->encode($thc);
  close OUT;
}

sub dump_ph_freqs {
  my $self = shift;
  return $self->_generic_dump_freqs(@_, 'Protein', 'Ha');
}

sub dump_ch_freqs {
  my $self = shift;
  return $self->_generic_dump_freqs(@_, 'CDS', 'Ha');
}

sub dump_pd_freqs {
  my $self = shift;
  return $self->_generic_dump_freqs(@_, 'Protein', 'Di');
}

sub dump_cd_freqs {
  my $self = shift;
  return $self->_generic_dump_freqs(@_, 'CDS', 'Di');
}

sub _generic_dump_freqs {
  my ($self, $thc, $type, $haplo_diplo) = @_;

  my $dir = $self->get_write_dir();

  my $tr_id = $self->required_param('transcript_stable_id');

  my $file = sprintf('%s/%s.%s_%splotypes.txt', $dir, $tr_id, lc($type), lc($haplo_diplo));
  open OUT, ">$file" or die "ERROR: Could not open file $file for writing\n";

  my $method = 'get_all_'.$type.$haplo_diplo.'plotypes';

  foreach my $ph(sort {$b->count <=> $a->count} @{$thc->$method}) {
    my $f  = $ph->get_all_population_frequencies;

    my $flags = $ph->can('get_all_flags') ? join(",", @{$ph->get_all_flags}) : '.';

    printf OUT "%s\t%s\t%s\t%s\t%i\t%.3g\t%.3g\t%.3g\t%.3g\t%.3g\t%.3g",
      $ph->_hex,
      $tr_id,
      $ph->name,
      $flags || '.',
      $ph->count,
      $ph->frequency,
      ($f->{AFR} || $f->{'1000GENOMES:phase_3:AFR'} || 0),
      ($f->{AMR} || $f->{'1000GENOMES:phase_3:AMR'} || 0),
      ($f->{EAS} || $f->{'1000GENOMES:phase_3:EAS'} || 0),
      ($f->{EUR} || $f->{'1000GENOMES:phase_3:EUR'} || 0),
      ($f->{SAS} || $f->{'1000GENOMES:phase_3:SAS'} || 0);

    if($ph->can('expected_frequency')) {
      my $ef = $ph->get_all_expected_population_frequencies;

      printf OUT "\t%.3g\t%.3g\t%.3g\t%.3g\t%.3g\t%.3g",
        $ph->expected_frequency,
        ($ef->{AFR} || $ef->{'1000GENOMES:phase_3:AFR'} || 0),
        ($ef->{AMR} || $ef->{'1000GENOMES:phase_3:AMR'} || 0),
        ($ef->{EAS} || $ef->{'1000GENOMES:phase_3:EAS'} || 0),
        ($ef->{EUR} || $ef->{'1000GENOMES:phase_3:EUR'} || 0),
        ($ef->{SAS} || $ef->{'1000GENOMES:phase_3:SAS'} || 0);
    }

    print OUT "\n";
  }

  close OUT;
}

sub dump_diff_freqs {
  my ($self, $thc) = @_;

  my $dir = $self->get_write_dir();

  my $tr_id = $self->required_param('transcript_stable_id');

  my $file = sprintf('%s/%s.diff_freqs.txt', $dir, $tr_id);
  open OUT, ">$file" or die "ERROR: Could not open file $file for writing\n";

  my $freqs = $thc->_protein_allele_frequencies();

  foreach my $pos(sort {$a <=> $b} keys %$freqs) {
    foreach my $a(keys %{$freqs->{$pos}}) {
      printf OUT "%i\t%s\t%g\t%g\t%g\t%g\t%g\t%g\n",
        $pos + 1, $a,
        $freqs->{$pos}->{$a}->{'_all'} || 0,
        $freqs->{$pos}->{$a}->{'1000GENOMES:phase_3:AFR'} || 0,
        $freqs->{$pos}->{$a}->{'1000GENOMES:phase_3:AMR'} || 0,
        $freqs->{$pos}->{$a}->{'1000GENOMES:phase_3:EAS'} || 0,
        $freqs->{$pos}->{$a}->{'1000GENOMES:phase_3:EUR'} || 0,
        $freqs->{$pos}->{$a}->{'1000GENOMES:phase_3:SAS'} || 0,
    }
  }

  close OUT;
}

sub get_write_dir {
  my $self = shift;

  my $tr_id = $self->required_param('transcript_stable_id');

  my $dir = $self->required_param('pipeline_dir').'/haplotypes/'.substr(md5_hex($tr_id), 0, 2);

  unless(-d $dir) {
    mkdir($dir) or die "ERROR: Could not make directory $dir\n";
  }

  return $dir;
}

1;
