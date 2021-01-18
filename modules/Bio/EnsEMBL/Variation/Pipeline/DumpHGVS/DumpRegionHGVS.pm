=head1 LICENSE
Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute
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

=head1 NAME

Bio::EnsEMBL::Variation::Pipeline::DumpHGVS::DumpRegionHGVS

=head1 DESCRIPTION

Dumps the HGVS for all variants in a region

=cut

package Bio::EnsEMBL::Variation::Pipeline::DumpHGVS::DumpRegionHGVS;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

use FileHandle;
use Bio::EnsEMBL::Registry;
use File::Path qw(make_path);
use Bio::EnsEMBL::Utils::Slice qw(split_Slices);


sub fetch_input {
  my $self = shift;
  my $hgvs_dir = $self->param_required('hgvs_dir');
  my $region = $self->param_required('region');
  my $bin_size = $self->param_required('bin_size');
  my $overlap = $self->param_required('region_overlap');
  $self->warning("Dumping ($region) to ($hgvs_dir)");
  $self->warning("bin_size ($bin_size) with overlap ($overlap)");
}

sub run {
  my $self = shift;
  my $hgvs_dir = $self->param('hgvs_dir');
  my $region = $self->param('region');
  my $bin_size = $self->param('bin_size');
  my $overlap = $self->param('region_overlap');

  $self->warning("region ($region)");

  my $core_dba = $self->get_species_adaptor('core');
  my $var_dba = $self->get_species_adaptor('variation');

  my $sa = $core_dba->get_SliceAdaptor();
  my $vfa = $var_dba->get_VariationFeatureAdaptor();
  $vfa->db->include_failed_variations(1);

  my $var_dbh = $var_dba->dbc->db_handle;
  my $core_dbh = $core_dba->dbc->db_handle;
  
  my ($seqname, $start, $end)  = split(/\:|\-/, $region);
  my $hgvs_file = $hgvs_dir . '/' . join('-', 'hgvs', $seqname, $start, $end) . '.tab';
  my $log_file = $hgvs_dir . '/' . join('-', 'hgvs', $seqname, $start, $end) . '.log';
  open(my $fh, '>', $hgvs_file) or die("Unable to open $hgvs_file : $!");
  open(my $lh, '>', $log_file) or die("Unable to open $log_file: $!");

  my $slice    = $sa->fetch_by_location($region, 'chromosome');

  # Split the slices into bin bin_size
  my $slice_pieces = split_Slices([$slice], $bin_size, $overlap);
  my $hgvs_count=0;

  # To do use an iterator
  for my $sub_slice (@$slice_pieces) {
    my $vfs = $vfa->fetch_all_by_Slice($sub_slice);
      foreach my $vf(@{$vfs}){
        my $hgvsg = $vf->hgvs_genomic($slice);
        for my $allele (sort keys %$hgvsg) {
            print $fh join("\t", $vf->get_Variation_dbID(), $vf->name(),
                             $vf->seq_region_name, $vf->seq_region_start(), $vf->seq_region_end(),
                             $hgvsg->{$allele}), "\n";
            $hgvs_count++;
        }
      }
  }
  $self->warning("number of records - $hgvs_count");
  print $lh "$region\t$hgvs_count\n";
  close($fh);
  close($lh);
}

sub write_output {
  my $self = shift;
}

1;
