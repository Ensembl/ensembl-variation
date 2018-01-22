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
 <http://www.ensembl.org/Help/Contact>.




=head1 NAME 

Bio::EnsEMBL::Variation::Pipeline::EquivalentAlleles::FinishEquivalentAlleles

=head1 DESCRIPTION

Final module for Equivalent Allele pipeline
What checks are useful?

=cut

package Bio::EnsEMBL::Variation::Pipeline::EquivalentAlleles::FinishEquivalentAlleles;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);


sub run {
   
  my $self = shift;

  my $dir = $self->required_param('pipeline_dir');
  open my $report, ">", "$dir/EA_report.txt" || die "Failed to open EA_report.txt : $!\n";


  my $var_dba   = $self->get_species_adaptor('variation');

  my $count_attrib_ext_sth  = $var_dba->dbc->prepare(qq[ select attrib_id, count(*)
                                                         from variation_attrib
                                                        ]);

  $count_attrib_ext_sth->execute();
   
  my $data = $count_attrib_ext_sth->fetchall_arrayref();
  print $report "Attrib counts:\n";
  foreach my $l (@{$data}){
    print $report "$l->[0]\t$l->[1]\n";
  }

  $self->update_meta();

}

## store the date the job is run

sub update_meta{

  my $self = shift;

  my $var_dba  = $self->get_species_adaptor('variation');

  my $var_dbh = $var_dba->dbc->db_handle;
    
  my $update_meta_sth = $var_dbh->prepare(qq[ insert ignore into meta 
                                              ( meta_key, meta_value) values (?,?)
                                            ]);


  $update_meta_sth->execute('EquivalentAlleles_run_date', $self->run_date() );

}
1; 
 

