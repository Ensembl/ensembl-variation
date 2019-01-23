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

=cut
package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::Constants;

use base qw(Exporter);
use Const::Fast;

our @EXPORT_OK = qw(RGD AnimalQTL ZFIN GWAS OMIA EGA Orphanet NONE species);

use constant RGD       => 'RGD';
use constant AnimalQTL => 'AnimalQTL';
use constant ZFIN => 'ZFIN';
use constant GWAS => 'GWAS';
use constant OMIA => 'OMIA';
use constant EGA  => 'EGA';
use constant Orphanet => 'Orphanet';
use constant NONE => 'NONE';

use constant species => ( 'RGD' => ['rattus_norvegicus'],
                  'AnimalQTL' => ['gallus_gallus','equus_caballus'],
                # 'AnimalQTL' => ['bos_taurus', 'gallus_gallus', 'equus_caballus', 'sus_scrofa', 'ovis_aries'],
                 'ZFIN' => ['danio_rerio'],
            #     'OMIA' => ['felis_catus','gallus_gallus','pan_troglodytes',
            #                'bos_taurus','canis_familiaris','equus_caballus',
            #                'macaca_mulatta','sus_scrofa','ovis_aries',
            #                'meleagris_gallopavo','danio_rerio'],
                 'OMIA' => ['gallus_gallus'],
                 'GWAS' => ['homo_sapiens'],
                 'EGA'  => ['homo_sapiens'],
                 'Orphanet' => ['homo_sapiens'],
                 );
#TODO: hash key source -> array of species: used by the modules to schedule jobs for each species: see vep dum
#TODO: @: how is it best to take input params? numbers as constants or string? eg. for this _conf takes RGD, while cmd line takes 1, latest: made const string
#cow bos_taurus

#TODO: RGD: look into source/evidence for mouse and human qtls, which currently are not imported


1;
