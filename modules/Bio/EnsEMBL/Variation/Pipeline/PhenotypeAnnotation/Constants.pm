=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

=head1 Constants

General constatns for the PhenotypeAnnotation pipeline.

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::Constants;

use base qw(Exporter);
use Const::Fast;

our @EXPORT_OK = qw(
            RGD
            AnimalQTL
            ZFIN
            GWAS
            OMIA
            EGA
            Orphanet
            MIMmorbid
            DDG2P
            CGC
            IMPC
            MGI
            NONE
            species);

use constant RGD       => 'RGD';
use constant AnimalQTL => 'AnimalQTL';
use constant ZFIN      => 'ZFIN';
use constant GWAS      => 'GWAS';
use constant OMIA      => 'OMIA';
use constant EGA       => 'EGA';
use constant Orphanet  => 'Orphanet';
use constant MIMmorbid => 'MIMmorbid';
use constant DDG2P     => 'DDG2P';
use constant CGC       => 'CGC';
use constant IMPC      => 'IMPC';
use constant MGI       => 'MGI';
use constant NONE      => 'NONE';

use constant species => ( 'RGD'       => ['rattus_norvegicus'],
                          'AnimalQTL' => ['bos_taurus', 'gallus_gallus', 'equus_caballus',
                                          'sus_scrofa', 'ovis_aries'],
                          'ZFIN'      => ['danio_rerio'],
                          'OMIA'      => ['felis_catus','gallus_gallus','capra_hircus',
                                          'bos_taurus','canis_familiaris','equus_caballus',
                                          'macaca_mulatta','sus_scrofa','ovis_aries'],

                          'GWAS'      => ['homo_sapiens'],
                          'EGA'       => ['homo_sapiens'],
                          'Orphanet'  => ['homo_sapiens'],
                          'MIMmorbid' => ['homo_sapiens'],
                          'DDG2P'     => ['homo_sapiens'],
                          'CGC'       => ['homo_sapiens'],

                          'IMPC'      => ['mus_musculus'],
                          'MGI'       => ['mus_musculus'],
                          'ontology'  => ['homo_sapiens', 'gallus_gallus',
                                          'sus_scrofa', 'bos_taurus',
                                          'equus_caballus', 'ovis_aries',
                                          'capra_hircus', 'canis_familiaris']
                        );

1;
