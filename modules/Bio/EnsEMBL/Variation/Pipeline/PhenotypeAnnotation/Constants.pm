=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     https://www.apache.org/licenses/LICENSE-2.0

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

General constants for the PhenotypeAnnotation pipeline.

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::Constants;

use base qw(Exporter);
use Const::Fast;

our @EXPORT_OK = qw(
            RGD
            ANIMALQTL
            ZFIN
            GWAS
            OMIA
            EGA
            ORPHANET
            MIMMORBID
            DDG2P
            CGC
            IMPC
            MGI
            MOUSE
            HUMAN
            HUMAN_VAR
            HUMAN_GENE
            ANIMALSET
            GROUP_RUN_TYPES
            SOURCES_IN_RUN_TYPES
            NONE
            SPECIES);

use constant {
  RGD       => 'RGD',
  ZFIN      => 'ZFIN',

  ANIMALQTL => 'ANIMALQTL',
  OMIA      => 'OMIA',
  ANIMALSET => 'ANIMALSET',

  GWAS      => 'GWAS',
  EGA       => 'EGA',
  ORPHANET  => 'Orphanet',
  MIMMORBID => 'MIMmorbid',
  DDG2P     => 'DDG2P',
  CGC       => 'CGC',
  HUMAN     => 'HUMAN',
  HUMAN_VAR => 'HUMAN_VAR', #perform all variants only imports
  HUMAN_GENE => 'HUMAN_GENE', #perform all gene phenotype only imports

  IMPC      => 'IMPC',
  MGI       => 'MGI',
  MOUSE     => 'MOUSE',

  NONE      => 'NONE',
};

use constant GROUP_RUN_TYPES => (ANIMALSET => ['OMIA','AnimalQTL'],
                                MOUSE     => ['IMPC', 'MGI'],
                                HUMAN     => ['GWAS', 'EGA',
                                              'Orphanet', 'MIMmorbid',
                                              'DDG2P', 'CGC'],
                                HUMAN_VAR => ['GWAS', 'EGA'],
                                HUMAN_GENE => ['Orphanet', 'MIMmorbid',
                                              'DDG2P', 'CGC'],
                                        );

use constant SOURCES_IN_RUN_TYPES => ( OMIA      => 'ANIMALSET',
                                       ANIMALQTL => 'ANIMALSET',
                                       IMPC      => 'MOUSE',
                                       MGI       => 'MOUSE',

                                       GWAS      => 'HUMAN',
                                       EGA       => 'HUMAN',
                                       ORPHANET  => 'HUMAN',
                                       MIMMORBID => 'HUMAN',
                                       DDG2P     => 'HUMAN',
                                       CGC       => 'HUMAN',
                                      );

use constant SPECIES => ( 'RGD'       => ['rattus_norvegicus'],
                          'ZFIN'      => ['danio_rerio'],

                          'ANIMALQTL' => ['bos_taurus', 'gallus_gallus', 'equus_caballus',
                                          'sus_scrofa', 'ovis_aries', 'ovis_aries_rambouillet'],
                          'OMIA'      => ['felis_catus','gallus_gallus','capra_hircus',
                                          'bos_taurus','canis_lupus_familiaris','equus_caballus',
                                          'macaca_mulatta','sus_scrofa','ovis_aries', 'ovis_aries_rambouillet',
                                          'meleagris_gallopavo', 'pan_troglodytes'],

                          'GWAS'      => ['homo_sapiens'],
                          'EGA'       => ['homo_sapiens'],
                          'ORPHANET'  => ['homo_sapiens'],
                          'MIMMORBID' => ['homo_sapiens'],
                          'DDG2P'     => ['homo_sapiens'],
                          'CGC'       => ['homo_sapiens'],
                          'HUMAN'     => ['homo_sapiens'],
                          'HUMAN_VAR' => ['homo_sapiens'],
                          'HUMAN_GENE'=> ['homo_sapiens'],


                          'IMPC'      => ['mus_musculus'],
                          'MGI'       => ['mus_musculus'],

                          'ontology'  => ['homo_sapiens', 'gallus_gallus',
                                          'sus_scrofa', 'bos_taurus',
                                          'equus_caballus', 'ovis_aries', 'ovis_aries_rambouillet',
                                          'capra_hircus', 'canis_lupus_familiaris']
                        );

1;
