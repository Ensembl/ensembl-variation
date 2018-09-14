=head1 LICENSE

Copyright [2018] EMBL-European Bioinformatics Institute

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


package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportZFIN;

use strict;
use warnings;

use File::Copy;
use File::Path qw(make_path remove_tree);
use Data::Dumper;

use Bio::EnsEMBL::Registry;

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

my %source_info;
my $workdir;
my $debug;
my $inputFile;

my $core_dba;
my $variation_dba;

#-source zfin

sub fetch_input {
    #create output folder structure and fetches input files 
    my $self = shift;


sub run {
    my $self = shift;
    
    
}

sub write_output {
  }


1;

