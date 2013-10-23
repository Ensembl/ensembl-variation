#!/usr/bin/env perl
# Copyright 2013 Ensembl
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.



=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut


#Script to generate the classes, population and sets tables in the documentation page "Data description".

use strict;
use warnings;

# Print the usage instructions if run without parameters
usage() unless (scalar(@ARGV));

my $species     = shift;
my $host        = shift;
my $db_version  = shift;
my $output_file = shift;

die ("Species, db_host and db_version must be specified") unless ($species && $host && $db_version);

open F, "> $output_file" or die $!;

# Generates the "Variation classes" table documentation
my $html = `perl generate_classes_table.pl $species $host $db_version`;
$html .= "\n\n\n\n\n\n";


# Generates the "Populations" table documentation
$html .= `perl generate_population_table.pl $species $host $db_version`;
$html .= "\n\n\n\n\n\n";


# Generates the "Variation sets" table documentation
$html .= `perl generate_variation_set_table.pl $species $host $db_version`;
$html .= "\n\n\n\n\n\n";


# Generates the "Clinical significance" tables documentation
$html .= `perl generate_clin_significance_tables.pl $species $host $db_version`;


print F $html;

sub usage {
    
    print STDOUT qq{
Usage:

  $0 SPECIES DB_HOST DB_VERSION OUTPUT_FILE_NAME
  
Description:

  Prints html code for the tables containing the available variation classes, populations and 
  variation sets for a species. The species has to be specified on the command line as the first 
  argument and the host database has to be specified as the second argument
         
};
    
    exit(0);
}
