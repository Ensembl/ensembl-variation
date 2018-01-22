#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut


use strict;
#use DBH;
use Bio::EnsEMBL::Registry;
use FindBin qw( $Bin );
use Getopt::Long;
use ImportUtils qw(dumpSQL debug load);


our ($input_file, $species, $TMP_DIR, $TMP_FILE);

GetOptions('species=s'    => \$species,
	   'input_file=s' => \$input_file,
	   'tmpdir=s'     => \$ImportUtils::TMP_DIR,
	   'tmpfile=s'    => \$ImportUtils::TMP_FILE
	  );

usage('-species argument is required') if(!$species);

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

my $registry_file;
$registry_file ||= $Bin . "/ensembl.registry";


Bio::EnsEMBL::Registry->load_all( $registry_file );

my $dbCore = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $dbVariation = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');

my $variation_adaptor = $dbVariation->get_VariationAdaptor;
my $gene_adaptor = $dbCore->get_GeneAdaptor;
my $dbentry_adaptor = $dbCore->get_DBEntryAdaptor();

open IN, "<$input_file" || die "Could not open inout file to parse: $!\n";
#let's parse the file
<IN>; #ignore first line, the header
my ($variation_name,$gene_xref,$p_value);
my ($variation,$genes);
my @row;
my $buffer = {}; #buffer to store the lines to be written
while (<IN>){
    chomp;
    @row = split;
    $variation_name = $row[0];
    $gene_xref = $row[2];
    $p_value = $row[8];
    $variation = $variation_adaptor->fetch_by_name($variation_name);
    if (!defined $variation){
	print "Variation $variation_name not in the database\n";
	next;
    }
    $genes = $gene_adaptor->fetch_all_by_external_name($gene_xref);
    if (!defined $genes->[0]){
	print "Gene $gene_xref not in the database\n";
	next;
    }
    print_buffered($buffer,"$TMP_DIR/$TMP_FILE",join("\t",$variation->dbID,$genes->[0]->dbID,"population",$p_value,'gene_expression',"citation") . "\n");

}
print_buffered($buffer);

close IN || die "Could not close file with association data: $!\n";
#import the file
load($dbVariation->dbc,"association","variation_id","gene_id","population","p_value","association_type","citation");

sub print_buffered {
    my $buffer = shift;
    my $filename = shift;
    my $text = shift;

    local *FH;

    if( ! $filename ) {
	# flush the buffer
	foreach my $file (keys %{$buffer}){
	    open( FH, ">>$file" ) or die;
	    print FH $buffer->{ $file };
	    close FH;
	}

    } else {
	$buffer->{ $filename } .= $text;
	if( length( $buffer->{ $filename } ) > 10_000 ) {
	    open( FH, ">>$filename" ) or die;
	    print FH $buffer->{ $filename };
	    close FH;
	    $buffer->{ $filename } = '';
	}
    }
}
