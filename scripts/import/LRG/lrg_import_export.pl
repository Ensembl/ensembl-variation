#!/usr/local/ensembl/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);
use Writer;
use IO::File;

#use Bio::EnsEMBL::Registry;
#use Bio::EnsEMBL::Mapper::RangeRegistry;
#use Bio::EnsEMBL::Variation::Utils::Sequence qw(ambiguity_code);

my ($species, $input_file,$xml_file);


GetOptions('xml_file=s' => \$xml_file,
           'species=s'   => \$species,
           'input_file=s' => \$input_file,
           );

print "opening $xml_file\n";

# ----------- XML settings ---------------------------------------#
my $infh = new IO::File("<${input_file}");
my $output = new IO::File(">${xml_file}");

my $writer = new XML::Writer(OUTPUT      => $output,
                             DATA_MODE   => 'true',
                             DATA_INDENT => 2,
			     #NEWLINES    => 'true',
                             );

while (<$infh>) {}
# -----------  XML printing ---------------------------------------#

my $date = "Today";
my $dna_sequence = "ACTG";
my @transcripts = ("t1","t2");
my $exon1 = {};
my $exon2 = {};
$exon1->{start} = 1;
$exon1->{end}  = 100;
$exon2->{start} = 200;
$exon2->{end} = 300;

my @exons;

push @exons,$exon1,$exon2;

my ($trans_count,$exon_count,$cds_count);
$trans_count = $exon_count = $cds_count = 0;
my $cds_start = 1;
my $cds_end = 100;
my $cds_sequence = "GYC";
my $cdna_sequence = "actg";
my $modification_date = "06/02/09";
my $source_name = "Osteogenesis Imperfecta Mutation Database";
my $url = "oi.gene.le.ac.uk/home.php?select_db=COL1A1";
my $contact_name = "Raymond Dalgleish";
my $contact_address = "Department of Genetics, University Road, LE1 7RH, United Kingdom";
my $email = "raymond.dalgleish\@le.ac.uk";

my $coding = 1;
my $cdna = 1;
my $translation=1;
my $other_exon_naming=1;
my $exon_url = "exon_url";
my $mapping_assembly=1;
my $amino_acid_url = "amino_acid_url";
my $amino_acid_mapping = 1;

$writer->xmlDecl('UTF-8');
$writer->pi('xml-stylesheet', 'type="text/xsl href="lrg2html.xsl"');
$writer->startTag('lrg');

$writer->startTag('fixed_annotation');
$writer->dataElement('id','LRG_1');
$writer->dataElement('organism', 'Homo sapiens', 'taxon' => '9606');


print_source_tag($source_name,$url,$contact_name,$contact_address,$email);

$writer->dataElement('mol_type','genomic DNA');

$writer->dataElement('creation_data',$date);

$writer->dataElement('sequence',$dna_sequence);


foreach my $transcript (@transcripts) {
  $writer->startTag('transcript', 'name' => "t".$trans_count);
  $trans_count++;
  foreach my $exon (sort{$a->{start}<=>$b->{start}} @exons) {
    $exon_count++;
    my $start = $exon->{start};
    my $end = $exon->{end};
    $writer->startTag('exon', 'lrg_number' => $exon_count);
    $writer->dataElement('start',$start);
    $writer->dataElement('end',$end);
    $writer->endTag('exon');
  }
  if ($coding) {
    $writer->startTag('coding_region');
    $writer->dataElement('cds_start',$cds_start);
    $writer->dataElement('cds_end',$cds_end);
    if ($translation) {
      $writer->startTag('translation');
      $writer->dataElement('sequence', $cds_sequence);
      $writer->endTag('translation');
    }
    $writer->endTag('coding_region');
  }

  if ($cdna) {
    $writer->startTag('cdna');
    $writer->dataElement('sequence',$cdna_sequence);
    $writer->endTag('cdna');
  }

  $writer->endTag('transcript');
}

$writer->endTag('fixed_annotation');

$writer->startTag('updatable_annotation');
$writer->dataElement('modification_data',$modification_date);

print_source_tag ($source_name,$url,$contact_name,$contact_address,$email);

if ($other_exon_naming) {
    $writer->startTag('other_exon_naming');
    $writer->startTag('source', 'description' => $exon_url);
    $writer->emptyTag('exon', 'lrg_number' => 1, 'other_name' => 2);
    $writer->endTag('source');
    $writer->endTag('other_exon_naming');
}

if ($mapping_assembly) {
    $writer->startTag('mapping','assembly' => 'NCBI36');
    $writer->emptyTag('align', 'chromosome' => 17, 'strand' => -1, 'lrg_start' => 1, 'lrg_end' => 229, 'start' => 12345, 'end' => 678910);
    $writer->endTag('mapping');
}

if ($amino_acid_mapping) {
    $writer->startTag('amino_acid_mapping');
    $writer->startTag('source', 'description' => $amino_acid_url);
    $writer->dataElement('align', 'lrg_start=1 lrg_end=100 start=1 end=100');
    $writer->endTag('source');
    $writer->endTag('amino_acid_mapping');
}

$writer->startTag('feature');
$writer->endTag('feature');

$writer->startTag('variation');
$writer->endTag('variation');

$writer->endTag('updatable_annotation');

sub print_source_tag {
  my ($source_name,$url,$contact_name,$contact_address,$email) = @_;

      $writer->startTag('source');
      $writer->dataElement('name',$source_name);
      $writer->dataElement('url',$url);
      $writer->startTag('contact');
      $writer->dataElement('name',$contact_name);
      $writer->dataElement('address',$contact_address);
      $writer->dataElement('email',$email);
      $writer->endTag('contact');
      $writer->endTag('source');
}

sub print_tag {

  my ($tag,$text,$attributes) = @_;

  $writer->startTag($tag);
  $writer->characters($text);
  $writer->endTag($tag);

}

