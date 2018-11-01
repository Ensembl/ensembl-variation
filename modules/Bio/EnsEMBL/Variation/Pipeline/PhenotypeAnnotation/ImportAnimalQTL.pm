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


package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportAnimalQTL;

use strict;
use warnings;

use File::Copy;
use File::Path qw(make_path remove_tree);
use PerlIO::gzip;
use Data::Dumper;

use Bio::EnsEMBL::Registry;

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

my %source_info;
my $workdir;
my $debug;
my $inputFile;

my $core_dba;
my $variation_dba;
my $ontology_db_adaptor;

my $pubmed_prefix = 'PMID:';


#TODO: figure out how to do this automatic in future. currently human download is needed, also make species assembly automatic based on variation assembly
#my $animalqtl_url = 'https://www.animalgenome.org/cgi-bin/QTLdb/index';
#my %animalQTL_species_url = ( 
#    gallus_gallus => 'https://www.animalgenome.org/cgi-bin/QTLdb/GG/download?file=gbpGG_5.0', #Gallus gallus
#    pig => 'https://www.animalgenome.org/cgi-bin/QTLdb/SS/download?file=gbpSS_11.1', #Sus scrofa 
#    sheep => 'https://www.animalgenome.org/cgi-bin/QTLdb/OA/download?file=gbpOAR_3.1',  # Ovis aries
#    cow => 'https://www.animalgenome.org/cgi-bin/QTLdb/BT/download?file=gbpUMD_3.1', #Bos taurus
#    horse => 'https://www.animalgenome.org/cgi-bin/QTLdb/EC/download?file=gbpEC_2.0', #Equus caballus  
#);

sub fetch_input {
    #create output folder structure and fetches input files 
    my $self = shift;

    my $pipeline_dir = $self->required_param('pipeline_dir');
    my $species      = $self->required_param('species');
    my $animalqtl_inputDir = $self->required_param('animalqtl_input_dir');
    my $animalqtl_version = $self->required_param('animalqtl_version');

    $debug = $self->param('debug_mode');

    $core_dba    = $self->get_species_adaptor('core');
    $variation_dba  = $self->get_species_adaptor('variation'); #TODO: why did the init -> Base class -> this not work?
    $ontology_db_adaptor = $self->get_adaptor('multi', 'ontology');

    # import specific constants
    %source_info = (source_description => 'The Animal Quantitative Trait Loci (QTL) database (Animal QTLdb) is designed to house all publicly available QTL and association data on livestock animal species',
                    source_url => 'http://www.animalgenome.org/cgi-bin/QTLdb/index',
                    object_type => 'QTL',
                    source_name => 'Animal_QTLdb', #TODO: is this used anywhere?
                    source_version => $self->required_param('animalqtl_version'),
                    source_status => 'germline',
                    source => 'Animal_QTLdb', #TODO: is this really needed?
                    threshold => $self->required_param('threshold_qtl'),

                #    source_mapped_attrib_type => 'Rat Genome Database', #for ontology mapping (attr_type_id 509) entry in phenotype_ontology_accession (attr_id 588)
                #    set => undef
                    ); 

    #create workdir folder
    $workdir = $pipeline_dir."/".$source_info{source_name}."/".$species;
    make_path($workdir);
    local (*STDOUT, *STDERR);
    open STDOUT, ">", $workdir."/".'log_import_AnimalQTL_out'; #TODO: what is best error/out log naming convention?
    open STDERR, ">", $workdir."/".'log_import_AnimalQTL_err';

    print "AnimalQTL import expects input folder with gff3 files: example format gallus_gallus.*.gff3  \n" if $debug;
    print "using input folder: $animalqtl_inputDir for species: $species \n" if $debug;


  


    #create folder + copy user input folder data here
    print "Fetching data for: ", $species, "\n" if $debug;
    opendir(INDIR, $animalqtl_inputDir);     #TODO: move this read dir into initPipe -> each job will have species and specific file
    my @files = readdir(INDIR);
    closedir(INDIR);
    
    my $ok = 0; 
    foreach my $file (@files){
      if ($file =~/^$species.*gff3$/ || $file =~/^$species.*gff3.gz$/){
        $inputFile = $file; 
        $ok = 1;
      }
    }
    die "Animal_QTLdb file not found for $species in inputDir ($animalqtl_inputDir)! " unless $ok;
    
    print "Found inputDir file: $inputFile ", "\n" if $debug;
    copy($animalqtl_inputDir."/".$inputFile, $workdir."/".$inputFile) unless -e $workdir."/".$inputFile;
    print "Found file (".$workdir."/".$inputFile."), will skip new fetch\n" if -e $workdir."/".$inputFile;
    
#getVErsion
    # check in dir 
    #todo: for each of the animalQTLdb:
    #- create a local workdir folder
    #- download in that folder the correct animalQTLdb file 
    #- for chicken https://www.animalgenome.org/cgi-bin/QTLdb/GG/download?file=gbpGG_5.0

    #fetch coreDB assembly: TODO: is there a nicer/ most reliable way of doing this?
}

sub run {
    my $self = shift;
    
    # get seq_region_ids
    my $seq_region_ids = $self->get_seq_region_ids($variation_dba);
    
    local (*STDOUT, *STDERR);
    open STDOUT, ">", $workdir."/".'log_import_out_'.$inputFile; #TODO: what is best error/out log naming convention?
    open STDERR, ">", $workdir."/".'log_import_err_'.$inputFile;
    
    # get phenotype data
    my $results = parse_animal_qtl($seq_region_ids, $workdir."/".$inputFile);
    print "Got ".(scalar @{$results->{'phenotypes'}})." phenotypes \n" if $debug ;

    # save phenotypes
    $self->save_phenotypes(\%source_info, $results, $core_dba, $variation_dba);
    
    my %param_source = (source_name => $source_info{source_name},
                        type => ['QTL']);
    $self->param('output_ids', { source => \%param_source, 
                                 species => $self->required_param('species')
                               });
}

sub write_output {
    my $self = shift;
    $self->dataflow_output_id($self->param('output_ids'), 1);
    print "Passing AnimalQTL import ($self->required_param('species')) for checks\n" if $self->param('debug_mode');
}

# AnimalQTL specific phenotype parsing method for gff3 files
sub parse_animal_qtl {
  my $seq_region_ids = shift;
  my $infile = shift;
 
  my $ontology_term_adaptor = $ontology_db_adaptor->get_OntologyTermAdaptor; 
  
  my @phenotypes;
  
  # Open the input file for reading
  if($infile =~ /gz$/) {
    open IN, "zcat $infile |" or die ("Could not open $infile for reading");
  }
  else {
    open(IN,'<',$infile) or die ("Could not open $infile for reading");
  }
  
  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;
    
    next if /^(\#|\s)/ || !$_;
    
    my @data = split /\t/, $_;
    
    # fix chr
    $data[0] =~ s/^chr(om)?\.?//i;
    
    if(!defined($seq_region_ids->{$data[0]})) {
      print STDERR "WARNING: Could not find seq_region_id for chromosome name $data[0]\n";
      next;
    }

    # parse "extra" GFF fields
    my $extra = {};
    
    foreach my $chunk(split /\;/, $data[8]) {
      my ($key, $value) = split /\=/, $chunk;
      next unless defined($key) && defined($value);
      $value =~ s/\"|\'//g;
      $extra->{$key} = $value;
    }
    
    if ($extra->{'Map_Type'} eq 'Linkage' || $data[3] eq '' || $data[3] == 0) {
      print STDERR "WARNING: Could not find a precise location for the QTL ".$extra->{QTL_ID}."\n";
      next;
    }

    if ($data[4] !~ /^\d+$/) {
      print STDERR "WARNING: Could not find a numeric seq_region_end for the QTL ".$extra->{QTL_ID}."\n";
      next;
    }

    # create phenotype hash

    my $phenotype = {
      'id' => $extra->{QTL_ID},
      'description' => $extra->{Name},
      'seq_region_id' => $seq_region_ids->{$data[0]},
      'seq_region_start' => $data[3] || 1,
      'seq_region_end' => $data[4],
      'seq_region_strand' => 1
    };
    
    my @accessions = ();
    # CM ontology
    my $CMO_name = $extra->{CMO_name};
    if ($CMO_name) {
      my $terms = $ontology_term_adaptor->fetch_all_by_name($CMO_name, 'CMO'); 
      foreach my $term (@$terms) {
        push @accessions, $term->accession;
      }
    }
    
    # VT ontology
    my $VTO_name = $extra->{VTO_name};
    if ($VTO_name) {
      my $terms = $ontology_term_adaptor->fetch_all_by_name($VTO_name, 'VT'); 
      foreach my $term (@$terms) {
        push @accessions, $term->accession;
      }
    }
    # Ontology
    if (scalar(@accessions) > 0) {
      $phenotype->{accessions} = \@accessions;
      $phenotype->{ontology_mapping_type} = 'is';
    }
    
    if ($phenotype->{'seq_region_start'} > $phenotype->{'seq_region_end'}) {
      my $tmp_end = $phenotype->{'seq_region_end'};
      $phenotype->{'seq_region_start'} = $phenotype->{'seq_region_end'};
      $phenotype->{'seq_region_end'}   = $tmp_end;
    }
    
    # add additional fields if found
    $phenotype->{'study'} = $pubmed_prefix.$extra->{'PUBMED_ID'} if defined($extra->{'PUBMED_ID'} && $extra->{'PUBMED_ID'} =~ /^\d+$/);
    $phenotype->{'p_value'} = $extra->{'P-value'} if defined($extra->{'P-value'});
    $phenotype->{'f_stat'} = $extra->{'F-stat'} if defined($extra->{'F-stat'});
    $phenotype->{'lod_score'} = $extra->{'LOD-score'} if defined($extra->{'LOD-score'});
    $phenotype->{'variance'} = $extra->{'Variance'} if defined($extra->{'Variance'});
    $phenotype->{'associated_gene'} = $extra->{'Candidate_Gene_Symble'} if defined($extra->{'Candidate_Gene_Symble'});
    
    push @phenotypes, $phenotype;
  }
  
  my %result = ('phenotypes' => \@phenotypes);
  return \%result;
}

1;

