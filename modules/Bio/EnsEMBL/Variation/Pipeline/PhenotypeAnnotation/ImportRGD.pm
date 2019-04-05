=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportRGD;

use warnings;
use strict;

use File::Path qw(make_path);
use LWP::Simple;
use Data::Dumper; #TODO: remove if not needed
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::Constants qw(RGD NONE);
#use Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation qw($variation_dba);

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

#TODO: should contain everything about RGD, including getting the file(s) and for all RGD species

my %source_info;
my $workdir;
my $core_dba;
my $variation_dba;
my $phenotype_dba;

my $debug;

sub fetch_input {
    my $self = shift;

    my $pipeline_dir = $self->required_param('pipeline_dir');
    my $species      = $self->required_param('species');
    
    $core_dba    = $self->get_species_adaptor('core');
    $variation_dba  = $self->get_species_adaptor('variation'); #TODO: why did the init -> Base class -> this not work?
    $phenotype_dba  = $variation_dba->get_PhenotypeAdaptor; 
    
    $debug        = $self->param('debug_mode');
  
    #TODO:info RAT species assembly used to be an input parameter, now is parsed from core db and used latest

    my $rgd_ftp_url_qtl = 'ftp://ftp.rgd.mcw.edu/pub/data_release/';
    my $rgd_ftp_url_gene= 'ftp://ftp.rgd.mcw.edu/pub/data_release/annotated_rgd_objects_by_ontology/with_terms/';

    %source_info = (source_description => 'QTLs from the Rat Genome Database (RGD)',
                    source_url => 'http://rgd.mcw.edu/',
                    object_type => 'QTL', #default type, import code will switch to Gene for the Gene-type ones
                    source_name => 'RGD',
                    source_mapped_attrib_type => 'Rat Genome Database', #for ontology mapping (attr_type_id 509) entry in phenotype_ontology_accession (attr_id 588)
                    source_status => 'germline',
                    threshold => $self->required_param('threshold_qtl'),
                    );
    #TODO: smart map between species and RGD files Currently only RAT is imported... could remove the rest from here
    my %rgd_names_qtl = (rattus_norvegicus =>'RAT',
                     rat => 'RAT',
                     );
    my %rgd_names_gene = (rattus_norvegicus =>'rattus',
                     rat => 'rattus',
                     );
    
    $workdir = $pipeline_dir."/".$source_info{source_name}."/".$species;
    #get input files RGD eQTL, GENE:
    my $file_qtl = 'QTLS_'.$rgd_names_qtl{$species}.'.txt';
    make_path($workdir);
    getstore($rgd_ftp_url_qtl.$file_qtl, $workdir."/".$file_qtl) unless -e $workdir."/".$file_qtl;
    print "Found files (".$workdir."/".$file_qtl.") and will skip new fetch\n" if -e $workdir."/".$file_qtl;
    my @files_todo =();
    foreach my $rgd_f (qw/_terms_mp _terms_rdo _terms_vt/){
      my $rgd_file = $rgd_names_gene{$species}.$rgd_f;
      getstore($rgd_ftp_url_gene.$rgd_file, $workdir."/".$rgd_file) unless -e $workdir."/".$rgd_file;
      push @files_todo, $rgd_file;
    }
    $self->param('qtl_file', $file_qtl);
    $self->param('gene_file', [@files_todo]);
    #TODO: could/should I(?) add a warning if file already exists? latest: if file there it skips fetch from ftp
    
    #fetch coreDB assembly: TODO: is there a nicer/ most reliable way of doing this?
    my $gc =  $core_dba->get_adaptor('GenomeContainer');
    if ($species eq 'rattus_norvegicus'){
      my @assemblyV = split('_',$gc->get_version); #Rnor_6.0
      $self->param('species_assembly', $assemblyV[1]);  #'6.0';
      print 'Found core species_assembly:'. $self->param('species_assembly'). "\n" if $debug;
    }

}

sub run {
  my $self = shift;

  # parameters from fetch input
  my $species_assembly = $self->required_param('species_assembly');
  
  #PROCESS QTLs file
  my $rgd_file = $self->required_param('qtl_file');   #GO through files and parse them in the correct format
  
  # get seq_region_ids
  my $seq_region_ids = $self->get_seq_region_ids($variation_dba);

  # parse phenotypes
  {
    local (*STDOUT, *STDERR);
    open STDOUT, ">", $workdir."/".'log_import_out_'.$rgd_file; #TODO: what is best error/out log naming convention?
    open STDERR, ">", $workdir."/".'log_import_err_'.$rgd_file;
    my ($results, $version) = parse_rgd_qtl($seq_region_ids, $rgd_file, $species_assembly);
    warn "Got ".(scalar @{$results->{'phenotypes'}})." phenotypes \n" if $debug ;
  
    # save phenotypes:
    $source_info{source} = 'rgd_qtl';
    $source_info{source_version} = $version;

    $self->save_phenotypes(\%source_info, $results, $core_dba, $variation_dba);
  }

  #TODO: the last read in file saves also the version; TODO:Q: should I only use the version from QTL file?
  #PROCESS GENEs files 
  foreach my $rgd_gene_file (@{$self->required_param('gene_file')}){
  #TODO: remove after test  $rgd_gene_file = 'rattus_terms_vt';
    local (*STDOUT, *STDERR);
    open STDOUT, ">", $workdir."/".'log_import_out_'.$rgd_gene_file; #TODO: what is best error/out log naming convention?
    open STDERR, ">", $workdir."/".'log_import_err_'.$rgd_gene_file;
    my ($results, $version) = parse_rgd_gene($seq_region_ids, $rgd_gene_file);
    warn "Got ".(scalar @{$results->{'phenotypes'}})." phenotypes from file $rgd_gene_file\n" if $debug ;
    $source_info{source} = 'rgd_gene';
    $source_info{object_type} = 'Gene'; # By default it is set to 'QTL

    $self->save_phenotypes(\%source_info, $results, $core_dba, $variation_dba);
  }

  my %param_source = (source_name => 'RGD',
                      type => ['QTL', 'Gene']);
  $self->param('output_ids', { source => \%param_source, 
                               species => $self->required_param('species'),
                             });
}

sub write_output {
  my $self = shift;

  $self->dataflow_output_id($self->param('output_ids'), 1);
  print "Passing RGD import for checks\n" if $self->param('debug_mode');
}


# RGD specific phenotype parsing method for GENE files
sub parse_rgd_gene {
  my $seq_region_ids = shift;
  my $infile = shift;
  
  my @phenotypes;

  my $gene_adaptor = $core_dba->get_GeneAdaptor;
  die("ERROR: Could not get gene adaptor") unless defined($gene_adaptor);

  # Open the input file for reading
  if($infile =~ /gz$/) {
    open IN, "zcat $workdir."/".$infile |" or die ("Could not open $infile for reading");
  }
  else {
    open(IN,'<',$workdir."/".$infile) or die ("Could not open $infile for reading");
  }
  
  my %rgd_coords;

  my (%headers, $line_num);

  # Read through the file and parse out the desired fields
  my $version;
  while (<IN>) {
    chomp;
    $line_num++;
    
    if (/^\#\s+GENERATED-ON:\s+(.*)/){
      $version = $1;
      $version =~ s/\///g;
      print "$infile source version: ", $version, "\n" if $debug;
    }

    next if /^(\#|\s)/ || !$_;

    my @data_line = split /\t/, $_;

    # header
    if(/^RGD_ID/) {
      $headers{$data_line[$_]} = $_ for 0..$#data_line;
    }
    else {
      die "ERROR: Couldn't find header data\n" unless %headers;

      my %data;
      $data{$_} = $data_line[$headers{$_}] for keys %headers;

      next unless ($data{'OBJECT_TYPE'} =~ /gene/i);

      my $symbol = $data{'OBJECT_SYMBOL'};

      if (!$rgd_coords{$symbol}) {

        my $gene = $gene_adaptor->fetch_by_display_label($symbol);

        if (!$gene) {
          my $genes = $gene_adaptor->fetch_all_by_external_name($symbol,'RGD');
          $gene = $genes->[0] if (scalar(@$genes) > 0);
        }

        if (!$gene) {
          print STDERR "Symbol '$symbol' not found in the Core database\n";
          next;
        }

        $rgd_coords{$symbol} = { 'gene'   => $gene->stable_id,
                                 'chr'    => $gene->slice->seq_region_name,
                                 'start'  => $gene->start,
                                 'end'    => $gene->end,
                                 'strand' => $gene->strand
                               };
      }

      my $phenotype = {
        id => $rgd_coords{$symbol}{'gene'},
        description => $data{'TERM_NAME'},
        seq_region_id => $seq_region_ids->{$rgd_coords{$symbol}{'chr'}},
        seq_region_start => $rgd_coords{$symbol}{'start'},
        seq_region_end => $rgd_coords{$symbol}{'end'},
        seq_region_strand => $rgd_coords{$symbol}{'strand'},
      };

      if ($data{'REFERENCES'} =~ /^(RGD\:\d+)\|PMID\:(\d+)/) {
        $phenotype->{external_id} = $1;
        $phenotype->{pubmed_id} = $2;
      }

      if ($data{'TERM_ACC_ID'} =~ /^MP\:\d+/) {
        $phenotype->{accessions} = [ $data{'TERM_ACC_ID'} ];
      }

      push @phenotypes, $phenotype;
    }
  }
  close IN;

  my %result = ('phenotypes' => \@phenotypes);
  return (\%result, $version);
}


# RGD specific phenotype parsing method for eQTLs
sub parse_rgd_qtl {
  my $seq_region_ids = shift;
  my $infile = shift;
  my $assembly = shift;

  my @phenotypes;

  # Open the input file for reading
  if($infile =~ /gz$/) {
    open IN, "zcat $workdir."/".$infile |" or die ("Could not open $infile for reading");
  }
  else {
    open(IN,'<',$workdir."/".$infile) or die ("Could not open $infile for reading");
  }

  my (%headers, $line_num);

  # Read through the file and parse out the desired fields
  my $version;
  while (<IN>) {
    chomp;
    $line_num++;
    
    if (/^\#\s+GENERATED-ON:\s+(.*)/){
      $version = $1;
      $version =~ s/\///g;
      print "source version: ", $version, "\n" if $debug;
    }
    next if (/^(\#|\s)/ || !$_) ;

    my @data = split /\t/, $_;
  

    if(/^QTL_RGD_ID/) {     # header
      $headers{$data[$_]} = $_ for 0..$#data;
    } else {
      die "ERROR: Couldn't find header data\n" unless %headers;

      my %data;
      $data{$_} = $data[$headers{$_}] for keys %headers;

      # check chromosome
      my $chr = $data{$assembly.'_MAP_POS_CHR'};

      if(!defined($chr) || !$chr) {
        print STDERR "WARNING: Could not get coordinates for assembly $assembly on line $line_num\n";
        next;
      }

      $chr =~ s/^chr(om)?\.?//i;

      if(!defined($seq_region_ids->{$chr})) {
        print STDERR "WARNING: Could not find seq_region_id for chromosome name $chr on line $line_num\n";
        next;
      }

      my $description = $data{TRAIT_NAME};
      my @description_acc;
      if ($description =~ /(.*)\s\((VT:.*)\)/) {
        $description = $1;
        push(@description_acc, $2);
      }

      my $pubmed_id = $data{CURATED_REF_PUBMED_ID};
      $pubmed_id =~ s/;/,/g;

      my $phenotype = {
        id => $data{QTL_SYMBOL},
        description => $description,
        accessions => \@description_acc,
        ontology_mapping_type => 'is',
        pubmed_id => $pubmed_id,
        external_id => $data{QTL_RGD_ID},
        seq_region_id => $seq_region_ids->{$chr},
        seq_region_start => $data{$assembly.'_MAP_POS_START'},
        seq_region_end => $data{$assembly.'_MAP_POS_STOP'},
        seq_region_strand => 1,
        lod_score => $data{LOD},
        p_value => $data{P_VALUE},
        variance => $data{VARIANCE},
        associated_gene => $data{CANDIDATE_GENE_SYMBOLS}
      };

      push @phenotypes, $phenotype;
    }
  }
  close IN;
  
  my %result = ('phenotypes' => \@phenotypes);
  return (\%result, $version);
}

1;