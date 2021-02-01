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


=head1 ImportRGD

This module imports RGD (Rat Genome Database) data.

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportRGD;

use warnings;
use strict;

use File::Path qw(make_path);
use LWP::Simple;

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

my %source_info;

sub fetch_input {
  my $self = shift;

  my $pipeline_dir = $self->required_param('pipeline_dir');
  my $species      = $self->required_param('species');
  my $threshold    = $self->param('threshold_qtl');
  my $run_type     = $self->required_param('run_type');

  $self->debug($self->param('debug_mode'));

  my $rgd_ftp_url_qtl = 'ftp://ftp.rgd.mcw.edu/pub/data_release/';
  my $rgd_ftp_url_gene= 'ftp://ftp.rgd.mcw.edu/pub/data_release/annotated_rgd_objects_by_ontology/with_terms/';

  %source_info = (source_description => 'QTLs from the Rat Genome Database (RGD)',
                  source_url => 'https://rgd.mcw.edu/',
                  object_type => 'QTL', #default type, import code will switch to Gene for the Gene-type phenotypes
                  #source_version  will be set based on the date in the fetched QTL input file  (year/month/day-> yyyymmdd)
                  source_mapped_attrib_type => 'Rat Genome Database', #for ontology mapping (attr_type_id 509) entry in phenotype_ontology_accession (attr_id 588)
                  source_status => 'germline',
                  threshold => $threshold,

                  source_name => 'RGD',        #source name in the variation db
                  source_name_short => 'RGD',  #source identifier in the pipeline
                  data_types => 'phenotype_feature',
                  );
  #NOTE: smart map between species and RGD files Currently only RAT is imported,
  #could be extended for the rest of species: human, mouse.
  my %rgd_names_qtl = (rattus_norvegicus =>'RAT',
                   rat => 'RAT',
                   );
  my %rgd_names_gene = (rattus_norvegicus =>'rattus',
                   rat => 'rattus',
                   );

  my $workdir = $pipeline_dir."/".$source_info{source_name}."/".$species;
  make_path($workdir) or die "Failed to create $workdir $!\n";
  $self->workdir($workdir);

  open(my $logFH, ">", $workdir."/".'log_import_out_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open(my $errFH, ">", $workdir."/".'log_import_err_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open(my $pipelogFH, ">", $workdir."/".'log_import_debug_pipe_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  $self->logFH($logFH);
  $self->errFH($errFH);
  $self->pipelogFH($pipelogFH);

  #get input files RGD eQTL, GENE:
  my $file_qtl = 'QTLS_'.$rgd_names_qtl{$species}.'.txt';
  getstore($rgd_ftp_url_qtl.$file_qtl, $workdir."/".$file_qtl) unless -e $workdir."/".$file_qtl;
  print $logFH  "Found files (".$workdir."/".$file_qtl.") and will skip new fetch\n" if -e $workdir."/".$file_qtl;
  my @files_todo =();
  foreach my $rgd_f (qw/_terms_mp _terms_rdo _terms_vt/){
    my $rgd_file = $rgd_names_gene{$species}.$rgd_f;
    getstore($rgd_ftp_url_gene.$rgd_file, $workdir."/".$rgd_file) unless -e $workdir."/".$rgd_file;
    push @files_todo, $rgd_file;
  }
  $self->param('qtl_file', $file_qtl);
  $self->param('gene_file', [@files_todo]);

  #fetch coreDB assembly:
  my $gc =  $self->core_db_adaptor->get_adaptor('GenomeContainer');
  if ($species eq 'rattus_norvegicus'){
    my @assemblyV = split('_',$gc->get_version); #Rnor_6.0
    $self->param('species_assembly', $assemblyV[1]);  #'6.0';
    print $logFH 'Found core species_assembly:'. $self->param('species_assembly'). "\n" if ($self->debug);
  }

}

sub run {
  my $self = shift;

  # parameters from fetch input
  my $species_assembly = $self->required_param('species_assembly');

  #Process QTLs file
  my $rgd_file = $self->required_param('qtl_file');   #GO through files and parse them in the correct format

  # dump and clean pre-existing phenotype features
  $self->dump_phenotypes($source_info{source_name}, 1);

  # get seq_region_ids
  my $seq_region_ids = $self->get_seq_region_ids();

  # parse phenotypes
  my ($results, $version) = $self->parse_input_file_qtl($seq_region_ids, $rgd_file, $species_assembly);
  $self->print_pipelogFH("Got ".(scalar @{$results->{'phenotypes'}})." phenotypes \n") if ($self->debug);

  # save phenotypes:
  $source_info{source} = 'rgd_qtl';
  $source_info{source_version} = $version;

  $self->save_phenotypes(\%source_info, $results);

  #Process GENEs files
  foreach my $rgd_gene_file (@{$self->required_param('gene_file')}){
    my ($results, $version) = $self->parse_input_file_gene($seq_region_ids, $rgd_gene_file);
    $self->print_pipelogFH("Got ".(scalar @{$results->{'phenotypes'}})." phenotypes from file $rgd_gene_file\n") if ($self->debug);
    $source_info{source} = 'rgd_gene';
    $source_info{object_type} = 'Gene'; # By default it is set to 'QTL

    $self->save_phenotypes(\%source_info, $results);
  }

  my %param_source = (source_name => 'RGD',
                      type => ['QTL', 'Gene']);
  $self->param('output_ids', { source => \%param_source, 
                               species => $self->required_param('species'),
                               run_type => $self->required_param('run_type'),
                             });
}

sub write_output {
  my $self = shift;

  $self->print_pipelogFH("Passing $source_info{source_name_short} import (".$self->required_param('species').") for checks (check_phenotypes)\n") if ($self->debug);
  close($self->logFH) if defined $self->logFH ;
  close($self->errFH) if defined $self->errFH ;
  close($self->pipelogFH) if defined $self->pipelogFH ;

  $self->dataflow_output_id($self->param('output_ids'), 2);
}


=head2 parse_input_file_gene

  Arg [1]    : arrayref $seq_region_ids
              The array of string seq_region_id(s)
  Arg [2]    : string $infile
               The input file name.
  Example    : ($results,$version) = $obj->parse_input_file_gene($seq_region_ids, $infile)
  Description: Parse phenotypes from RGD Gene input file, uses gene symbols lookup in core
  Returntype : hashref with results (key 'phenotypes') and string date_version
  Exceptions : none

=cut

sub parse_input_file_gene {
  my ($self, $seq_region_ids, $infile) = @_;

  my $errFH1;
  open($errFH1, ">", $self->workdir."/".'log_import_err_'.$infile) ;

  my @phenotypes;

  my $gene_adaptor = $self->core_db_adaptor->get_GeneAdaptor;
  die("ERROR: Could not get gene adaptor\n") unless defined($gene_adaptor);

  # Open the input file for reading
  if($infile =~ /gz$/) {
    open(IN, "zcat ".$self->workdir."/$infile |") || die ("Could not open $infile for reading\n");
  }
  else {
    open(IN,'<',$self->workdir."/".$infile) || die ("Could not open $infile for reading\n");
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
      print $errFH1 "$infile source version: ", $version, "\n" if ($self->debug);
    }

    next if /^(\#|\s)/ || !$_;

    my @data_line = split /\t/, $_;

    # header
    if(/^RGD_ID/) {
      $headers{$data_line[$_]} = $_ for 0..$#data_line;
    }
    else {
      die ("ERROR: Could not find header data\n") unless %headers;

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
          print $errFH1 "Symbol '$symbol' not found in the Core database\n";
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
  close ($errFH1);

  my %result = ('phenotypes' => \@phenotypes);
  return (\%result, $version);
}


=head2 parse_input_file_qtl

  Arg [1]    : arrayref $seq_region_ids
               The array of string seq_region_id(s)
  Arg [2]    : string $infile
               The input file name.
  Arg [3]    : string $assembly
               Specific assembly number for fetching the eQTL coordinates.
  Example    : ($results,$version) = $obj->parse_input_file_gene($seq_region_ids, $infile)
  Description: Parse phenotypes from RGD eQTLs input file
  Returntype : hashref with results (key 'phenotypes') and string date_version
  Exceptions : none

=cut

sub parse_input_file_qtl {
  my ($self, $seq_region_ids, $infile, $assembly)  = @_ ;

  my $errFH1;
  open($errFH1, ">", $self->workdir."/".'log_import_err_'.$infile) ;

  my @phenotypes;

  # Open the input file for reading
  if($infile =~ /gz$/) {
    open(IN, "zcat ".$self->workdir."/$infile |") || die ("Could not open $infile for reading\n");
  }
  else {
    open(IN,'<',$self->workdir."/".$infile) || die ("Could not open $infile for reading\n");
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
      print $errFH1 "source version: ", $version, "\n" if ($self->debug);
    }
    next if (/^(\#|\s)/ || !$_) ;

    my @data = split /\t/, $_;
  

    if(/^QTL_RGD_ID/) {     # header
      $headers{$data[$_]} = $_ for 0..$#data;
    } else {
      die ("ERROR: Could not find header data\n") unless %headers;

      my %data;
      $data{$_} = $data[$headers{$_}] for keys %headers;

      # check chromosome
      my $chr = $data{$assembly.'_MAP_POS_CHR'};

      if(!defined($chr) || !$chr) {
        print $errFH1 "WARNING: Could not get coordinates for assembly $assembly on line $line_num\n";
        next;
      }

      $chr =~ s/^chr(om)?\.?//i;

      if(!defined($seq_region_ids->{$chr})) {
        print $errFH1 "WARNING: Could not find seq_region_id for chromosome name $chr on line $line_num\n";
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
  close ($errFH1);

  my %result = ('phenotypes' => \@phenotypes);
  return (\%result, $version);
}

1;
