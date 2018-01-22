#!/usr/bin/env perl

=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

=head1 NAME

veplotyper.pl

by Will McLaren (wm2@ebi.ac.uk)
=cut

use strict;
use Getopt::Long;
use FileHandle;
use JSON;
use HTTP::Tiny;
use Digest::MD5 qw(md5_hex);
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Individual;
use Bio::EnsEMBL::Variation::Population;
use Bio::EnsEMBL::Variation::IndividualGenotypeFeature;
use Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer;
use Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor;
use Bio::EnsEMBL::Variation::Utils::VEP qw(
  prefetch_transcript_data
  parse_line
  load_dumped_transcript_cache
);

# set output autoflush for progress bars
$| = 1;

# configure from command line opts
my $config = configure(scalar @ARGV);

# run the main sub routine
main($config);

sub configure {
  my $args = shift;
  
  my $config = {};
  
  GetOptions(
    $config,
    'help|h',                  # displays help message
    'quiet|q',                 # no status output
    
    'transcript|t=s',          # transcript ID to fetch
    
    'database',                # indicate we want to use the database
    'cache',                   # indicate we want to use the VEP cache
    
    'vcf_file|v=s',                 # vcf file
    'collection|n=s',          # JSON config file for a VCFCollection
    'output_file|o=s',         # output file
    'force_overwrite',         # force overwrite output
    
    'rest|r=s',                # URL of REST server
    
    'individual=s',            # individual name to limit analysis to
    'panel|p=s',               # panel file containing individual population
    'individual_prefix|x=s',   # prefix to add to individual names if looking up in database
    
    'host=s',                  # DB options
    'user=s',
    'pass=s',
    'port=i',
    'version=i',
    'registry=s',
    'species=s',
    
    'cache_version|c=s',       # cache version if looking up transcripts from VEP cache
  ) or die "ERROR: Failed to parse command-line flags\n";
  
  # print usage message if requested or no args supplied
  if(defined($config->{help}) || !$args) {
    &usage;
    exit(0);
  }
  
  
  debug($config, 'Starting');
  
  # defaults
  $config->{rest}              ||= 'http://rest.ensembl.org';
  $config->{vcf_file}          ||= 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr###CHR###.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz';
  $config->{cache_region_size} ||= 1000000;
  $config->{chunk_size}        ||= 50000;
  $config->{compress}          ||= 'zcat';
  $config->{terms}             ||= 'SO';
  $config->{cache}             ||= 1;
  $config->{format}            ||= 'vcf';
  $config->{polyphen_analysis} ||= 'humvar';
  $config->{sift}              ||= 1;
  $config->{polyphen}          ||= 1;
  $config->{species}           ||= 'homo_sapiens';
  $config->{individual}        ||= 'all';
  $config->{output_file}       ||= 'veplotyper_out.xml';
  
  $config->{individual} = [split(',', $config->{individual})];
  
  # get adaptors
  $config->{vfa}  = Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor->new_fake($config->{species});
  $config->{tva}  = Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor->new_fake($config->{species});
  
  
  $config->{reg} = 'Bio::EnsEMBL::Registry';
  
  debug($config, 'Testing connections');
  
  # DB connection?
  if(defined($config->{registry})) {
    $config->{reg}->load_all($config->{registry});
  }
  elsif(defined($config->{host})) {
    $config->{reg}->load_registry_from_db(
      -host       => $config->{host},
      -user       => $config->{user},
      -pass       => $config->{password},
      -port       => $config->{port},
      -db_version => $config->{version},
      -species    => $config->{species} =~ /^[a-z]+\_[a-z]+/i ? $config->{species} : undef,
      #-no_cache   => 1,
    );
    
    $config->{db} = $config->{reg}->get_adaptor($config->{species}, 'variation', 'variation')->db;
    $config->{pfpma} = $config->{reg}->get_adaptor($config->{species}, 'variation', 'proteinfunctionpredictionmatrix');
  }
  
  # cache
  $config->{dir} ||= $config->{dir_cache} || join '/', ($ENV{'HOME'}, '.vep');  
  $config->{dir} .= '/'.(
    join '/', (
      $config->{species},
      $config->{cache_version} || $config->{reg}->software_version
    )
  );
  
  # die("ERROR: cache directory ".$config->{dir}." not found\n") unless -d $config->{dir};
  
  # check tabix
  die "ERROR: tabix does not seem to be in your path\n" unless `which tabix 2>&1` =~ /tabix$/;
  
  debug($config, "Checking VCF file(s)");
  
  # check VCF
  my $vcf_file = $config->{vcf_file};
  
  if($vcf_file) {
    # check for ###CHR###
    $vcf_file =~ s/###CHR###/1/;
  
    # remote
    if($vcf_file =~ /tp\:\/\//) {
    
      my $remote_test = `tabix -f $vcf_file 1:1-1 2>&1`;
      if($remote_test =~ /fail/) {
        die "$remote_test\nERROR: Could not find file or index file for remote annotation file $vcf_file\n";
      }
    }
    # local
    else {
      die "ERROR: Custom annotation file $vcf_file not found\n" unless -e $vcf_file;
      die "ERROR: Tabix index file $vcf_file\.tbi not found - perhaps you need to create it first?\n" unless -e $vcf_file.'.tbi';
    }
  }
  
  # output file
  die("ERROR: output file ".$config->{output_file}." already exists") if -e $config->{output_file} && !defined($config->{force_overwrite});
  my $oh = FileHandle->new(">".$config->{output_file});
  $config->{out_handle} = $oh;
  
  return $config;
}

sub main {
  my $config = shift;
  
  # fetch transcript
  debug($config, 'Fetching transcript data');
  my $tr = fetch_transcript($config, $config->{transcript});
  
  # now fetch genotypes that overlap
  debug($config, 'Fetching variant data');
  my $gts = fetch_genotypes($config, $tr->{chr}, $tr->{start}, $tr->{end});
  
  $DB::single = 1;
  
  debug($config, 'Mutating sequences');
  my $c = Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer->new($tr, $gts, $config->{db});
  
  $DB::single = 1;
  
  foreach my $p(sort {$b->count <=> $a->count} @{$c->get_all_ProteinHaplotypes}) {
    
    my $f = $p->get_all_population_frequencies;
    
    printf(
      "AFR:%.2f\tAMR:%.2f\tASN:%.2f\tEUR:%.2f\t%.2f\%\tSIFT:%.3f\t%s\n",
      100 * ($f->{AFR} || $f->{'1000GENOMES:phase_1_AFR'}),
      100 * ($f->{AMR} || $f->{'1000GENOMES:phase_1_AMR'}),
      100 * ($f->{ASN} || $f->{'1000GENOMES:phase_1_ASN'}),
      100 * ($f->{EUR} || $f->{'1000GENOMES:phase_1_EUR'}),
      100 * $p->frequency,
      $p->mean_sift_score,
      $p->name
    );
  }
  
  $DB::single = 1;
  
  debug($config, 'Writing output');
  
  # write to JSON
  write_JSON($config, $c);
  
  debug($config, 'All done');
  
  1;
}

sub write_JSON {
  my $config = shift;
  
  my $c = shift;
  my $name = $c->transcript->stable_id;
  
  use JSON;
  my $json = JSON->new;
  
  my $oh = $config->{out_handle};
  print $oh $json->allow_blessed->convert_blessed->pretty->encode($c);
  #  {
  #    transcript_id => $name,
  #    populations => $config->{pop_counts},
  #    cds_haplotypes => [sort {$b->{count} <=> $a->{count}} grep {$_->{type} eq 'cds'} values %$alt_seqs],
  #    peptide_haplotypes => [sort {$b->{count} <=> $a->{count}} grep {$_->{type} eq 'pep'} values %$alt_seqs],
  #  }
  #);
}

sub fetch_transcript {
  my $config = shift;
  my $tr_id = shift;
  
  # use API
  if(defined($config->{host})) {
    return fetch_transcript_api($config, $tr_id);
  }
  else {
    return fetch_transcript_cache($config, $tr_id);
  }
}

sub fetch_transcript_api {
  my $config = shift;
  my $tr_id = shift;
  
  $config->{ta} ||= $config->{reg}->get_adaptor($config->{species}, 'core', 'transcript');
  my $tr = $config->{ta}->fetch_by_stable_id($tr_id);
  
  #$DB::single = 1;
  die("ERROR: Could not fetch transcript\n") unless defined($tr);
  
  prefetch_transcript_data($config, $tr);
  
  # cache chromosome name
  $tr->{chr} = $tr->seq_region_name;
  
  return $tr;
}

sub fetch_transcript_cache {
  my $config = shift;
  my $tr_id = shift;
  
  # get transcript coords from REST API
  my $tr_rest = rest_fetch($config, 'lookup/id', $tr_id, 'expand=0');
  
  my $c = $tr_rest->{seq_region_name};
  my $s = int($tr_rest->{start} / $config->{cache_region_size});
  my $e = $s + 1;
  $s  = ($s * $config->{cache_region_size}) + 1;
  $e *= $config->{cache_region_size};
  
  my $q = $config->{quiet};
  $config->{quiet} = 1;
  my $tr_cache = load_dumped_transcript_cache($config, $c, "$s-$e");
  $config->{quiet} = $q;
  
  my ($tr) = grep {$_->stable_id eq $tr_id} @{$tr_cache->{$c}};
  
  # cache chromosome name
  $tr->{chr} = $c;
  
  return $tr;
}

sub fetch_genotypes {
  my $config = shift;
  
  #if(1) {
  #  return fetch_genotypes_db($config, @_);
  #}
  if(defined($config->{collection})) {
    return fetch_genotypes_collection($config, @_);
  }
  else {
    return fetch_genotypes_vcf($config, @_);
  }
}

sub fetch_genotypes_collection {
  my $config = shift;
  
  $Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor::CONFIG_FILE = $config->{collection};
  
  my ($ca, $sa, $slice);
  
  if(defined($config->{host})) {
    $ca = $config->{reg}->get_adaptor($config->{species}, 'variation', 'vcfcollection');
    $sa = $config->{reg}->get_adaptor($config->{species}, 'core', 'slice');
    $slice = $sa->fetch_by_region('chromosome', @_);
  }
  
  else {
    $ca = Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor->new();
    $slice = Bio::EnsEMBL::Slice->new_fast({
      seq_region_name => $_[0],
      start => $_[1],
      end => $_[2],
    });
  }
  
  $_->use_db(0) for @{$ca->fetch_all};
  
  # add populations to individuals
  foreach my $c(@{$ca->fetch_all}) {
    my $hash = $c->_get_Population_Individual_hash;
    my $pops = $c->get_all_Populations();
    my $inds = $c->get_all_Individuals();
    
    my %pop_dbID_map = map {$_->dbID => $_} @$pops;
    my %ind_dbID_map = map {$_->dbID => $_} @$inds;
    
    foreach my $pop_id(keys %$hash) {
      my $pop = $pop_dbID_map{$pop_id};
      
      foreach my $ind_id(keys %{$hash->{$pop_id}}) {  
        my $ind = $ind_dbID_map{$ind_id};
        
        push @{$ind->{populations}}, $pop;
      }
    }
  }
  
  $DB::single = 1;
  
  my @gts = map {@{$_->get_all_IndividualGenotypeFeatures_by_Slice($slice, undef, 1)}} @{$ca->fetch_all};
  my @return;
  
  # filter out ref homozygotes
  foreach my $gt(@gts) {
    my $ref = (split '/', $gt->variation_feature->allele_string)[0];
    my @gt = @{$gt->genotype};
    next if scalar (grep {$_ eq $ref} @gt) == scalar @gt;
    push @return, $gt;
  }
  
  return \@return;
}

sub fetch_genotypes_db {
  my $config = shift;
  
  # get a slice
  my $sa = $config->{reg}->get_adaptor($config->{species}, 'core', 'slice');
  my $slice = $sa->fetch_by_region('chromosome', @_);
  
  # get variation features
  my %vfs = map {$_->{_variation_id} => $_} @{$slice->get_all_VariationFeatures};
  
  # get phased genotypes
  my $gta = $config->{reg}->get_adaptor($config->{species}, 'variation', 'individualgenotypefeature');
  return [
    grep {$_->{variation_feature}}
    map {$_->{variation_feature} = $vfs{$_->{_variation_id}}; $_}
    #grep {$_->phased}
    @{$gta->fetch_all_by_Slice($slice)}
  ];
}
  
sub fetch_genotypes_vcf {
  my $config = shift;
  
  my ($c, $s, $e) = @_;
  
  my $file = $config->{vcf_file};
  $file =~ s/###CHR###/$c/;
  
  get_vcf_individuals($config, $file);
  
  # get populations
  get_populations($config);
  
  open VARS, "tabix -f $file $c:$s-$e |"
    or die "\nERROR: Could not open tabix pipe for $file\n";
  
  my @vars;
  
  VAR: while(<VARS>) {
    chomp;
    push @vars, grep {!defined($_->{hom_ref}) && !$_->isa('Bio::EnsEMBL::Variation::StructuralVariationFeature')} @{parse_line($config, $_)};
  }
  
  close VARS;
  
  my $prefix = $config->{individual_prefix} || '';
  
  # convert vars to IndividualGenotypeFeature objects, keep only phased ones
  @vars = map {
    Bio::EnsEMBL::Variation::IndividualGenotypeFeature->new_fast({
      variation_feature => $_,
      genotype => $_->{genotype},
      individual => Bio::EnsEMBL::Variation::Individual->new_fast({
        name => $_->{individual},
        populations => [map {$config->{pops}->{$_}} keys %{$config->{ind_pops}->{$_->{individual}}}],
      })
    })
  } grep {$_->{phased}} @vars;
  
  return \@vars;
}

sub get_vcf_individuals {
  my $config = shift;
  my $vcf = shift;
  
  if(!defined($config->{ind_cols})) {
    open IN, "tabix -hf $vcf 1:1-1 |";
    while(<IN>) {
      if(/^\#CHROM/) {
        chomp;
        my @split = split /\s+/;
        
        # no individuals
        die("ERROR: No individual data found in VCF\n") if scalar @split <= 9;
        
        my $prefix = $config->{individual_prefix} || '';
        
        # get individual column indices
        my %ind_cols = map {$prefix.$split[$_] => $_} (9..$#split);
        
        # all?
        if(scalar @{$config->{individual}} == 1 && $config->{individual}->[0] =~ /^all$/i) {
          $config->{ind_cols} = \%ind_cols;
        }
        else {
          my %new_ind_cols;
          
          # check we have specified individual(s)
          foreach my $ind(@{$config->{individual}}) {
            die("ERROR: Individual named \"$ind\" not found in VCF\n") unless defined $ind_cols{$ind};
            $new_ind_cols{$ind} = $ind_cols{$ind};
          }
          
          $config->{ind_cols} = \%new_ind_cols;
        }
      }
    }
    close IN;
  }
  
  return $config->{ind_cols};
}

sub get_populations {
  my $config = shift;
  
  if(!defined($config->{ind_pops})) {
    
    my %populations = ();
    
    # get from panel file
    if(defined($config->{panel})) {
      open IN, $config->{panel} or die("ERROR: Could not read from panel file ".$config->{panel}."\n");
      while(<IN>) {
        chomp;
        my @split = split /\s+/, $_;
        my $i = shift @split;
        #next unless defined($config->{ind_cols}->{$i});
        
        # could be more than one
        foreach my $p(@split) {
          $populations{$i}{$p} = 1;
          $config->{pop_counts}->{$p}++;
          $config->{pops}->{$p} ||= Bio::EnsEMBL::Variation::Population->new_fast({name => $p});
        }
      }
      close IN;
    }
    
    # get from database
    elsif(defined($config->{host})) {
      
      # look up individual objects
      my $ia = $config->{reg}->get_adaptor($config->{species}, 'variation', 'individual');
      my $prefix = $config->{individual_prefix} || '';
      my %ind_names_by_dbID = map {$_->dbID => $_->name} @{$ia->fetch_all_by_name_list([keys %{$config->{ind_cols}}])};
      
      my $pa = $config->{reg}->get_adaptor($config->{species}, 'variation', 'population');
      my $hash = $pa->_get_individual_population_hash([keys %ind_names_by_dbID]);
      my @pops = @{$pa->fetch_all_by_dbID_list([keys %$hash])};
      my %pop_names_by_dbID = map {$_->dbID => $_->name} @pops;
      %{$config->{pops}} = map {$_->name => $_} @pops;
      
      # "invert" the hash
      foreach my $p_id(keys %$hash) {
        my $p = $pop_names_by_dbID{$p_id};
        $populations{$_}{$p} = 1 for map {$ind_names_by_dbID{$_}} keys %{$hash->{$p_id}};
      }
      
      # log counts
      %{$config->{pop_counts}} = map {$pop_names_by_dbID{$_} => scalar keys %{$hash->{$_}}} keys %$hash;
    }
    
    $config->{ind_pops} = \%populations;
  }
  
  return $config->{ind_pops};
}

sub rest_fetch {
  my $config = shift;
  my $endpoint = shift;
  my $value = shift;
  my $extra = shift;
  
  my $url =
    $config->{rest}.'/'.
    $endpoint.'/'.
    $value.
    '?content-type=application/json;'.
    $extra;
  
  my $response = HTTP::Tiny->new()->get($url);
  
  die "ERROR: Failed to fetch transcript $config->{transcript} from $url\n" unless $response->{success};
  
  return decode_json($response->{content});
}

sub debug {
  my $config = shift;
  return if defined($config->{quiet});
  my $msg = shift;
  $msg .= "\n" unless $msg =~ /\n$/;
  print STDERR get_time().' - '.$msg;
}

sub get_time() {
  my @time = localtime(time());

  # increment the month (Jan = 0)
  $time[4]++;

  # add leading zeroes as required
  for my $i(0..4) {
    $time[$i] = "0".$time[$i] if $time[$i] < 10;
  }

  # put the components together in a string
  my $time =
    ($time[5] + 1900)."-".
    $time[4]."-".
    $time[3]." ".
    $time[2].":".
    $time[1].":".
    $time[0];

  return $time;
}

sub usage {
  print qq{#---------------#
# veplotyper.pl #
#---------------#

By Will McLaren (wm2\@ebi.ac.uk)

Usage:
perl veplotyper.pl [arguments]
  
--help               -h   Print usage message and exit
};
}
