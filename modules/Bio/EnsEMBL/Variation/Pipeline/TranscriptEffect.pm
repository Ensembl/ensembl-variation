=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Variation::Pipeline::TranscriptEffect;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::TranscriptVariation;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(MAX_DISTANCE_FROM_TRANSCRIPT overlap);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);

use ImportUtils qw(load);
use FileHandle;
use Fcntl qw(:flock SEEK_END);
use Digest::MD5 qw(md5_hex);

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

my $DEBUG   = 0;

sub run {
  my $self = shift;

  my $gene_id = $self->required_param('gene_stable_id'); 

  my $disambiguate_sn_alleles = 
    $self->param('disambiguate_single_nucleotide_alleles');

  my $mtmp = $self->param('mtmp_table');
  
  my $variations_to_include;
  
  # if (my $vars = $self->param('variations_to_include')) {
  #   # turn the list of variation names into a hash to speed up checking
  #   $variations_to_include = { map { $_ => 1 } @$vars };
  # }

  # set max distance
  my $opt_max_distance = $self->param('max_distance');
  my $max_distance;

  if(defined($opt_max_distance)) {
    $max_distance = $opt_max_distance;
    $Bio::EnsEMBL::Variation::Utils::VariationEffect::UPSTREAM_DISTANCE = $opt_max_distance;
    $Bio::EnsEMBL::Variation::Utils::VariationEffect::DOWNSTREAM_DISTANCE = $opt_max_distance;
  }
  else {
    $max_distance = MAX_DISTANCE_FROM_TRANSCRIPT;
  }

  # clear the registry here
  # this hopefully prevents any sequence caching issues
  # overhanging from previous jobs executed in the same hive process
  Bio::EnsEMBL::Registry->clear();

  my $core_dba = $self->get_species_adaptor('core');
  my $var_dba = $self->get_species_adaptor('variation');
  
  my $ga = $core_dba->get_GeneAdaptor;
  my $sa = $core_dba->get_SliceAdaptor;
  
  my $tva = $var_dba->get_TranscriptVariationAdaptor;

  # we need to include failed variations
  $tva->db->include_failed_variations(1);

  if((my $fasta = $self->param('fasta')) && !$Bio::EnsEMBL::Slice::_fasta_redefined && !$Bio::EnsEMBL::Slice::fasta_db) {

    # we need to find the assembly version to tell it about PARs
    my ($highest_cs) = @{$core_dba->get_CoordSystemAdaptor->fetch_all()};
    my $assembly = $highest_cs->version();

    setup_fasta(-FASTA => $fasta, -ASSEMBLY => $assembly);
  }
  else {
    # set seq cache size higher
    # this prevents repeated DB lookups for sequence for HGVS e.g. in long introns
    # to do this we fetch a sequence adaptor
    my $seq_ad = $core_dba->get_SequenceAdaptor;

    # then reset its cache
    # the first param passed to this method is the "chunk power"
    # essentially the length in bp of cached sequence will be 2 ^ chunk_power
    $seq_ad->_init_seq_instance($seq_ad->chunk_power + 2);
  }

  print STDERR "Fetching gene $gene_id\n" if $DEBUG;

  my $gene = $ga->fetch_by_stable_id($gene_id) 
    or die "failed to fetch gene for stable id: $gene_id";

  my $gene_name;
  $gene_name = $gene->display_xref->display_id if defined $gene->display_xref();

  my $slice = $sa->fetch_by_gene_stable_id(
    $gene->stable_id, 
    $max_distance
  ) or die "failed to get slice around gene: ".$gene->stable_id;
  
  # call seq here to help cache
  $slice->seq();

  $gene = $gene->transfer($slice);

  print STDERR "Getting variation features\n" if $DEBUG;

  my $vfa = $var_dba->get_VariationFeatureAdaptor();

  my @vfs = (
    @{ $vfa->fetch_all_by_Slice_SO_terms($slice) },
    @{ $vfa->fetch_all_somatic_by_Slice_SO_terms($slice) }
  );

  my $table_files_dir = $self->get_table_files_prefix($gene_id);

  print STDERR "STUB $table_files_dir\n";

  # write VFs to table file for web search indexes
  if($gene_name) {
    my $fh = FileHandle->new();
    $fh->open(">".$table_files_dir."_variation_genename.txt") or die "Cannot open dump file ".$table_files_dir."_variation_genename.txt: $!";
    print $fh "$_\t$gene_name\n" for map {$_->get_Variation_dbID()} @vfs;
    $fh->close();
  }

  # get a fh for the hgvs file too
  my $hgvs_fh = FileHandle->new();
  $hgvs_fh->open(">".$table_files_dir."_variation_hgvs.txt") or die "Cannot open dump file ".$table_files_dir."_variation_hgvs.txt: $!";


  my @write_data;
  my $t0;

  my $hgvs_by_var = {};
  my $var_count = 0;

  # initialise a hash of files
  my $files = {
    transcript_variation      => { 'cols' => [$tva->_write_columns],      },
    MTMP_transcript_variation => { 'cols' => [$tva->_mtmp_write_columns], },
    # variation_genename        => { 'cols' => [qw(variation_id gene_name)] },
    # variation_hgvs            => { 'cols' => [qw(variation_id hgvs_name)] },
  };

  # create filenames, file handles etc
  my $tmpdir = '/tmp';
  $ImportUtils::TMP_DIR = $tmpdir;

  # create file handles
  for my $table(keys %$files) {
    
    my $hash = $files->{$table};

    $hash->{filename} = sprintf('%i_%s_%s.txt', $$, $gene->stable_id, $table);
    $hash->{filepath} = sprintf('%s/%s', $tmpdir, $hash->{filename});

    my $fh = FileHandle->new();
    $fh->open(">".$hash->{filepath}) or throw("ERROR: Could not write to ".$hash->{filepath});
    $hash->{fh} = $fh;
  }

  for my $transcript (@{ $gene->get_all_Transcripts }) {

    for my $vf(@vfs) {

      if (defined $variations_to_include) {
        next unless $variations_to_include->{$vf->variation_name};
      }

      next unless overlap($vf->start, $vf->end, $transcript->start - $max_distance, $transcript->end + $max_distance);

      my $tv = Bio::EnsEMBL::Variation::TranscriptVariation->new(
        -transcript     => $transcript,
        -variation_feature  => $vf,
        -adaptor      => $tva,
        -disambiguate_single_nucleotide_alleles => $disambiguate_sn_alleles,
        -no_transfer    => 1,
      );

      # if the variation has no effect on the transcript $tv will be undef
      if ($tv) {#} && ( scalar(@{ $tv->consequence_type }) > 0) ) {

        # store now or save to store later? Uncomment out the behaviour you want
        # save to store later uses more memory but means you don't have to sort human TV after the run
        
        ## BEHAVIOUR 1: store to DB now
        ## comment: slow as MySQL does INSERT DELAYED
        # $tva->store($tv, $mtmp);
        ## end block

        ## BEHAVIOUR 2: store to memory
        ## comment: seems like MySQL doesn't work very quickly using table locks required and also still does one INSERT per row
        # push @write_data, @{$tva->_get_write_data($tv)};
        ## end block
        
        ## BEHAVIOUR 3: store to tmp file
        ## comment: seems to be the fastest as it uses LOAD DATA which inherently locks anyway and does only one MySQL command per gene
        my $data = $tva->_get_write_data($tv);
        my $tv_fh = $files->{transcript_variation}->{fh};
        print $tv_fh join("\t", map {defined($_) ? $_ : '\N'} @$_)."\n" for @$data;

        if($mtmp) {
          my $mtmp_data = $tva->_get_mtmp_write_data_from_tv_write_data($data);
          my $mtmp_fh = $files->{MTMP_transcript_variation}->{fh};
          print $mtmp_fh join("\t", map {defined($_) ? $_ : '\N'} @$_)."\n" for @$mtmp_data;
        }
        ## end block

      
        ## populate tables for website index building
        my $var_id = $vf->get_Variation_dbID();

        for my $allele (@{ $tv->get_all_alternate_TranscriptVariationAlleles }) {

          next unless defined $allele->hgvs_transcript();
          my $hgvs_transcript = (split/\:/, $allele->hgvs_transcript())[1];
          $hgvs_by_var->{$var_id}->{$hgvs_transcript} = 1;
          
          next unless defined $allele->hgvs_protein();
          my $hgvs_protein  = (split/\:/, $allele->hgvs_protein())[1];
          $hgvs_by_var->{$var_id}->{$hgvs_protein} = 1 if defined $hgvs_protein && $hgvs_protein =~/^p/; ## don't store synonymous
        }

        # dump out these hashes periodically to stop memory exploding
        # will lead to more duplicates in the output file but better that than job failing
        if(++$var_count > 100000) {
          $self->dump_hgvs_var($hgvs_by_var, $hgvs_fh);
          $var_count = 0;
        }
      }                       
    }
  }

  ## uncomment this if using BEHAVIOUR 2 above
  # if(@write_data) {
  #   $var_dba->dbc->do("LOCK TABLES transcript_variation WRITE");
  #   $tva->_store_write_data(\@write_data, 1);
  #   $var_dba->dbc->do("UNLOCK TABLES");

  #   if($mtmp) {
  #     $var_dba->dbc->do("LOCK TABLES MTMP_transcript_variation WRITE");
  #     $tva->_store_mtmp_write_data($tva->_get_mtmp_write_data_from_tv_write_data(\@write_data), 1) if $mtmp;
  #     $var_dba->dbc->do("UNLOCK TABLES");
  #   }
  # }
  ## end block

  ## uncomment this if using BEHAVIOUR 3 above
  foreach my $table(keys %$files) {
    print STDERR "Importing data to $table\n" if $DEBUG;
    $ImportUtils::TMP_FILE = $files->{$table}->{filename};
    $files->{$table}->{fh}->close();
    load($var_dba->dbc, ($table, @{$files->{$table}->{cols}}));
  }
  ## end block

  # dump HGVS stubs to file for web index
  $self->dump_hgvs_var($hgvs_by_var, $hgvs_fh);

  $hgvs_fh->close();

  print STDERR "All done\n" if $DEBUG;

  return;
}

sub write_output {
  my $self = shift;
}

sub get_table_files_prefix {
  my $self = shift;
  my $id = shift;

  my $dir = $self->required_param('pipeline_dir').'/table_files/'.substr(md5_hex($id), 0, 2);

  unless(-d $dir) {
    mkdir($dir) or die "ERROR: Could not make directory $dir\n";
  }

  return $dir.'/'.$id;
}

sub dump_hgvs_var {
  my ($self, $hgvs_by_var, $fh) = @_;

  if($hgvs_by_var) {
    for my $var_id(keys %$hgvs_by_var) {
      print $fh "$var_id\t$_\n" for keys %{$hgvs_by_var->{$var_id}};
    }
  }

  $hgvs_by_var = {};
}

1;
