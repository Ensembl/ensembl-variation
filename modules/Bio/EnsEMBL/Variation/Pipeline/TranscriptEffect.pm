=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

  print STDERR "Fetching gene $gene_id\n" if $DEBUG;

  my $gene = $ga->fetch_by_stable_id($gene_id) 
    or die "failed to fetch gene for stable id: $gene_id";

  my $gene_name;
  $gene_name = $gene->display_xref->display_id if defined $gene->display_xref();

  my $slice = $sa->fetch_by_gene_stable_id(
    $gene->stable_id, 
    MAX_DISTANCE_FROM_TRANSCRIPT
  ) or die "failed to get slice around gene: ".$gene->stable_id;
  
  # call seq here to help cache
  $slice->seq();

  $gene = $gene->transfer($slice);

  print STDERR "Getting variation features\n" if $DEBUG;

  my @vfs = (
    @{ $slice->get_all_VariationFeatures },
    @{ $slice->get_all_somatic_VariationFeatures }
  );

  # write VFs to table file for web search indexes
  if($gene_name) {
    my $fh = FileHandle->new();
    $fh->open(">>".$self->param('pipeline_dir')."/variation_genename.txt") or die "Cannot open dump file ".$self->param('pipeline_dir')."/variation_genename.txt: $!";

    print STDERR "Getting flock for variation_hgvs\n" if $DEBUG;
    flock($fh, LOCK_EX) or die "Cannot lock - $!\n";
    # and, in case someone appended while we were waiting...
    seek($fh, 0, SEEK_END) or die "Cannot seek - $!\n";

    print STDERR "Writing to variation_hgvs file\n" if $DEBUG;
    print $fh "$_\t$gene_name\n" for map {$_->get_Variation_dbID()} @vfs;

    print STDERR "Releasing variation_hgvs flock\n" if $DEBUG;
    flock($fh, LOCK_UN) or die "Cannot unlock - $!\n";

    $fh->close();
  }


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

      next unless overlap($vf->start, $vf->end, $transcript->start - MAX_DISTANCE_FROM_TRANSCRIPT, $transcript->end + MAX_DISTANCE_FROM_TRANSCRIPT);

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
          $self->dump_hgvs_var($hgvs_by_var);
          $var_count = 0;
        }
      }                       
    }
  }

  close OUT;

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
  $self->dump_hgvs_var($hgvs_by_var);

  print STDERR "All done\n" if $DEBUG;

  return;
}

sub write_output {
  my $self = shift;
}

sub dump_hgvs_var {
  my ($self, $hgvs_by_var) = @_;

  if($hgvs_by_var) {
    my $fh = FileHandle->new();
    $fh->open(">>".$self->param('pipeline_dir')."/variation_hgvs.txt") or die "Cannot open dump file ".$self->param('pipeline_dir')."/variation_hgvs.txt: $!";

    print STDERR "Getting flock for variation_hgvs\n" if $DEBUG;
    flock($fh, LOCK_EX) or die "Cannot lock - $!\n";
    # and, in case someone appended while we were waiting...
    seek($fh, 0, SEEK_END) or die "Cannot seek - $!\n";

    print STDERR "Writing to variation_hgvs file\n" if $DEBUG;

    for my $var_id(keys %$hgvs_by_var) {
      print $fh "$var_id\t$_\n" for keys %{$hgvs_by_var->{$var_id}};
    }

    print STDERR "Releasing variation_hgvs flock\n" if $DEBUG;

    flock($fh, LOCK_UN) or die "Cannot unlock - $!\n";

    $fh->close();
  }

  $hgvs_by_var = {};
}

1;
