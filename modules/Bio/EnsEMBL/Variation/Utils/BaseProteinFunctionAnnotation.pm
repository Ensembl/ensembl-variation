=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

=head1 DESCRIPTION

This module contains functionality for adding scores and predictions
for amino acid changes in a translation from file.

=cut

package Bio::EnsEMBL::Variation::Utils::BaseProteinFunctionAnnotation;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::Tools::CodonTable;
use Bio::DB::HTS::Tabix;
use Bio::EnsEMBL::Variation::Utils::BaseDatabaseUtils;
use FileHandle;
use Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix qw(@ALL_AAS);

our @ISA = ('Bio::EnsEMBL::Variation::Utils::BaseDatabaseUtils');

my $LOW_QUALITY = 0;

=head2 new

  Arg [-working_dir] :
    string - location of the working directory when running the pipeline. The directory is used for storing debug information.
  Arg [-annotation_file] :
    string - location of dbNSFP or CADD file
  Arg [-annotation_file_version] :
    string - version of annotation file
  Arg [-assembly] :
    string - assembly version e.g. GRCh37 or GRCh38
  Arg [-pipeline_mode] :
    boolean - If set to 1 run in pipeline mode and store results in the database. The default value is 1.
  Arg [-debug_mode] :
    boolean - If set to 1 write debug information to the working directory. The default value is 0.
  Arg [-write_mode] :
    boolean - If set to 1 write error file which reports missmatch between protein sequence and translated protein sequence as a result of annotation.
=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);
  my ($working_dir, $annotation_file, $annotation_file_version, $assembly, $pipeline_mode, $debug_mode, $write_mode) = rearrange([qw(WORKING_DIR ANNOTATION_FILE ANNOTATION_FILE_VERSION ASSEMBLY PIPELINE_MODE DEBUG_MODE WRITE_MODE)], @_);
  $self->{'working_dir'} = $working_dir;
  $self->{'annotation_file'} = $annotation_file;
  $self->{'annotation_file_version'} = $annotation_file_version;
  $self->{'assembly'} = $assembly;
  $self->{'pipeline_mode'} = (!defined $pipeline_mode) ? 1 : $pipeline_mode; # default set to 1
  $self->{'write_mode'} = (!defined $write_mode) ? 1 : $write_mode; # default set to 1
  $self->{'debug_mode'} = $debug_mode;

  if (! grep {$_ eq $self->assembly} ('GRCh37', 'GRCh38')) {
    die "Assembly $assembly is not supported.";
  }

  return $self;
}

=head2 run
  Arg 1      : String $translation_md5  
  Arg 2      : (optional) Hashref of translation mappings for testing.
  Example    : $self->run('4d08f77b4cb14259684ce086ba089565');
               $self->run('4d08f77b4cb14259684ce086ba089565', {'4d08f77b4cb14259684ce086ba089565' => 'ENSP00000435699'});
  Description: Runs protein function prediction annotations for a protein translation.
  Returntype : none
  Exceptions : throws on missing argument
               throws on undefined $translation_stable_id
  Caller     : Bio::EnsEMBL::Variation::Utils::CADDProteinFunctionAnnotation::run() 
               Bio::EnsEMBL::Variation::Utils::DbNSFPProteinFunctionAnnotation::run()
  Status     : 
=cut

sub run {
  my $self = shift;
  my $translation_md5 = shift;
  my $translation_mappings = shift;

  throw("Translation_md5 string expected") if (!defined $translation_md5);
  my $translation_stable_id = (defined $translation_mappings) ? $translation_mappings->{$translation_md5} : $self->get_stable_id_for_md5($translation_md5);
  throw("No translation_stable_id for translation_md5 $translation_md5") if (!defined $translation_stable_id);

  my $translation = $self->get_translation($translation_stable_id);
  my $translation_seq = $translation->seq;
  my $transcript = $translation->transcript;
  $self->reverse($transcript->strand < 0);
  my $transcript_stable_id = $transcript->stable_id;

  $self->init_protein_matrix($translation, $translation_md5);

  $self->init_header;

  my $all_triplets = $self->get_triplets($translation_stable_id);

  $self->load_predictions_for_triplets($all_triplets);

  $self->store_protein_matrix($translation_stable_id, $translation_md5) if ($self->{'pipeline_mode'});

  if ($self->{'write_mode'}) {
    if ($translation_seq ne join('', @{$self->amino_acids})) {
      my $fh = FileHandle->new($self->working_dir. "/$translation_stable_id", 'w');
      print $fh "$transcript_stable_id\n$translation_seq\n";
      print $fh join('', @{$self->amino_acids}), "\n";
      $fh->close;
    }
  }
}

=head2 amino_acids
  Arg [1]    : String $aa Amino acid  
  Example    : $runCADDAnnotation->amino_acids('K');
  Description: Holds reference to an array of string. This method is used to collect each
               annotated amino acid. The final array will be compared against the input translation.
  Returntype : Arrayref of string
  Exceptions : None
  Caller     : General
  Status     : At Risk
=cut
sub amino_acids {
  my $self = shift;
  my $aa = shift;
  if (defined $aa) {
    push @{$self->{'amino_acids'}}, $aa;
  } else {
    return $self->{'amino_acids'};
  }
}

=head2 analysis
  Arg 1      : Arrayref of string analysis (optional)  
  Example    : $self->analysis([qw/dbnsfp_revel dbnsfp_meta_lr dbnsfp_mutation_assessor/]);
  Description: Set and get available analysis. 
  Returntype : 
  Exceptions : None
  Caller     : General
  Status     : At Risk
=cut
sub analysis {
  my $self = shift;
  return $self->{'analysis'} = shift if(@_);
  return $self->{'analysis'};
}

=head2 annotation_file
  Description: Get annotation file
  Returntype : String $annotation_file
  Exceptions : None
  Caller     : General
  Status     : At Risk
=cut
sub annotation_file {
  my $self = shift;
  return $self->{'annotation_file'};
}

=head2 annotation_file_version
  Description: Get annotation file version
  Returntype : String $annotation_file_version
  Exceptions : None
  Caller     : General
  Status     : At Risk
=cut
sub annotation_file_version {
  my $self = shift;
  return $self->{'annotation_file_version'};
}

=head2 assembly
  Description: Get assembly 
  Returntype : String $assembly
  Exceptions : None
  Caller     : General
  Status     : At Risk
=cut
sub assembly {
  my $self = shift;
  return $self->{'assembly'};
}

=head2 reverse
  Arg 1      : (optional) Boolean $reverse  
  Description: Get and set reverse flag 
  Returntype : Boolean $reverse
  Exceptions : None
  Caller     : General
  Status     : At Risk
=cut
sub reverse {
  my $self = shift;
  return $self->{'reverse'} = shift if(@_);
  return $self->{'reverse'};
}

=head2 working_dir
  Description: Get working directory
  Returntype : String $working_dir 
  Exceptions : None
  Caller     : General
  Status     : At Risk
=cut
sub working_dir {
  my $self = shift;
  return $self->{'working_dir'};
}

=head2 header
  Description: - Get arrayref of header row 
               - Initialise at first use
  Returntype : Arrayref $header 
  Exceptions : None
  Caller     : Bio::EnsEMBL::Variation::Utils::CADDProteinFunctionAnnotation::get_CADD_row 
               Bio::EnsEMBL::Variation::Utils::DbNSFPProteinFunctionAnnotation::get_dbNSFP_row
  Status     : At Risk
=cut
sub header {
  my $self = shift;
  return $self->{'header'} = shift if(@_);
  $self->init_header if (!defined $self->{'header'});
  return $self->{'header'};
}

=head2 init_header
  Description: Initialise header row as arrayref 
  Returntype : None 
  Exceptions : None
  Caller     : header()
  Status     : 
=cut
sub init_header {
  my $self = shift;
  my $header;
  my $annotation_file = $self->annotation_file;
  open HEAD, "tabix -fh $annotation_file 1:1-1 2>&1 | ";
  while(<HEAD>) {
    next unless /^\#/;
    chomp;
    $header = [split];
  }
  close HEAD;

  $self->header($header);
}

=head2 parser
  Example    : $self->parser->query("$chrom:$start-$end"); 
  Description: Get parser for annotation file 
  Returntype : Bio::DB::HTS::Tabix $parser 
  Exceptions : none
  Caller     : get_tabix_iterator() 
  Status     : 
=cut
sub parser {
  my $self = shift;
  if (!defined $self->{'parser'}) {
    my $annotation_file = $self->annotation_file;
    $self->{'parser'} = Bio::DB::HTS::Tabix->new(filename => $annotation_file);
  }
  return $self->{'parser'};
}

=head2 get_tabix_iterator
  Arg 1      : Int $chrom 
  Arg 2      : Int $start  
  Arg 3      : Int $end 
  Example    : my $iter = $self->get_tabix_iterator($chrom, $triplet_start, $triplet_end); 
  Description: Get iterator over provided region $chrom, $start, $end 
  Returntype : Bio::DB::HTS::Tabix::Iterator $iter 
  Exceptions : none
  Caller     : load_predictions_for_triplets() 
  Status     : 
=cut
sub get_tabix_iterator {
  my ($self, $chrom, $start, $end) = @_;
  return $self->parser->query("$chrom:$start-$end");
}

=head2 codon_table
  Description: Get Bio::Tools::CodonTable object
  Returntype : Bio::Tools::CodonTable $codon_table
  Exceptions : None
  Caller     : get_triplets() 
  Status     : 
=cut
sub codon_table {
  my $self = shift;
  if (!defined $self->{'codon_table'}) {
    $self->{'codon_table'} = Bio::Tools::CodonTable->new();
  }
  return $self->{'codon_table'};
}

=head2 get_stable_id_for_md5
  Arg 1      : String $translation_md5 
  Description: Get translation stable for md5 from translation_mapping table
  Returntype : String $translation_stable_id  
  Exceptions : none
  Caller     : run() 
  Status     : 
=cut
sub get_stable_id_for_md5 {
  my ($self, $md5) = @_;
  my $var_dba = $self->get_species_adaptor('variation');
  my $get_stable_id_sth = $var_dba->dbc->prepare(qq{
    SELECT  stable_id
    FROM    translation_mapping
    WHERE   md5 = ?
  });

  $get_stable_id_sth->execute($md5);
  my ($stable_id) = $get_stable_id_sth->fetchrow_array;
  return $stable_id;
}

=head2 get_translation
  Arg 1      : String $translation_stable_id 
  Description: - Get translation object 
               - Choose core database based on translation prefix
  Returntype : Bio::EnsEMBL::Translation $translation
  Exceptions : throw on missing value
  Caller     : run() 
  Status     : 
=cut
sub get_translation {
  my $self = shift;
  my $translation_stable_id = shift;
  throw("Translation_stable_id string expected") if (!defined $translation_stable_id);
  my $core_type = ($translation_stable_id =~ /^NP|XP/) ? 'otherfeatures' : 'core';
  my $cdba = $self->get_species_adaptor($core_type);
  my $translation_adaptor = $cdba->get_TranslationAdaptor or die "Failed to get translation adaptor";
  my $translation = $translation_adaptor->fetch_by_stable_id($translation_stable_id);
  return $translation;
}

=head2 mutate
  Arg 1      : String $triplet  
  Arg 2      : Boolean $reverse
  Example    : 
  Description: Mutate all positions of the input triplet sequence 
  Returntype : Hashref $triplet with all possible mutation e.g. for triplet ATG:
                { 'ATG' => {
                           '1' => {
                                    'A' => 'AAG',
                                    'T' => 'ATG',
                                    'C' => 'ACG',
                                    'G' => 'AGG'
                                  },
                           '0' => {
                                    'A' => 'ATG',
                                    'T' => 'TTG',
                                    'C' => 'CTG',
                                    'G' => 'GTG'
                                  },
                           '2' => {
                                    'A' => 'ATA',
                                    'T' => 'ATT',
                                    'C' => 'ATC',
                                    'G' => 'ATG'
                                  }
                         }
                  }
  Exceptions : None
  Caller     : get_triplets() 
  Status     : 
=cut
sub mutate {
  my $self = shift;
  my $triplet = shift;
  my $reverse = shift;
  my @nucleotides = split('', $triplet);
  my $new_triplets;
  foreach my $i (0 .. $#nucleotides) {
    $new_triplets = $self->get_mutated_triplets($triplet, $i, $new_triplets, $reverse);
  }
  return $new_triplets;
}

=head2 get_mutated_triplets
  Arg 1      : String $triplet   
  Arg 2      : Int $nucleotide_position 
  Arg 3      : Hashref $new_triplets 
  Arg 4      : Boolean reverse
  Example    : 
  Description: Generate all mutations for a given triplet position 
  Returntype : 
  Exceptions : None
  Caller     : mutate()
  Status     : At Risk
=cut
sub get_mutated_triplets {
  my $self = shift;
  my $triplet = shift;
  my $position = shift;
  my $new_triplets = shift;
  my $reverse = shift;
  my $mutations = ['A', 'G', 'C', 'T'];
  foreach my $mutation (@$mutations) {
    my $update_triplet = $triplet;
    if ($reverse) {
      my $reverse_mutation = $mutation;
      reverse_comp(\$reverse_mutation);
      substr($update_triplet, $position, 1, $reverse_mutation);
    } else {
      substr($update_triplet, $position, 1, $mutation);
    }
    $new_triplets->{$triplet}->{$position}->{$mutation} = $update_triplet;
  }
  return $new_triplets;
}

=head2 get_triplets
  Arg 1      : String $translation_stable_id  
  Example    : 
  Description: - get translation object for translation_stable_id
               - get corresponding transcript object for translation object
               - iterate over each amino acid of the translation object and
                 use the transcript_mapper to get pep2genomic coordinates     
               - foreach set of coordinates for an amino acid build the 
                 triplet sequence and store the genomic start and end
                 coordinates of the triplet 
               - create an entry of coords, aa_position, chrom and triplet_seq
                 for each triplet
               - add the amino acid (aa) to each entry by using the codon table to
                 translate the triplet to amino acid
               - add all possible mutations (new_triplets) for a triplet using mutate() 
  Returntype : arrayref of hashes. Each hash has the following keys: triplet_seq, coords,
               aa_position, aa_position, new_triplets, aa 
  Exceptions : None
  Caller     : run()
  Status     : At Risk
=cut
sub get_triplets {
  my $self = shift;
  my $translation_stable_id = shift;
  my $translation = $self->get_translation($translation_stable_id);
  my $core_type = ($translation_stable_id =~ /^NP|XP/) ? 'otherfeatures' : 'core';
  my $cdba = $self->get_species_adaptor($core_type);
  my $slice_adaptor = $cdba->get_SliceAdaptor or die "Failed to get slice adaptor";
  my $transcript = $translation->transcript;
  my $chrom = $transcript->seq_region_name;
  my $start = $transcript->seq_region_start;
  my $end = $transcript->seq_region_end;
  my $strand = $transcript->seq_region_strand;
  my $slice = $slice_adaptor->fetch_by_region('toplevel', $chrom,  $start, $end);
  my $transcript_mapper = $transcript->get_TranscriptMapper();

  my $codon_table = $self->codon_table;
  my @all_triplets = ();
  foreach my $i (1 .. $translation->length) {
    my @pep_coordinates = $transcript_mapper->pep2genomic($i, $i);
    my $triplet = '';
    my @coords = ();
    foreach my $coord (@pep_coordinates) {
      my $coord_start = $coord->start;
      my $coord_end = $coord->end;
      next if ($coord_start <= 0);
      my $new_start = $coord_start - $start + 1;
      my $new_end   = $coord_end   - $start + 1;
      my $subseq = $slice->subseq($new_start, $new_end, $strand);
      $triplet .= $subseq;
      push @coords, [$coord_start, $coord_end];
    }
    my $entry = {
      coords => \@coords,
      aa_position => $i,
      chrom => $chrom,
      triplet_seq => $triplet,
    };
    my $aa = $codon_table->translate($triplet);
    if (!$aa) {
      $entry->{aa} = 'X';
    } else {
      $entry->{aa} = $aa;
      my $reverse = ($strand < 0);
      my $new_triplets = $self->mutate($triplet, $reverse);
      $entry->{new_triplets} = $new_triplets;
    }
    push @all_triplets, $entry;
  } 
  return \@all_triplets;
}

=head2 init_protein_matrix
  Arg 1      : Bio::EnsEMBL::Translation $translation
  Arg 2      : String $translation_md5
  Example    : 
  Description: Initialise ProteinFunctionPredictionMatrix for each analysis.
  Returntype : 
  Exceptions : none
  Caller     : run()
  Status     : 
=cut
sub init_protein_matrix {
  my $self = shift;
  my $translation = shift;
  my $translation_md5 = shift;
  my $pred_matrices = {};
  foreach my $analysis (@{$self->analysis}) {
    my $pred_matrix = Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix->new(
      -analysis       => $analysis,
      -peptide_length   => $translation->length,
      -translation_md5  => $translation_md5,
    );
    $pred_matrices->{$analysis} = $pred_matrix;
  }

  $self->{pred_matrices} = $pred_matrices;
}

=head2 load_predictions_for_triplets
  Arg 1      : Arrayref of triplet hashes
  Description: Get and store score and prediction for each mutated triplet
               Override in child class
  Returntype : None
  Exceptions : None
  Caller     : run()
  Status     :
=cut
sub load_predictions_for_triplets {
  return;
}

=head2 add_predictions
  Arg 1      : Hashref mapping header column to row value
  Arg 2      : Int $i amino acid position
  Arg 3      : String $mutated_aa mutated amino acid
  Description: Store score and prediction for mutated amino acid
               if score exists.
               Override in child class 
  Returntype : None
  Exceptions : None
  Caller     : load_predictions_for_triplets()
  Status     :
=cut
sub add_predictions {
  return;
}

=head2 add_prediction
  Arg 1      : Int - location of amino acid in translation 
  Arg 2      : String - amino acid 
  Arg 3      : String - analysis e.g. dbnsfp_mutation_assessor, cadd
  Arg 4      : Double - score
  Arg 5      : String - prediction 
  Example    : 
  Description: Collect all prediction results. 
  Returntype : 
  Exceptions : None
  Caller     : Bio::EnsEMBL::Variation::Utils::RunCADDAnnotationUtils::add_predictions() 
               Bio::EnsEMBL::Variation::Utils::RunDbNSFPAnnotationUtils::add_predictions() 
  Status     : 
=cut
sub add_prediction {
  my ($self, $i, $mutated_aa, $analysis, $score, $prediction) = @_;
  $self->{pred_matrices}->{$analysis}->add_prediction(
    $i,
    $mutated_aa,
    $prediction,
    $score,
    $LOW_QUALITY,
  );

  $self->{results_available}->{$analysis} = 1;
  $self->{debug_data}->{$analysis}->{$i}->{$mutated_aa}->{$prediction} = $score if ($self->{'debug_mode'});
}

=head2 store_protein_matrix
  Arg 1      : String translation_stable_id
  Arg 2      : String translation_md5
  Example    : $self->store_protein_matrix('ENSP00000435699', '4d08f77b4cb14259684ce086ba089565');
  Description: Store ProteinFunctionPredictionMatrix for each analysis and translation if results are available 
  Returntype : none
  Exceptions : none
  Caller     : run()
  Status     : 
=cut
sub store_protein_matrix {
  my $self = shift;
  my $translation_stable_id = shift;
  my $translation_md5 = shift;
  my $pred_matrices = $self->{pred_matrices};

  my $vdba = $self->get_species_adaptor('variation');
  my $pfpma = $vdba->get_ProteinFunctionPredictionMatrixAdaptor or die "Failed to get matrix adaptor";

  foreach my $analysis (keys %$pred_matrices) {
    my $pred_matrix = $pred_matrices->{$analysis};
    if ($self->{results_available}->{$analysis}) {
      $pfpma->store($pred_matrix);
      if ($self->{'debug_mode'} && $self->{'write_mode'}) {
        my $write_file = $self->working_dir. "/$analysis\_$translation_stable_id";
        my $fh = FileHandle->new($write_file, 'w') or die "cannot open file $write_file: $!";
        my $matrix = $pfpma->fetch_by_analysis_translation_md5($analysis, $translation_md5);
        my $debug_data = $self->{debug_data};
        foreach my $i (sort keys %{$debug_data->{$analysis}}) {
          foreach my $aa (keys %{$debug_data->{$analysis}->{$i}}) {
            next if ($aa eq '*');
            foreach my $prediction (keys %{$debug_data->{$analysis}->{$i}->{$aa}}) {
              my ($new_pred, $new_score) = $matrix->get_prediction($i, $aa);
              print $fh join(' ', $analysis, $i, $aa, $prediction, $debug_data->{$analysis}->{$i}->{$aa}->{$prediction}, $new_pred, $new_score), "\n";
            }
          }
        }
        $fh->close;
      }
    }
  }
}

1;
