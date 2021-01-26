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

 This module assigns ancestral alleles.
=cut


package Bio::EnsEMBL::Variation::Utils::AncestralAllelesUtils;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw);
=head2 new

  Arg [-fasta_db] :
    Fasta DB object - Bio::DB::HTS::Faidx and Bio::DB::Fasta are supported.

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my ($fasta_db) = rearrange([qw(FASTA_DB)], @_);
  my $self = bless {
    'fasta_db' => $fasta_db,
  }, $class;

  return $self;
}

=head2 fasta_db

  Example    : $self->fasta_db;
  Description: Getter for fasta db object
  Returntype : Bio::DB::HTS::Faidx or Bio::DB::Fasta
  Exceptions : None
  Caller     : sequence_id_mappings, assign 

=cut

sub fasta_db {
  my $self = shift;
  return $self->{'fasta_db'}; 
}

=head2 sequence_id_mappings

  Example    : $self->sequence_id_mappings;
  Description: Getter for sequence id mappings. The sequence ids from the ancestral
               fasta file look like ANCESTOR_for_chromosome:GRCh38:4:1:190214555:1.
               We need to extract the chromosome name and store the mapping between
               fasta file sequence id and chromosome name for look up in the assign
               ancestral allele step.
  Returntype : Hashref
  Exceptions : Throw error if the fasta_db object is neither of type
               Bio::DB::HTS::Faidx or Bio::DB::Fasta.
               Throw error if the sequence ids have changed and don't follow the 
               expected pattern of colon separated values.
  Caller     : get_fasta_sequence_id 

=cut

sub sequence_id_mappings {
  my $self = shift;
  if (!defined $self->{'sequence_id_mappings'}) {
    my $fasta_db = $self->fasta_db;
    my @sequence_ids = ();
    if ($fasta_db->isa('Bio::DB::HTS::Faidx')) {
      @sequence_ids = $fasta_db->get_all_sequence_ids;
    } elsif ($fasta_db->isa('Bio::DB::Fasta')) {
      @sequence_ids = $fasta_db->get_all_ids;
    } else {
      throw("ERROR: Couldn't get sequence ids from ".ref($fasta_db)."\n");
    }

    my $sequence_id_2_chr_number;

    foreach my $sequence_id (@sequence_ids) {
      throw("ERROR: sequence ids have changed and don't follow the expected pattern of colon separated values\n") if ($sequence_id !~ m/:/);
      # expected something like ANCESTOR_for_chromosome:GRCh38:10:1:133797422:1
      my @split = split(/:/, $sequence_id);
      throw("ERROR: sequence ids have changed and don't follow the expected pattern of 6 colon separated values\n") if (scalar @split != 6);
      my $chrom = $split[2];
      throw("ERROR: undefined chromosome in $sequence_id\n") if (!$chrom);
      $sequence_id_2_chr_number->{$chrom} = $sequence_id;
    }
    $self->{'sequence_id_mappings'} = $sequence_id_2_chr_number;
  }
  return $self->{sequence_id_mappings};
}

=head2 get_fasta_sequence_id

  Arg [1]    : String $chromosome_name
  Example    : my $sequence_id = $self->get_fasta_sequence_id(4);
  Description: Getter for sequence fasta id for a given chromosome name
  Returntype : String
  Exceptions : None
  Caller     : assign

=cut

sub get_fasta_sequence_id {
  my $self = shift;
  my $chrom = shift;
  my $sequence_id_mappings = $self->sequence_id_mappings;
  return $sequence_id_mappings->{$chrom};  
}

=head2 assign

  Arg [1]    : String $chromosome_name
  Arg [2]    : String $start
  Arg [3]    : String $end
  Example    : my $ancestral_allele = $ancestral_alleles_utils->assign($chrom, $start, $end);
  Description: Assigns the ancestral allele for the given chromosome, start and end.
               Assigns ancestral allele if:
                - defined region is smaller than or equal to 50bp
                - defined region is not describing an insertion with start > end
                - returned ancestral allele only contains ACGT characters
                - chromosome is not an alternate sequence and contained in the 
                  ancestral fasta file
  Returntype : String $ancestral_allele if ancestral allele could be assigned
               else undef.
  Exceptions : Throw error if the fasta_db object is neither of type
               Bio::DB::HTS::Faidx or Bio::DB::Fasta.
  Caller     : ensembl-variation/scripts/import/dbSNP_v2/load_dbsnp.pl
               Bio::EnsEMBL::Variation::Pipeline::AncestralAlleles::Assign 

=cut

sub assign {
  my ($self, $chrom, $start, $end) = @_;

  return undef if (!($start && $end));

  # insertion
  return undef if ($start > $end);
  # allele size limit to 50bp
  return undef if (($end - $start) >= 50); 

  # alternative sequences are not represented in the ancestral fasta file 
  my $fasta_sequence_id = $self->get_fasta_sequence_id($chrom);
  return undef if (!$fasta_sequence_id);

  my $ancestral_allele = undef;

  my $fasta_db = $self->fasta_db;
  if ($fasta_db->isa('Bio::DB::HTS::Faidx') ) {
    $ancestral_allele = $fasta_db->get_sequence_no_length("$fasta_sequence_id:$start-$end");
  } elsif ($fasta_db->isa('Bio::DB::Fasta')) {
    $ancestral_allele = $fasta_db->seq("$fasta_sequence_id:$start,$end");
  } else {
    throw("ERROR: Don't know how to fetch sequence from a ".ref($fasta_db)."\n");
  }

  return undef unless ($ancestral_allele && $ancestral_allele =~ m/^[ACGT]+$/i);
  
  $ancestral_allele = uc $ancestral_allele;

  return $ancestral_allele;

}

1;
