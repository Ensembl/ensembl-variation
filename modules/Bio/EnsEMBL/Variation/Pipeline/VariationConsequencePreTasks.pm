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

=cut

package Bio::EnsEMBL::Variation::Pipeline::VariationConsequencePreTasks;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);

use File::Path qw(mkpath rmtree);

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

sub fetch_input {
   
    my $self = shift;

    my $mtmp = $self->param('mtmp_table');
    

    # check for out of date seq_regions in variation database
    my $sequences_ok = $self->check_seq_region();
    die "Seq region ids are not compatible. Run ensembl-variation/scripts/misc/update_seq_region_ids.pl\n" unless $sequences_ok == 1;

    my $core_dba = $self->get_species_adaptor('core');
    my $var_dba = $self->get_species_adaptor('variation');
    
    my $dbc = $var_dba->dbc();

    $dbc->do(qq{
      ALTER TABLE variation_feature
      ORDER BY seq_region_id, seq_region_start, seq_region_end
    }) if $self->param('sort_variation_feature');


      # truncate the table because we don't want duplicates
      $dbc->do("TRUNCATE TABLE transcript_variation");

      # disable the indexes on the table we're going to insert into as
      # this significantly speeds up the TranscriptEffect process

      $dbc->do("ALTER TABLE transcript_variation DISABLE KEYS");

      # truncate tables incase TranscriptVariation is being updated for a pre-existing database
      $dbc->do("TRUNCATE TABLE variation_hgvs");
      $dbc->do("ALTER TABLE variation_hgvs DISABLE KEYS");

      # remove temporary files if they exist
      my $dir = $self->param('pipeline_dir');
      unless(-d $dir) {
        mkpath($dir) or die "ERROR: Could not create directory $dir (required for dump files)\n";
      }

      foreach my $folder_name (qw/web_index transcript_effect load_log/) {
        rmtree($dir.'/'.$folder_name.'_files');
        mkdir($dir.'/'.$folder_name.'_files') or die "ERROR: Could not make directory $dir/$folder_name\_files\n";
      }

      my @rebuild = qw(transcript_variation variation_hgvs);

      # set up MTMP table
      if($mtmp) {
        my @exclude = qw(transcript_variation_id hgvs_genomic hgvs_protein hgvs_transcript somatic codon_allele_string);
        my ($source_table, $table) = qw(transcript_variation MTMP_transcript_variation);

        # drop existing MTMP
        $dbc->do(qq{DROP TABLE IF EXISTS $table});

        my $sth = $dbc->prepare(qq{
          SHOW CREATE TABLE $source_table
        });
        $sth->execute();

        my $create_sth = $sth->fetchall_arrayref->[0]->[1];
        $sth->finish;

        # convert set to enum
        $create_sth =~ s/^set/enum/;

        # rename table
        $create_sth =~ s/TABLE \`$source_table\`/TABLE \`$table\`/;

        # filter out some columns
        $create_sth =~ s/\`?$_.+?,// for @exclude;

        # filter out some indices
        $create_sth =~ s/AUTO_INCREMENT=\d+//;
        $create_sth =~ s/somatic_feature_idx/feature_idx/;
        $create_sth =~ s/$_.+,// for ('PRIMARY KEY', 'KEY `somatic', 'KEY `cons');
        $create_sth =~ s/,\`somatic\`//;

        # remove final comma
        $create_sth =~ s/,(\s+\))/$1/;

        $dbc->do($create_sth);
        $dbc->do("ALTER TABLE $table DISABLE KEYS");

        push @rebuild, $table;

        # setup fasta
        if (my $fasta = $self->param('fasta')) {

          # run this here as it generates the index
          # don't want competing jobs writing to it
          setup_fasta(-FASTA => $fasta);
        }
    }
    $self->param(
        'rebuild_indexes', [{
            tables => \@rebuild,
        }]
    );

}

sub write_output {
  my $self = shift;
  $self->dataflow_output_id($self->param('rebuild_indexes'), 2);
  return;
}

## check for out of date seq_regions in variation database
## human patches can change between releases. such differences break TranscriptEffect
sub check_seq_region{

   my $self = shift;

   my $stmt = qq[ select seq_region_id, name from seq_region];

   my $core_dba = $self->get_species_adaptor('core');
   my $var_dba  = $self->get_species_adaptor('variation');

   my $core_seq_sth = $core_dba->dbc->prepare($stmt);
   $core_seq_sth->execute();
   my $core_ids = $core_seq_sth->fetchall_arrayref();

   my %expected_ids;
   foreach my $l(@{$core_ids}){
       $expected_ids{$l->[0]} = $l->[1];
   }

   my $var_seq_sth = $var_dba->dbc->prepare($stmt);
   $var_seq_sth->execute();
   my $var_ids = $var_seq_sth->fetchall_arrayref();

   my $OK = 1;
   foreach my $l(@{$var_ids}){
       unless (defined $expected_ids{$l->[0]} ){
           $self->warning( 'Seq_region_id in variation db is not in core: '.  $l->[0]. ' ' . $l->[1]);
           $OK = 0;
       }
   }
   return $OK;
}

1;
