=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2023] EMBL-European Bioinformatics Institute

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
package Bio::EnsEMBL::Variation::Pipeline::UpdateVariationFeature;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

sub run {

    my $self = shift;
    my $file = $self->param('update_diff');
    my $var_dba = $self->get_species_adaptor('variation');
    $var_dba->dbc->reconnect_when_lost(1);   
    my $dbc = $var_dba->dbc;
    my $core_dba = $self->get_species_adaptor('core');
    my $sa = $core_dba->get_SliceAdaptor;
    
    # first set the default consequence type - SKIPPED if running update mode
    $dbc->do(qq{
        UPDATE  variation_feature
        SET     consequence_types = 'intergenic_variant'
    }) or die "Failed to reset consequence_types on variation_feature" if (!-e $self->param('update_diff'));

    # create a temp table (dropping it if it exists)
    my $temp_table = 'variation_feature_with_tv';

    $dbc->do(qq{DROP TABLE IF EXISTS $temp_table})
        or die "Failed to drop pre-existing temp table";
    
    $dbc->do(qq{CREATE TABLE $temp_table LIKE variation_feature})
        or die "Failed to create temp table";

    # remove unneccessary non-null columns (for EGenomes)
    $dbc->do(qq{ALTER TABLE $temp_table  drop seq_region_id,       drop variation_id ,
                                         drop seq_region_start,    drop seq_region_end,
                                         drop seq_region_strand,   drop source_id,
                                         drop map_weight })
        or die "Failed to alter temp table";

   
    # in here if we're running in update mode
    if (-e $self->param('update_diff')){
        my %genes;
        my @update_transcripts = ();
        #Grab all the genes that are new and updated
        open (DIFF, $file) or die "Can't open file $file: $!";
        while (<DIFF>){
            chomp;
            next if /^transcript_id/;
            next if /deleted/g;
            my ($transcript_id, $status, $gene_id) = split(/\t/);
            $genes{$gene_id} = 1;
        }
        close(DIFF);

        #Get slice around the genes with updated / new transcripts (go 50k up and downstream to make sure capture everything) and then store all the transcripts that are on that slice
        # and pass to update mode so all variation feature overlapping these transcripts get updated
        foreach my $gene_id (keys %genes) {
            my $slice = $sa->fetch_by_gene_stable_id($gene_id, 50000);
                foreach my $gene ( @{ $slice->get_all_Genes() } ) {
                    foreach my $transcript ( @{ $gene->get_all_Transcripts() } ) {
                        push @update_transcripts, $transcript->stable_id;
                }
            }
        }

        #Batch update variation feature for all new and updated transcripts (these have been run through transcript effect)
        #Update is slightly different than normal mode - limits update to specific feature_stable_ids - but does include all applicable transcripts, not just those that are new / updated.
        while (my @batch = splice(@update_transcripts, 0, 500) ) {
            my $joined_ids = '"' . join('", "', @batch) . '"';

            $dbc->do(qq{
                    INSERT IGNORE INTO $temp_table (variation_feature_id, consequence_types)
                    SELECT  variation_feature_id, GROUP_CONCAT(DISTINCT(consequence_types)) 
                    FROM    transcript_variation 
                    WHERE   feature_stable_id IN ($joined_ids)
                    GROUP BY variation_feature_id
                }) or die "Populating temp table failed";
        }

        # Load VF ids that overlap genes determined as deleted, then either set consequence to intergenic if VF has 
        # no other overlapping genes, or set consequences to those from overlapping non-deleted genes.
        my @old_genes;
        my $del_file = $self->param('pipeline_dir') . '/del_log/vf_affected_by_removed_transcripts.txt';
        if (! -z $del_file) {
            open (DEL, $del_file) or die "Can't open file $del_file: $!\n";
            chomp (@old_genes = <DEL>);
            close(DEL);
        }

        if(@old_genes) {
            foreach my $vf_id (@old_genes) {

                my $sth = $dbc->prepare(qq[
                        SELECT DISTINCT(consequence_types)
                        FROM transcript_variation
                        WHERE variation_feature_id = $vf_id
                     ]);
                $sth->execute();
                my $overlap = $sth->fetchall_arrayref();
                
                if (defined($overlap->[0]) ) {
                    my $consequences='';
                    for my $i (0 .. $#$overlap) {
                        my $part = $overlap->[$i];
                        for my $j (0 .. $#$part ) {
                            $consequences .= $overlap->[$i][$j] . ',';
                        }
                    }
                    my @csqs = split',',$consequences;
                    my %get_unique = map{$_ => 1 } @csqs;
                    my @unique_csqs = keys %get_unique;
                    my $csq_unique = join(",", @unique_csqs);

                    $dbc->do(qq{
                            INSERT INTO $temp_table (variation_feature_id, consequence_types)
                            VALUES ($vf_id, '$csq_unique')
                    }) or die "Populating temp table failed";
                }
        
                else{
                    $dbc->do(qq{
                            INSERT INTO $temp_table (variation_feature_id, consequence_types)
                            VALUES ($vf_id, 'intergenic_variant')
                    }) or die "Populating temp table failed";
                }
            }
        }
  }

  #Else run the standard way - i.e. full update mode
  else {
    $dbc->do(qq{
            INSERT INTO $temp_table (variation_feature_id, consequence_types)
            SELECT  variation_feature_id, GROUP_CONCAT(DISTINCT(consequence_types)) 
            FROM    transcript_variation 
            GROUP BY variation_feature_id
        }) or die "Populating temp table failed";
  }

  # update variation_feature
  $dbc->do(qq{
        UPDATE  variation_feature vf, $temp_table tvf
        SET     vf.consequence_types = tvf.consequence_types
        WHERE   vf.variation_feature_id = tvf.variation_feature_id
  }) or die "Failed to update vf table";

  # and get rid of our temp table
  $dbc->do(qq{DROP TABLE $temp_table})
      or die "Failed to drop temp table";
}

1;