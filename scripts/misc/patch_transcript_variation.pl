# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

use warnings;

use Getopt::Long;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Pipeline::TranscriptEffect;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(MAX_DISTANCE_FROM_TRANSCRIPT);

my $species;
my $reg_file;
my $transcripts_file;
my $variations_file;
my $disamb = 0;
my $verbose;
my $help;

GetOptions(
    'transcripts=s' => \$transcripts_file,
    'variations=s'  => \$variations_file,
    'species=s'     => \$species,
    'registry=s'    => \$reg_file,
    'disamb'        => \$disamb,
    'verbose|v'     => \$verbose,
    'help|h'        => \$help,
);

unless ($species && $reg_file && ($transcripts_file || $variations_file)) {
    warn "Missing required argument\n" unless $help;
    $help = 1;
}

if ($help) {
    print "Usage: $0 --species <species> --registry <file> --transcripts <file> --variations <file> --disamb --help --verbose\n";
    exit(0);
}

my $reg = 'Bio::EnsEMBL::Registry';

$reg->load_all($reg_file);

my $vdba = $reg->get_DBAdaptor($species, 'variation') or die "Failed to get variation DBA\n";

my $dbh = $vdba->dbc->db_handle;

my $remove_existing_by_tran_sth = $dbh->prepare(qq{
    DELETE FROM transcript_variation
    WHERE feature_stable_id = ?
});

my $get_affected_vfs_sth = $dbh->prepare(qq{
    SELECT  variation_feature_id
    FROM    transcript_variation
    WHERE   feature_stable_id = ?
});

my $update_vfs_sth = $dbh->prepare(qq{
    UPDATE  variation_feature vf
    SET     consequence_type = (
        SELECT  GROUP_CONCAT(tv.consequence_types)
        FROM    transcript_variation tv
        WHERE   tv.variation_feature_id = vf.variation_feature_id
    )
    WHERE vf.variation_feature_id = ?
});

my $transcript_stable_ids;

if ($transcripts_file) {
    open TRANSCRIPTS, "<$transcripts_file" or die "Failed to open $transcripts_file";

    while (<TRANSCRIPTS>) {
        next if /^\s*$/;
        chomp;
        $transcript_stable_ids->{$_} = 1;
    }
    
    close TRANSCRIPTS;
}

my $variations_by_transcript;

if ($variations_file) {
    
    my $va  = $reg->get_adaptor($species, 'variation', 'variation');
    my $vfa = $reg->get_adaptor($species, 'variation', 'variationfeature');
    my $ta  = $reg->get_adaptor($species, 'core', 'transcript');

    open VARIATIONS, "<$variations_file" or die "Failed to open $variations_file";

    while (<VARIATIONS>) {
        next if /^\s*$/;
        chomp;

        my $v = $va->fetch_by_name($_);
        
        unless ($v) {
            warn "No variation found for name $_\n";
            next;
        }

        for my $vf (@{ $vfa->fetch_all_by_Variation($v) }) {
            
            my $chr_slice = $vf->feature_Slice->expand(MAX_DISTANCE_FROM_TRANSCRIPT, MAX_DISTANCE_FROM_TRANSCRIPT);

            my @transcripts = @{ $ta->fetch_all_by_Slice($chr_slice) };

            for my $seg (@{ $chr_slice->project('lrg') }){
                push @transcripts, @{ $ta->fetch_all_by_Slice($seg->to_Slice) };
            }

            for my $t (@transcripts) {

                if ($transcripts_file) {
                    next unless $transcript_stable_ids->{$t->stable_id};
                }
                else {
                    $transcript_stable_ids->{$t->stable_id} = 1;
                }

                push @{ $variations_by_transcript->{$t->stable_id} ||= [] }, $vf;
            }
        }
    }

    close VARIATIONS;

    if ($transcripts_file && (keys %$variations_by_transcript == 0)) {
        die "Both transcripts and variations supplied, but there are no overlaps?\n";
    }
}

my $num_transcripts = scalar(keys %$transcript_stable_ids);

print "Calculating consequences for $num_transcripts transcripts\n" if $verbose;

my $affected_vf_ids;

my $num_done = 0;

for my $transcript_stable_id (keys %$transcript_stable_ids) {

    my $te = Bio::EnsEMBL::Variation::Pipeline::TranscriptEffect->new();

    $te->param('species', $species);
    $te->param('ensembl_registry', $reg_file);
    $te->param('disambiguate_single_nucleotide_alleles', $disamb);
    $te->param('transcript_stable_id', $transcript_stable_id);
    
    my $vars = $variations_by_transcript->{$transcript_stable_id};

    if ($vars) {

        my $vf_id_str = join ',', map { $_->dbID } @$vars;

        $dbh->do(qq{
            DELETE FROM transcript_variation
            WHERE feature_stable_id = '$transcript_stable_id'
            AND variation_feature_id IN ($vf_id_str)
        });

        my $vars_to_include = [ map {  $_->variation_name } @$vars ];

        $te->param('variations_to_include', $vars_to_include);
    }
    else {
        $remove_existing_by_tran_sth->execute($transcript_stable_id);
    }

    $te->run;
    
    if ($vars) {
        map { $affected_vf_ids->{$_->dbID} = 1 } @$vars;
    }
    else {
        $get_affected_vfs_sth->execute($transcript_stable_id);
    
        while (my ($vf_id) = $get_affected_vfs_sth->fetchrow_array) {
            $affected_vf_ids->{$vf_id} = 1;
        }
    }

    $num_done++;

    print "Processed $transcript_stable_id ($num_done/$num_transcripts)\n" if $verbose;
}

print "Updating ".scalar(keys %$affected_vf_ids)." affected variation features\n" if $verbose;

$num_done = 0;

for my $vf_id (keys %$affected_vf_ids) {
    $update_vfs_sth->execute($vf_id);
    $num_done++;
    if ($num_done % 1000 == 0) {
        print "Updated $num_done variations\n" if $verbose;
    }
}

