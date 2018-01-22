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
package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::InitJobs;

use strict;

use Bio::EnsEMBL::Registry;

use Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::Constants qw(FULL UPDATE NONE);

use Digest::MD5 qw(md5_hex);

use File::Path qw(make_path);

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

sub fetch_input {
   
    my $self = shift;
    
    my $sift_run_type   = $self->required_param('sift_run_type');
    my $pph_run_type    = $self->required_param('pph_run_type');
    my $include_lrg     = $self->param('include_lrg');
    
    $self->update_meta ;


    my $core_dba = $self->get_species_adaptor('core');
    my $var_dba  = $self->get_species_adaptor('variation');

    # fetch all the transcripts from the core DB

    my @transcripts;

    if ($self->param('debug_mode')) {
        my $ga = $core_dba->get_GeneAdaptor or die "Failed to get gene adaptor";
        
        @transcripts = grep { $_->translation } @{ $ga->fetch_all_by_external_name('BRCA1')->[0]->get_all_Transcripts };
    }
    else {
        my $sa = $core_dba->get_SliceAdaptor or die "Failed to get slice adaptor";
        
        for my $slice (@{ $sa->fetch_all('toplevel', undef, 1, undef, ($include_lrg ? 1 : undef)) }) {
            for my $gene (@{ $slice->get_all_Genes(undef, undef, 1) }) {
                for my $transcript (@{ $gene->get_all_Transcripts }) {
                    if (my $translation = $transcript->translation) {
                        push @transcripts, $transcript;
                    }
                }
            }
        }
    }

    push @transcripts, @{$self->get_refseq_transcripts} if $self->param('include_refseq');

    # store a table mapping each translation stable ID to its corresponding MD5

    $var_dba->dbc->do(qq{DROP TABLE IF EXISTS translation_mapping});

    $var_dba->dbc->do(qq{
        CREATE TABLE translation_mapping (
            stable_id   VARCHAR(255),
            md5         CHAR(32),
            PRIMARY KEY (stable_id),
            KEY md5_idx (md5)
        )
    });

    my $add_mapping_sth = $var_dba->dbc->prepare(qq{
        INSERT IGNORE INTO translation_mapping (stable_id, md5) VALUES (?,?)
    });

    my $add_md5_sth = $var_dba->dbc->prepare(qq{
        INSERT IGNORE INTO translation_md5 (translation_md5) VALUES (?)
    });

    # build a hash mapping MD5s for our set of translations to their peptide sequences
    # and also write each stable ID - MD5 mapping to the database

    my %unique_translations;

    for my $tran (@transcripts) {
        
        my $tl = $tran->translation;

        my $seq = $tl->seq;
        
        my $md5 = md5_hex($seq);
        
        $unique_translations{$md5} = $seq;

        $add_mapping_sth->execute($tl->stable_id, $md5);
    }

    # work out which translations we need to run for which analysis
    
    my @translation_md5s = keys %unique_translations;
    
    my @sift_md5s;
    my @pph_md5s;

    # if we're doing full runs, then we just use all the translations

    @sift_md5s = @translation_md5s if $sift_run_type == FULL;
    @pph_md5s  = @translation_md5s if $pph_run_type == FULL;

    # if we're updating we need to check which translations already have predictions
        
    if ($sift_run_type == UPDATE || $pph_run_type == UPDATE) {
        
        my $var_dbh = $var_dba->dbc->db_handle;
    
        my $get_existing_sth = $var_dbh->prepare(qq{
            SELECT  t.translation_md5, a.value, p.prediction_matrix
            FROM    translation_md5 t, attrib a, protein_function_predictions p
            WHERE   t.translation_md5_id = p.translation_md5_id 
            AND     a.attrib_id = p.analysis_attrib_id
        });

        # store the set of existing MD5s in a hash

        $get_existing_sth->execute;

        my $existing_md5s;

        while ( my ($md5, $analysis, $preds) = $get_existing_sth->fetchrow_array ) {
            # there are 2 polyphen analyses, but we only want to track one
            $analysis = 'pph' if $analysis =~ /polyphen/;
            # just record true if we already have predictions for each translation
            $existing_md5s->{$md5}->{$analysis} = length($preds) > 0;   
        }

        # exclude any translations MD5s we find in the database

        @sift_md5s = grep { ! $existing_md5s->{$_}->{sift} } @translation_md5s if $sift_run_type == UPDATE;
        @pph_md5s  = grep { ! $existing_md5s->{$_}->{pph} } @translation_md5s if $pph_run_type == UPDATE;
    }
    
    # create a FASTA dump of all the necessary translation sequences

    # work out the set of all unique MD5s for sift and polyphen

    my %required_md5s = map { $_ => 1 } (@sift_md5s, @pph_md5s);

    my $fasta = $self->required_param('fasta_file');

    my @dir = split('/', $fasta);
    pop @dir;
    make_path(join('/', @dir));

    open my $FASTA, ">$fasta" or die "Failed to open $fasta for writing";

    # get rid of any existing index file

    if (-e "$fasta.fai") {
        unlink "$fasta.fai" or die "Failed to delete fasta index file";
    }
    
    # dump out each peptide sequence

    for my $md5 (keys %required_md5s) {
        my $seq = $unique_translations{$md5};
        $seq =~ s/(.{80})/$1\n/g;
        # get rid of any trailing newline
        chomp $seq;
        print $FASTA ">$md5\n$seq\n";
    }

    close $FASTA;

    # set up our list of output ids

    $self->param('pph_output_ids',  [ map { {translation_md5 => $_} } @pph_md5s ]);
    $self->param('sift_output_ids', [ map { {translation_md5 => $_} } @sift_md5s ]);
}

## hold code & protein database version in meta table if new complete run
sub update_meta{
    my $self = shift;

    my $var_dba  = $self->get_species_adaptor('variation');

    my $var_dbh = $var_dba->dbc->db_handle;
    
    my $update_meta_sth = $var_dbh->prepare(qq{
            insert ignore into meta ( meta_key, meta_value) values (?,?)
        });

    if ($self->required_param('sift_run_type')  == FULL){

	my @code =split/\//, $self->required_param('sift_dir');
	my $sift_version = pop @code;

	$update_meta_sth->execute('sift_version', $sift_version);

	unless($self->param('use_compara')){

	    my @db =split/\//, $self->required_param('blastdb');
	    my $db_version = pop @db;

	    $update_meta_sth->execute('sift_protein_db_version', $db_version);
	}

    }
    if ($self->required_param('pph_run_type')  == FULL){

	my @code =split/\//,$self->required_param('polyphen_dir');
	my $polyphen_version = pop @code;   

	$update_meta_sth->execute('polyphen_version', $polyphen_version);
    }

}


sub write_output {
    my $self = shift;
    
    unless ($self->param('pph_run_type') == NONE) {
        $self->dataflow_output_id($self->param('pph_output_ids'), 2);
    }

    unless ($self->param('sift_run_type') == NONE) {
        $self->dataflow_output_id($self->param('sift_output_ids'), 3);
    }
}


sub get_refseq_transcripts {
    my $self = shift;

    my $of_dba = $self->get_species_adaptor('otherfeatures');
    my $sa = $of_dba->get_SliceAdaptor or die "Failed to get slice adaptor";
    my $slices = $sa->fetch_all('toplevel', undef, 1, undef, ($self->param('include_lrg') ? 1 : undef));

    my $vep_obj;
    if(my $bam = $self->param('bam')) {
        
        die("ERROR: ensembl-vep and Bio::DB::HTS modules must be installed to add edited RefSeq transcripts\n") unless eval q{
            use Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript;
            use Bio::EnsEMBL::VEP::Config;
            use Bio::DB::HTS;
            1;
        };

        $vep_obj = Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript->new({
            config      => Bio::EnsEMBL::VEP::Config->new(),
            source_type => 'refseq',
            bam         => $bam,
            dir         => '/dev/null',
            valid_chromosomes => [map {$_->seq_region_name} @$slices],
        }) or die "Failed to get VEP AnnotationSource object";

        # we need synonyms as BAM file won't have same chr names as DB
        my $syn_file = $self->param('species_dir')."/chr_synonyms.txt" or die "Failed to write to synonyms file";
        open SYN, ">$syn_file";

        my $srsa = $of_dba->get_SeqRegionSynonymAdaptor() or die "Failed to get synonyms adaptor";

        # synonyms can be indirect i.e. A <-> B <-> C
        # and there may not be a direct link between A <-> C in the DB
        # so let's allow for one level of indirection
        my $tree = {};
        foreach my $syn(@{$srsa->fetch_all}) {
          my $syn_slice = $sa->fetch_by_seq_region_id($syn->seq_region_id);
          next unless $syn_slice;
          my ($a, $b) = sort ($syn_slice->seq_region_name, $syn->name);
          $tree->{$a}->{$b} = 1;
          $tree->{$_}->{$b} = 1 for keys %{$tree->{$a} || {}};
          $tree->{$_}->{$a} = 1 for keys %{$tree->{$b} || {}};
        }

        # now create uniq A <-> B / B <-> A pairs
        my %uniq;
        foreach my $a(keys %$tree) {
          foreach my $b(grep {$a ne $_} keys %{$tree->{$a}}) {
            $uniq{join("\t", sort ($a, $b))} = 1;
          }
        }

        print SYN "$_\n" for keys %uniq;
        close SYN;

        # now load them to the VEP obj
        $vep_obj->chromosome_synonyms($syn_file);

        $Bio::EnsEMBL::VEP::AnnotationType::Transcript::DEBUG = 1 if $self->param('debug_mode');
    }

    my @transcripts;
    
    my $ga = $of_dba->get_GeneAdaptor or die "Failed to get gene adaptor";
    
    if($self->param('debug_mode')) {    
        @transcripts =
            grep { $_->translation }
            map {$vep_obj->apply_edits($_); $_}
            @{ $ga->fetch_all_by_external_name('BRCA2')->[0]->get_all_Transcripts };
    }

    else {
        for my $slice (@{$slices}) {
            for my $gene (@{ $slice->get_all_Genes(undef, undef, 1) }) {
                for my $transcript (grep {$_->stable_id =~ /^NM_/} @{ $gene->get_all_Transcripts }) {
                    $vep_obj->apply_edits($transcript) if $vep_obj;

                    if (my $translation = $transcript->translation) {
                        push @transcripts, $transcript;
                    }
                }
            }
        }
    }

    return \@transcripts;
}

1;
