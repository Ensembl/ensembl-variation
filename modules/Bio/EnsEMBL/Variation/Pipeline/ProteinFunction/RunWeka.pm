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
package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::RunWeka;

use strict;

use File::Copy;
use File::Path qw(make_path remove_tree);

use Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix;

use base ('Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::BaseProteinFunction');

sub run {
    my $self = shift;

    my $translation_md5 = $self->required_param('translation_md5');
    
    my $feature_file = $self->required_param('feature_file');

    my $pph_dir = $self->required_param('pph_dir');

    my $humdiv_model = $self->required_param('humdiv_model');
    my $humvar_model = $self->required_param('humvar_model');

    # copy stuff to /tmp to avoid lustre slowness

    my ($output_dir, $feature_filename) = $feature_file =~ /(.+)\/([^\/]+)$/;

    my $tmp_dir = "/tmp/weka_${translation_md5}";

    make_path($tmp_dir);
    
    copy($feature_file, $tmp_dir);

    my $input_file = "${tmp_dir}/${feature_filename}";

    # unzip the file if necessary

    if ($input_file =~ /\.gz$/ && -e $input_file) {    
        system("gunzip -f $input_file") == 0 or die "Failed to gunzip input file: $input_file";
    }

    $input_file =~ s/.gz$//;

    chdir $output_dir or die "Failed to chdir to $output_dir";

    my @to_delete;

    my $var_dba = $self->get_species_adaptor('variation');

    my $pfpma = $var_dba->get_ProteinFunctionPredictionMatrixAdaptor
        or die "Failed to get matrix adaptor";

    for my $model ($humdiv_model, $humvar_model) {

        next unless $model;

        my $model_name = $model eq $humdiv_model ? 'humdiv' : 'humvar';

        my $output_file = "${tmp_dir}/${model_name}.txt";

        my $error_file  = "${model_name}.err";

        my $cmd = "$pph_dir/bin/run_weka.pl -l $model $input_file 1> $output_file 2> $error_file";

        system($cmd) == 0 or die "Failed to run $cmd: $?";

        if (-s $error_file) {
            warn "run_weka.pl STDERR output in $error_file\n";
        }
        else {
            push @to_delete, $error_file;
        }

        open (RESULT, "<$output_file") or die "Failed to open output file: $!";

        my @fields;

        my $peptide = $self->get_transcript_file_adaptor->get_translation_seq($translation_md5);

        my $pred_matrix = Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix->new(
            -analysis           => 'polyphen',
            -sub_analysis       => $model_name,
            -peptide_length     => length($peptide),
            -translation_md5    => $translation_md5,
        );

        while (<RESULT>) {
            if (/^#/) {
                s/#//g;
                @fields = split /\s+/;
                next;
            }

            die "No header line in result file $output_file?" unless @fields; 

            my @values = split /\t/;

            # trim whitespace

            map { $_ =~ s/^\s+//; $_ =~ s/\s+$// } @values; 
            
            # parse the results into a hash

            my %results = map { $fields[$_] => $values[$_] } (0 .. @fields-1);
            
            my $alt_aa      = $results{o_aa2};
            my $prediction  = $results{prediction};
            my $prob        = $results{pph2_prob};
            my $position    = $results{o_pos};

            next unless $position && $alt_aa;

            $pred_matrix->add_prediction(
                $position,
                $alt_aa,
                $prediction, 
                $prob,
            );
        }

        # save the predictions to the database

        $pfpma->store($pred_matrix);
    }

    remove_tree($tmp_dir);
    
    unlink @to_delete;
}

1;
