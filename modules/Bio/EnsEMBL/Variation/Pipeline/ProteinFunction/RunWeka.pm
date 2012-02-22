package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::RunWeka;

use strict;

use Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix;

use base ('Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::BaseProteinFunction');

sub run {
    my $self = shift;

    my $translation_md5 = $self->required_param('translation_md5');
    
    my $feature_file = $self->required_param('feature_file');

    my $pph_dir = $self->required_param('pph_dir');

    my $humdiv_model = $self->required_param('humdiv_model');
    my $humvar_model = $self->required_param('humvar_model');

    if ($feature_file =~ /\.gz$/ && -e $feature_file) {    
        system("gunzip -f $feature_file") == 0 or die "Failed to gunzip feature_file: $feature_file";
    }

    $feature_file =~ s/.gz$//;

    my ($output_dir, $feature_filename) = $feature_file =~ /(.+)\/([^\/]+)$/;

    chdir $output_dir or die "Failed to chdir to $output_dir";

    my @to_delete;

    for my $model ($humdiv_model, $humvar_model) {

        next unless $model;

        my $model_name = $model eq $humdiv_model ? 'humdiv' : 'humvar';

        my $output_file = "${model_name}.txt";

        push @to_delete, $output_file;

        my $error_file  = "${model_name}.err";

        my $cmd = "$pph_dir/bin/run_weka.pl -l $model $feature_filename 1> $output_file 2> $error_file";

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
            -analysis       => 'polyphen',
            -sub_analysis   => $model_name,
            -peptide_length => length($peptide),
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

        $self->save_predictions($pred_matrix);
    }

    system("gzip -f $feature_file") == 0 or warn "Failed to gzip $feature_file: $?"; 
    
    unlink @to_delete;
}

1;
