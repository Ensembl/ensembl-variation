package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::RunWeka;

use strict;

use Bio::EnsEMBL::Hive::Utils qw(stringify);
use Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix;

use base ('Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::BaseProteinFunction');

my $DEBUG = 0;

sub run {
    my $self = shift;

    my $translation_stable_id = $self->required_param('translation_stable_id');

    my $feature_file = $self->required_param('feature_file');

    my $pph_dir = $self->required_param('pph_dir');

    if ($feature_file =~ /\.gz$/ && -e $feature_file) {    
        system("gunzip -f $feature_file") == 0 or die "Failed to gunzip feature_file: $feature_file";
    }

    $feature_file =~ s/.gz$//;

    my ($output_dir) = $feature_file =~ /(.+)\/${translation_stable_id}.features$/;

    my @to_delete;

    my $output_file = "${output_dir}/${translation_stable_id}.out";
    my $error_file  = "${output_dir}/${translation_stable_id}.weka_stderr";

    push @to_delete, $feature_file, $output_file;
    
    my $cmd = "$pph_dir/bin/run_weka.pl $feature_file 1> $output_file 2> $error_file";

    if ($DEBUG) {
        $cmd = "cp $feature_file $output_file";
    }

    system($cmd) == 0 or die "Failed to run $cmd: $?";

    if (-s $error_file) {
        warn "run_weka.pl STDERR output in $error_file\n";
    }
    else {
        push @to_delete, $error_file;
    }

    open (RESULT, "<$output_file") or die "Failed to open output file: $!";

    my @fields;

    my $peptide = $self->get_transcript_file_adaptor->get_translation_seq($translation_stable_id);

    my $pred_matrix = Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix->new(
        -analysis       => 'polyphen',
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

    system("gzip -f $feature_file") == 0 or warn "Failed to gzip $feature_file: $?"; 
    
    unlink @to_delete;
}

1;
