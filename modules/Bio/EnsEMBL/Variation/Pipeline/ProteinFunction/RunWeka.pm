package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::RunWeka;

use strict;

use Bio::EnsEMBL::Hive::Utils qw(stringify);

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');

my $DEBUG = 0;

sub run {
    my $self = shift;

    my $transcript_stable_id = $self->param('transcript_stable_id');

    my $feature_file = $self->param('feature_file');

    my $pph_dir = $self->param('pph_dir'),
    
    my $var_dba = $self->get_species_adaptor('variation');

    my $dbh = $var_dba->dbc->db_handle;

    my $aa = $var_dba->get_AttributeAdaptor;

    my $save_sth = $dbh->prepare_cached(qq{
        INSERT INTO polyphen_prediction (
            protein_position_id,
            amino_acid,
            prediction_attrib_id,
            probability
        ) 
        VALUES (?,?,?,?)
    }) or die "DB error: ".$dbh->errstr;
    
    my $save_extra_sth = $dbh->prepare_cached(qq{
        INSERT INTO polyphen_supplementary_data (
            polyphen_prediction_id,
            compressed_result_hash
        ) 
        VALUES (?,COMPRESS(?))
    }) or die "DB error: ".$dbh->errstr;

    if ($feature_file =~ /\.gz$/ && -e $feature_file) {    
        system("gunzip -f $feature_file") == 0 or die "Failed to gunzip feature_file: $feature_file";
    }

    $feature_file =~ s/.gz$//;

    my ($output_dir) = $feature_file =~ /(.+)\/${transcript_stable_id}.features$/;

    my @to_delete;

    my $output_file = "${output_dir}/${transcript_stable_id}.out";
    my $error_file  = "${output_dir}/${transcript_stable_id}.weka_stderr";

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

    my $get_pos_sth = $dbh->prepare_cached(qq{
        SELECT  pp.protein_position_id
        FROM    protein_position pp, protein_info pi
        WHERE   pp.protein_info_id = pi.protein_info_id
        AND     pi.transcript_stable_id = ?
        AND     pp.position = ?
    }) or die "DB error: ".$dbh->errstr;

    open (RESULT, "<$output_file") or die "Failed to open output file: $!";

    my @fields;

    my @output_ids;

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
        
        # fetch and delete information we store in columns
        
        my $tran_ver    = delete $results{o_acc};
        my $alt_aa      = delete $results{o_aa2};
        my $prediction  = delete $results{prediction};
        my $prob        = delete $results{pph2_prob};
        my $position    = delete $results{o_pos};

        # delete results we don't need
        
        for my $val (qw{
                o_snp_id
                o_pos
                o_aa1
                snp_id
                acc
                pos
                aa1
                aa2
                nt1
                nt2
                based_on
                effect
            }) {
            delete $results{$val};
        }

        # get rid of any fields with no results

        for my $key (keys %results) {
            delete $results{$key} unless length $results{$key};
        }

        # serialize the hash (if anything remains in it)

        my $result_string = keys %results ? stringify(\%results) : undef;

        # fetch the relevant protein_position_id

        my ($transcript_stable_id_from_file, $transcript_version) = split /\./, $tran_ver;

        die "Mismatching transcript stable ids in $feature_file" 
            unless $transcript_stable_id_from_file eq $transcript_stable_id;

        $get_pos_sth->execute(
            $transcript_stable_id, 
            $position
        );

        my ($pos_id) = $get_pos_sth->fetchrow_array;

        die "No protein_position for $transcript_stable_id pos $position" unless $pos_id;
        
        # store the results in the database

        my $pred_attrib_id = $aa->attrib_id_for_type_value('polyphen_prediction', $prediction);

        die "No attrib_id for polyphen prediction: '$prediction'!" unless defined $pred_attrib_id;

        $save_sth->execute(
            $pos_id,
            $alt_aa,
            $pred_attrib_id,
            $prob
        );

        my $pph_pred_id = $dbh->last_insert_id(undef, undef, undef, undef);

        $save_extra_sth->execute(
            $pph_pred_id,
            $result_string
        );
    }

    system("gzip -f $feature_file") == 0 or warn "Failed to gzip $feature_file: $?"; 

    unlink @to_delete;
}

sub write_output {
    my $self = shift;
}

1;
