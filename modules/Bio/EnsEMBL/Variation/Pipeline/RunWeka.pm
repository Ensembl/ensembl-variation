=head1 LICENSE

 Copyright (c) 1999-2011 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

package Bio::EnsEMBL::Hive::RunnableDB::nsSNP::RunWeka;

use strict;

use Bio::EnsEMBL::Hive::Utils qw(stringify);

use base ('Bio::EnsEMBL::Hive::RunnableDB::nsSNP::BaseRunnable');

my $DEBUG = 0;

sub _run {
    my $self = shift;

    my $gene_id = $self->param('gene_stable_id') 
        or die "gene_stable_id is a required parameter";

    my $feature_file = $self->param('feature_file')
        or die "feature_file is a required parameter";

    my $pph_dir = $self->param('pph_dir')
        or die "pph_dir is a required parameter";

    if ($feature_file =~ /\.gz$/) {    
        system("gunzip $feature_file") == 0 or die "Failed to unzip feature_file";
        $feature_file =~ s/.gz$//;
    }

    my ($output_dir) = $feature_file =~ /(.+)\/${gene_id}.features$/;

    #print "OUTPUT_DIR: $output_dir\nFEATURE_FILE: $feature_file\nGENE_ID: $gene_id\n";

    my $output_file = "${output_dir}/${gene_id}.out";
    my $error_file  = "${output_dir}/${gene_id}.weka_stderr";

    $self->delete_after($feature_file, $output_file);
    
    my $cmd = "$pph_dir/bin/run_weka.pl $feature_file 1> $output_file 2> $error_file";

    if ($DEBUG) {
        $cmd = "cp $feature_file $output_file";
    }

    system($cmd) == 0 or die "Failed to run $cmd: $?";

    if (-s $error_file) {
        warn "run_weka.pl STDERR output in $error_file\n";
    }
    else {
        $self->delete_after($error_file);
    }

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

        my %results = map { $fields[$_] => $values[$_] } (0 .. @fields-1);

        my $result_hash = {
            name        => $results{o_snp_id},
            protein     => $results{o_acc},
            position    => $results{pos},
            analysis    => 'polyphen2',
            prediction  => $results{prediction},
            result_hash => stringify(\%results),
        };

        push @output_ids, $result_hash;
    }


    $self->param('output_ids', \@output_ids);
}

sub _write_output {
    my $self = shift;
    
    if (my $output_ids = $self->param('output_ids')) {
        $self->dataflow_output_id($output_ids, 3);
    }
}

1;
