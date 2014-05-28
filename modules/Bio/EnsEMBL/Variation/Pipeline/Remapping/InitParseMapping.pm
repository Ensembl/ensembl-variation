package Bio::EnsEMBL::Variation::Pipeline::Remapping::InitParseMapping;

use strict;

use base ('Bio::EnsEMBL::Hive::Process');

use File::Path qw(make_path);

sub run {
    my $self = shift;
}

sub write_output {
    my $self = shift;

    my $bam_files_dir   = $self->param('bam_files_dir');
    my $fasta_files_dir = $self->param('fasta_files_dir');
    my $mapping_results_dir = $self->param('mapping_results_dir');
    my @input;

    if ($self->param('mode') eq 'remap_read_coverage') {
       opendir (IND_DIR, $fasta_files_dir) or die $!;
        while (my $individual_dir = readdir(IND_DIR)) {
            next if ($individual_dir =~ /^\./);
            make_path("$mapping_results_dir/$individual_dir") or die "Failed to create dir $mapping_results_dir/$individual_dir $!";;
            
            my $file_number = 1;
            opendir(my $dh, "$fasta_files_dir/$individual_dir") or die "Not a directory $fasta_files_dir/$individual_dir";
            my $max_file_count = scalar(grep {$_ =~ m/\.fa$/} readdir($dh));
            closedir($dh);

            while ($file_number <= $max_file_count) {
                my $params = {};
                die "Fasta file $fasta_files_dir/$individual_dir/$file_number.fa" unless (-e "$fasta_files_dir/$individual_dir/$file_number.fa");
                die "Bam file $bam_files_dir/$individual_dir/$file_number.bam missing" unless (-e "$bam_files_dir/$individual_dir/$file_number.bam");
                push @input, {
                    'individual_id' => $individual_dir,
                    'file_number'   => $file_number,
                    'fasta_file'    => "$fasta_files_dir/$individual_dir/$file_number.fa",
                    'bam_file'      => "$bam_files_dir/$individual_dir/$file_number.bam",
                };
                $file_number++;
            }
        }
        closedir(IND_DIR); 

    } else {

        my $file_number = 1;
        opendir(my $dh, $fasta_files_dir) or die "Not a directory $fasta_files_dir";
        my $max_file_count = scalar(grep {$_ =~ m/\.fa$/} readdir($dh));
        closedir($dh);

        while ($file_number <= $max_file_count) {
            my $params = {};
            die "Fasta file $fasta_files_dir/$file_number.fa" unless (-e "$fasta_files_dir/$file_number.fa");
            die "Bam file $bam_files_dir/$file_number.bam missing" unless (-e "$bam_files_dir/$file_number.bam");
            push @input, {
                'file_number' => $file_number,
                'fasta_file'  => "$fasta_files_dir/$file_number.fa",
                'bam_file'    => "$bam_files_dir/$file_number.bam",
            };
            $file_number++;
        }
    }
    $self->dataflow_output_id(\@input, 2);
    1;
}


1;
