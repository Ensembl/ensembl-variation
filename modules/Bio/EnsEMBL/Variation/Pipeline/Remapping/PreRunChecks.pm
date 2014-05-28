package Bio::EnsEMBL::Variation::Pipeline::Remapping::PreRunChecks;

use strict;
use warnings;

use File::Path qw(make_path remove_tree);
use FileHandle;
use Bio::EnsEMBL::Registry;

use base ('Bio::EnsEMBL::Hive::Process');


sub fetch_input {
    my $self = shift;
    1;	
}

sub run {
    my $self = shift;

    # 1. check all folders are created 
    my $pipeline_dir = $self->param('pipeline_dir');
    die "$pipeline_dir doesn't exist" unless (-d $pipeline_dir);		

    foreach my $folder (qw/bam_files_dir filtered_mappings_dir load_features_dir mapping_results_dir statistics_dir/) {
        my $dir = $self->param($folder);
        if (-d $dir) {
            remove_tree($dir);
        }
        make_path($dir);
    } 
    unless ($self->param('use_fasta_files')) {
        foreach my $folder (qw/dump_features_dir fasta_files_dir/) {
            my $dir = $self->param($folder);
            if (-d $dir) {
                remove_tree($dir);
            }
            make_path($dir);
        } 
    } else {
        my $dir = $self->param('fasta_files_dir');
        my $count = $self->count_files($dir, '.fa');
        if ($count == 0) {
            die ("There are no fasta_files. Set parameter 'generate_fasta_files' to 1 in the conf file.");
        }				
        # remove index files
        if ($self->param('mode') eq 'remap_read_coverage') {
            opendir (IND_DIR, $dir) or die $!;
            while (my $individual_dir = readdir(IND_DIR)) {
                next if ($individual_dir =~ /^\./);
                $self->run_cmd("rm $dir/$individual_dir/*.fai");
            }
            closedir (IND_DIR);
        } else {
            $self->run_cmd("rm $dir/*.fai");
        }
    }

    foreach my $folder (qw/old_assembly_fasta_file_dir new_assembly_fasta_file_dir/) {
        my $dir = $self->param($folder);
        die "$dir for $folder doesn't exist" unless (-d $dir);		
    }

    # 2. check new_assembly_fasta_file is indexed
    my $dir = $self->param('new_assembly_fasta_file_dir');
    foreach my $file_type (('.fa.amb', '.fa.ann', '.fa.bwt', '.fa.pac', '.fa.sa')) {
        unless ($self->count_files($dir, $file_type)) {
            die("New assembly file is not indexed. $file_type is missing.");
        }
    }
    # 3. check bwa and samtools are working

    1;
}

sub is_empty {
    my $self = shift;
    my $dir = shift;
    opendir(my $dh, $dir) or die "Not a directory $dir";
    my $count =  scalar(grep { $_ ne "." && $_ ne ".." } readdir($dh)) == 0;
    closedir($dh);
    return $count;
}

sub count_files {
    my $self = shift;
    my $dir = shift;
    my $file_type = shift;
    opendir(my $dh, $dir) or die "Not a directory $dir";
    my $count = scalar(grep { $_ =~ m/\Q$file_type$/ } readdir($dh)) == 0;
    closedir($dh);
    return $count;
}

sub run_cmd {
    my $self = shift;
    my $cmd = shift;
    if (my $return_value = system($cmd)) {
        $return_value >>= 8;
        die "system($cmd) failed: $return_value";
    }
}

sub write_output {
    my $self = shift;
    1;
}



1;
