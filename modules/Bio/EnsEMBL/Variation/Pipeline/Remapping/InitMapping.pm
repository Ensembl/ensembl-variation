#!/usr/bin/env perl
# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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




=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut
package Bio::EnsEMBL::Variation::Pipeline::Remapping::InitMapping;

use strict;
use warnings;

use FileHandle;
use Bio::DB::Fasta;
use Bio::EnsEMBL::Registry;

use base ('Bio::EnsEMBL::Hive::Process');

sub fetch_input {
    my $self = shift;
    my $registry = 'Bio::EnsEMBL::Registry';
    $registry->load_all($self->param('registry_file'));
    my $cdba = $registry->get_DBAdaptor($self->param('species'), 'core');
    my $vdba = $registry->get_DBAdaptor($self->param('species'), 'variation');

    $self->param('registry', $registry);
    $self->param('cdba', $cdba);
    $self->param('vdba', $vdba);
}

sub run {
    my $self = shift;
    if ($self->param('mode') eq 'remap_multi_map') {
        $self->dump_multi_map_features();
    } else {
        if ($self->param('generate_fasta_files')) {
            $self->dump_features();
            $self->generate_mapping_input();
        }
    }
}

sub write_output {
    my $self = shift;
    # initialise mapping jobs
    my $file_count      = $self->param('file_count');
    my $fasta_files_dir = $self->param('fasta_files_dir');
    my $bam_files_dir   = $self->param('bam_files_dir');
    my @jobs;
    my $i = 1;
    while ($i <= $file_count) {
        push @jobs, {
            'file_number'   => $i,
            'bam_files_dir' => $bam_files_dir,
            'fasta_file'    => "$fasta_files_dir/$i.fa",
            'sam_file'      => "$bam_files_dir/$i.sam",
            'bam_file'      => "$bam_files_dir/$i.bam",
            'err_file'      => "$bam_files_dir/$i.err",
            'out_file'      => "$bam_files_dir/$i.out",
        };
        $i++;
    }
    $self->dataflow_output_id(\@jobs, 2);
}

sub generate_mapping_input {
    my $self = shift;

    my $old_assembly_fasta_file_dir = $self->param('old_assembly_fasta_file_dir');
    my $fasta_db = Bio::DB::Fasta->new($old_assembly_fasta_file_dir, -reindex => 1);
    $self->param('fasta_db', $fasta_db);

    # store end-coordinates for all seq_regions to check that variation_location + flank_seq_length < slice_end
    my $cdba = $self->param('cdba');
    my $sa = $cdba->get_SliceAdaptor;
    # don't include asm exceptions, fetch the full length of the Y chromosome
    my $slices = $sa->fetch_all('toplevel', undef, 0, 1);
    my $seq_regions = {};
    foreach my $slice (@$slices) {
        my $end             = $slice->end;
        my $seq_region_name = $slice->seq_region_name;
        if ($end > $seq_regions->{$seq_region_name}->{end}) {
            $seq_regions->{$seq_region_name}->{end} = $end;
        }
    }
    foreach my $seq_region_name (keys %$seq_regions) {
        my $end = $seq_regions->{$seq_region_name}->{end};
    }

    $self->param('seq_regions', $seq_regions);

    my $dump_features_dir = $self->param('dump_features_dir');
    my $fasta_files_dir   = $self->param('fasta_files_dir');
    my $pipeline_dir      = $self->param('pipeline_dir');
    my $fh_allele_length = FileHandle->new("$pipeline_dir/qc_allele_string_length.txt", 'w');
    my $fh_ref_seq       = FileHandle->new("$pipeline_dir/qc_ref_seq.txt", 'w');

    $self->param('qc_allele_string_length', $fh_allele_length);
    $self->param('qc_ref_seq', $fh_ref_seq);

    my $variants_with_multi_map = {};

    opendir(DIR, $dump_features_dir) or die $!;
    while (my $file = readdir(DIR)) {
        if ($file =~ /^(.+)\.txt$/) {
            my $file_number = $1;
            my $fh = FileHandle->new("$dump_features_dir/$file", 'r');
            my $fh_fasta_file = FileHandle->new("$fasta_files_dir/$file_number.fa", 'w');
            while (<$fh>) {
                chomp;
                my $data = $self->read_line($_);
                my $seq_region_name = $data->{seq_region_name},
                my $feature_id      = $data->{variation_feature_id};
                my $start           = $data->{seq_region_start};
                my $end             = $data->{seq_region_end};
                my $strand          = $data->{seq_region_strand};
                my $allele_string   = $data->{allele_string};
                my $map_weight      = $data->{map_weight};
                my $variation_name  = $data->{variation_name};
                if ($map_weight > 1) {
                    $variants_with_multi_map->{$variation_name}++;
                    next;
                }
                my ($flank_start, $upstream_flank_length, $downstream_flank_length, $flank_end, $variant_length) = @{$self->flank_coordinates($seq_region_name, $start, $end, $strand)};

                # variant surrounded by flank sequences
                my $query_sequence = $self->get_query_sequence($seq_region_name, $flank_start, $flank_end, $strand);

                # replace empty space with underscores <_>
                $allele_string =~ s/\s/_/g;
                $self->qc_alleles($query_sequence, $upstream_flank_length, $variant_length, $allele_string, $variation_name);
#                $allele_string = $self->qc_alleles($config, $query_sequence, $upstream_flank_length, $variant_length, $allele_string, $variation_name);

                my $name = "$seq_region_name:$start:$end:$strand:$variation_name";
                my $id = ">$feature_id-$upstream_flank_length-$variant_length-$downstream_flank_length-$name";
                print $fh_fasta_file "$id\n$query_sequence\n";
            } # end while (read feature file for seq_region)
            $fh->close();
            $fh_fasta_file->close();
        }
    }

    foreach my $file_type (qw/qc_allele_string_length qc_ref_seq/) {
        my $fh = $self->param($file_type);
        $fh->close();
    }

    # store multi map
    my $fh_qc_multi_map = FileHandle->new("$pipeline_dir/multi_map.txt", 'w');

    foreach my $variation_name (keys %$variants_with_multi_map) {
        my $count = $variants_with_multi_map->{$variation_name};
        print $fh_qc_multi_map "$variation_name\t$count\n";
    }
    $fh_qc_multi_map->close();
}

sub flank_coordinates {
    my ($self, $seq_name, $start, $end, $strand) = @_;

    my $flank_length   = $self->param('flank_seq_length');
    my $flank_start    = 0;
    my $flank_end      = 0;
    my $variant_length = 0;

    my $seq_regions    = $self->param('seq_regions');
    my $slice_end      = $seq_regions->{$seq_name}->{end};

    my $add_to_end = 0;
    my $add_to_start = 0;

    my $upstream_flank_length   = $flank_length;
    my $downstream_flank_length = $flank_length;

    if ($start < $flank_length) {
        $flank_start             = 1;
        $flank_end               = $end + $flank_length + ($flank_length - $start + 1);
        $upstream_flank_length   = $start - 1;
        $downstream_flank_length = $flank_length + ($flank_length - $start + 1);
    } elsif (($end + $flank_length) > $slice_end) {
        my $add_to_start = ($end + $flank_length) - $slice_end;
        $flank_start             = $start - $flank_length - $add_to_start;
        $flank_end               = $slice_end;
        $upstream_flank_length   = $flank_length + $add_to_start;
        $downstream_flank_length = $flank_length - $add_to_start;
    } else {
        $flank_start = $start - $flank_length;
        $flank_end   = $end + $flank_length;
    }

    $variant_length = $end - $start + 1;

    # if insertion
    if ($start > $end) {
        $variant_length = 0;
    }

    return [$flank_start, $upstream_flank_length, $downstream_flank_length, $flank_end, $variant_length];
}

sub get_query_sequence {
    my ($self, $seq_name, $start, $end, $strand) = @_;
    my $fasta_db = $self->param('fasta_db');
    if ($strand == -1) {
        return $fasta_db->seq("$seq_name:$end,$start");
    }
    return $fasta_db->seq("$seq_name:$start,$end");
}

sub qc_alleles {
    my ($self, $query_sequence, $upstream_flank_length, $variant_length, $allele_string, $variation_name) = @_;

    my $upstream   = substr $query_sequence, 0, $upstream_flank_length;
    my $ref        = substr $query_sequence, $upstream_flank_length, $variant_length;
    my $downstream = substr $query_sequence, $upstream_flank_length + $variant_length;

    my @alleles = split('/', $allele_string);
    my $ref_allele = $alleles[0];
    if ($ref_allele eq '-') {
        $ref_allele = '';
    }
    if ($ref ne $ref_allele) {
        my $fh = $self->param('qc_ref_seq');
        print $fh "$variation_name\t$ref_allele\t$ref\n";
    }
}

sub dump_multi_map_features {
    my $self = shift;

    # parse fasta files and store: rsname to fasta_file_number
    # Query_name in fasta_file >14086.1-110-1-497-G:C/G:rs708635
    # parse variation_features on ref_sequence: store variation_feature info in dump_feature/file_number
    # make a note if feature was stored

    my $fasta_files_dir = $self->param('fasta_files_dir');
    my $rs_id_to_file_number = {};
    my $file_numbers = {};
    opendir(DIR, $fasta_files_dir) or die $!;
    while (my $file = readdir(DIR)) {
        if ($file =~ /^(.+)\.fa$/) {
            my $file_number = $1;
            $file_numbers->{$file_number} = 1;
            my $fh = FileHandle->new("$fasta_files_dir/$file", 'r');
            while (<$fh>) {
                chomp;
                my $line = $_;
                if ($line =~ /^>/) {
                    my @query_name_components = split(':', $line);
                    my $rs_id = pop @query_name_components;
                    $rs_id_to_file_number->{$rs_id} = $file_number;
                }
            }
        }
    }
    closedir(DIR);

    my $dump_features_dir = $self->param('dump_features_dir');
    my $fh_hash = {};
    foreach my $file_number (keys %$file_numbers) {
        my $fh = FileHandle->new("$dump_features_dir/$file_number.txt", 'w');
        $fh_hash->{$file_number} = $fh;
    }

    my $file_count = scalar keys %$file_numbers;
    $self->param('file_count', $file_count);

    my $cdba = $self->param('cdba');
    my $vdba = $self->param('vdba');

    my $dbname = $vdba->dbc->dbname();
    my $feature_table = $self->param('feature_table');

    my $dbh = $vdba->dbc->db_handle;
    my $sth = $dbh->prepare(qq{
        SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS
        WHERE TABLE_SCHEMA = '$dbname'
        AND TABLE_NAME = '$feature_table';
    });
    $sth->execute();

    # QC that all necessary columns are there: e.g. seq_region_id, ...
    my @column_names = ();
    while (my @name = $sth->fetchrow_array) {
        push @column_names, $name[0];
    }
    $sth->finish();
    @column_names = sort @column_names;
    my $column_names_concat = join(',', @column_names);
    $self->param('sorted_column_names', $column_names_concat);
    $sth->finish();

    my $sa = $cdba->get_SliceAdaptor;
    # don't include asm exceptions
    my $slices = $sa->fetch_all('toplevel', undef, 0, 1);

    my $seq_region_ids = {};
    foreach my $slice (@$slices) {
        my $seq_region_name = $slice->seq_region_name;
        next if ($seq_region_name =~ /PATCH/);
        my $seq_region_id = $slice->get_seq_region_id;
        $seq_region_ids->{$seq_region_id} = $seq_region_name;
    }

    $sth = $dbh->prepare(qq{
        SELECT variation_name,map_weight,$column_names_concat FROM $feature_table WHERE seq_region_id = ?;
    }, {mysql_use_result => 1});

    my $vf_info_is_stored = {};

    foreach my $seq_region_id (keys %$seq_region_ids) {
        my $seq_region_name = $seq_region_ids->{$seq_region_id};
        $sth->execute($seq_region_id);
        while (my $row = $sth->fetchrow_arrayref) {
            my @values = map { defined $_ ? $_ : '\N' } @$row;
            my $variation_name = shift @values;
            my $map_weight     = shift @values;
            next if ($map_weight == 1);
            next if ($vf_info_is_stored->{$variation_name});
            my $file_number = $rs_id_to_file_number->{$variation_name};
            unless ($file_number) {
                $self->warning("No sequence information for $variation_name from dbSNP");
                next;
            }
            my $fh = $fh_hash->{$file_number};
            my @pairs = ();
            for my $i (0..$#column_names) {
                push @pairs, "$column_names[$i]=$values[$i]";
            }
            push @pairs, "seq_region_name=$seq_region_name";
            print $fh join("\t", @pairs), "\n";
        }
        $sth->finish();
    }

    foreach my $file_number (keys %$fh_hash) {
        my $fh = $fh_hash->{$file_number};
        $fh->close();
    }

}

sub dump_features {
    my $self = shift;
    my $cdba = $self->param('cdba');
    my $vdba = $self->param('vdba');

    my $dbname = $vdba->dbc->dbname();
    my $feature_table = $self->param('feature_table');

    my $dbh = $vdba->dbc->db_handle;
    my $sth = $dbh->prepare(qq{
        SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS
        WHERE TABLE_SCHEMA = '$dbname'
        AND TABLE_NAME = '$feature_table';
    });
    $sth->execute();

    # QC that all necessary columns are there: e.g. seq_region_id, ...
    my @column_names = ();
    while (my @name = $sth->fetchrow_array) {
        push @column_names, $name[0];
    }
    $sth->finish();
    @column_names = sort @column_names;
    my $column_names_concat = join(',', @column_names);
    $self->param('sorted_column_names', $column_names_concat);
    $sth->finish();

    my $dump_features_dir = $self->param('dump_features_dir');
    my $sa = $cdba->get_SliceAdaptor;
    # don't include asm exceptions
    my $slices = $sa->fetch_all('toplevel', undef, 0, 1);

    my $seq_region_ids = {};
    foreach my $slice (@$slices) {
        my $seq_region_name = $slice->seq_region_name;
        next if ($seq_region_name =~ /PATCH/);
        my $seq_region_id = $slice->get_seq_region_id;
        $seq_region_ids->{$seq_region_id} = $seq_region_name;
    }

    $sth = $dbh->prepare(qq{
        SELECT $column_names_concat FROM $feature_table WHERE seq_region_id = ?;
    }, {mysql_use_result => 1});

    my $file_count = 1;
    my $entries_per_file = $self->param('entries_per_file');
    my $count_entries = 0;
    my $fh = FileHandle->new("$dump_features_dir/$file_count.txt", 'w');

   foreach my $seq_region_id (keys %$seq_region_ids) {
        my $seq_region_name = $seq_region_ids->{$seq_region_id};
        $sth->execute($seq_region_id);
        while (my $row = $sth->fetchrow_arrayref) {
            my @values = map { defined $_ ? $_ : '\N' } @$row;
            my @pairs = ();
            for my $i (0..$#column_names) {
                push @pairs, "$column_names[$i]=$values[$i]";
            }
            push @pairs, "seq_region_name=$seq_region_name";
            if ($count_entries >= $entries_per_file) {
                $fh->close();
                $file_count++;
                $fh = FileHandle->new("$dump_features_dir/$file_count.txt", 'w');
                $count_entries = 0;
            }
            $count_entries++;
            print $fh join("\t", @pairs), "\n";
        }
        $sth->finish();
    }
    $fh->close();
    $self->param('file_count', $file_count);
}

sub read_line {
    my $self = shift;
    my $line = shift;
    my @key_values = split("\t", $line);
    my $mapping = {};
    foreach my $key_value (@key_values) {
        my ($table_name, $value) = split('=', $key_value, 2);
        $mapping->{$table_name} = $value;
    }
    return $mapping;
}

1;
