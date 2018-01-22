#!/usr/bin/env perl
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



=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut


#####
# Script for dumping out the genomic reference sequence for each variation 
#

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use FindBin qw( $Bin );
use Getopt::Long;

# Define the options to the script
my %options = ();
my @option_defs = (
  'species=s',
  'output_dir=s',
  'seq_region=s',
  'parallelize=i',
  'registry_file=s',
  'joblimit=i'
);

# Parse the command line
GetOptions(\%options,@option_defs);

#ÊGet the options or die/set default options
my $registryfile = $options{'registry_file'} || die("A registry file must be specified");
my $species = $options{'species'} || die("A species name must be specified");
my $output_dir = $options{'output_dir'} || '.';
my $joblimit = $options{'joblimit'};
$joblimit ||= 10;

#ÊThe seq_region_id can be explicitly specified, or it will correspond to the LSB_JOBINDEX variable
my $seq_region_id = $options{'seq_region'};
$seq_region_id ||= $ENV{'LSB_JOBINDEX'}; 

#ÊLoad the registry from the configuration file
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all($registryfile);

my $dba = $registry->get_DBAdaptor("human","variation");
my $stmt;

#ÊIf we are asked to parallelize, submit the parallel jobs, one for each seq_region in the variation database, and then exit when they have finished
if ($options{'parallelize'}) {
    
    # Unset the parallelize option so we won't pass it to the subjobs
    delete($options{'parallelize'});
    
    # Get all the seq_region_ids from the variation database
    $stmt = qq{
        SELECT DISTINCT
            seq_region_id
        FROM
            variation_feature
        ORDER BY
            seq_region_id
    };
    my @seq_region_ids = map {$_->[0]} @{$dba->dbc->db_handle->selectall_arrayref($stmt)};
    
    # A comma-separated list of seq_region_ids to use as job array elements
    my $jobarray = join(",",@seq_region_ids);
    
    # The bsub command
    my $bsub_cmd = "bsub -J ref_seq_export[$jobarray]\%$joblimit -o $output_dir/ref_seq_export.\%J.\%I.out -e $output_dir/ref_seq_export.\%J.\%I.err perl $Bin/export_variation_reference_sequence.pl";
    #ÊAdd the parameters
    while (my ($key,$value) = each(%options)) {
        $bsub_cmd .= " -$key $value";
    }
    print $bsub_cmd;
    # Submit to the farm 
    system($bsub_cmd);
    exit(0);
}

# First, create the necessary tables if they do not already exist
$stmt = qq{
    CREATE TABLE IF NOT EXISTS
        tmp_ref_allele_seq (
          ref_allele_seq_id INTEGER NOT NULL AUTO_INCREMENT,
          ref_allele VARCHAR(25000) NOT NULL,
          PRIMARY KEY (ref_allele_seq_id)
        )
};
$dba->dbc->do($stmt);

$stmt = qq{
    CREATE TABLE IF NOT EXISTS 
        tmp_ref_allele (
          variation_feature_id INTEGER NOT NULL,
          ref_allele_seq_id INTEGER NOT NULL,
          seq_region_strand TINYINT NOT NULL,
          PRIMARY KEY (variation_feature_id)
        )
};
$dba->dbc->do($stmt);

# Write the data to a temp file
my $tempfile = "$output_dir/tmp_refseq_data_$seq_region_id\.txt";
open(DATA,">",$tempfile) or die ("Could not open tempfile $tempfile for writing");

# Fetch the seq_region slice from the core database
my $slice_adaptor = $registry->get_adaptor($species,"core","slice");
my $slice = $slice_adaptor->fetch_by_seq_region_id($seq_region_id) or die ("Could not fetch slice for seq_region_id $seq_region_id");

# A hash to hold allele -> allele_id mapping
my %allele_id;

# Prepare a statement for looking up seqs
$stmt = qq{
    SELECT
        ref_allele_seq_id
    FROM
        tmp_ref_allele_seq
    WHERE
        ref_allele = ?
};
my $fetch_sth = $dba->dbc->prepare($stmt);

# Insert statement
$stmt = qq{
    INSERT INTO
        tmp_ref_allele_seq (
            ref_allele 
        )
    VALUES (
        ?
    )
};
my $insert_sth = $dba->dbc->prepare($stmt);

# Get the variation_feature_id and position on the seq_region for all features on the seq_region_id
$stmt = qq{
    SELECT
        variation_feature_id,
        seq_region_start,
        seq_region_end,
        seq_region_strand
    FROM
        variation_feature FORCE INDEX (pos_idx)
    WHERE
        seq_region_id = ?
};
my $sth = $dba->dbc->prepare($stmt);
$sth->execute($seq_region_id);
my ($vf_id,$start,$end,$strand,$seq,$id);
$sth->bind_columns(\$vf_id,\$start,\$end,\$strand);

#ÊLoop over the rows, print the data to the tempfile
while ($sth->fetch()) {
    
    # If an insertion, the allele should be '-'
    if ($start <= $end) {
        $seq = $slice->subseq($start,$end,$strand);
    }
    else {
        $seq = '-';
    }
    
    # Get the allele id for the sequence
    $id = $allele_id{$seq};
    #ÊIf not defined, look it up (or insert) in the database
    unless (defined($id)) {
        $id = lookup_id($seq,$fetch_sth,$insert_sth,$dba);
        $allele_id{$seq} = $id;
    }
    
    # Push data onto array
    print DATA join("\t",($vf_id,$id,$strand)) . "\n";
}

close(DATA);

# Load the data into the table, disable keys first
$stmt = qq{
    ALTER TABLE
        tmp_ref_allele
    DISABLE KEYS
};
$dba->dbc->do($stmt);

$stmt = qq{
    LOAD DATA LOCAL INFILE
        '$tempfile'
    INTO TABLE
        tmp_ref_allele (
            variation_feature_id,
            ref_allele_seq_id,
            seq_region_strand
        )
};
$dba->dbc->do($stmt);

# unlink($tempfile);

sub lookup_id {
    my $seq = shift;
    my $fetch_sth = shift;
    my $insert_sth = shift;
    my $dba = shift;
    
    my $id;
    $fetch_sth->execute($seq);
    $fetch_sth->bind_columns(\$id);
    $fetch_sth->fetch();
    
    # If id could not be found, insert it
    unless (defined($id)) {
        $insert_sth->execute($seq);
        $id = $dba->dbc->db_handle->selectall_arrayref(qq{SELECT LAST_INSERT_ID()})->[0][0]; 
    }
    
    return $id;
}

 


 
