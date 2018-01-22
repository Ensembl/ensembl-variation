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


use warnings;
use strict;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp expand);
use Getopt::Long;
use Fcntl qw( LOCK_SH LOCK_EX );
use Progress;

#ÊA hard-coded hash containing the subroutines to call for each check
my %ALLELE_PREDICATE = (
    4 => \&novariation_alleles,
    13 => \&illegal_character_alleles,
    14 => \&ambiguous_alleles
);

my %SUBSNP_PREDICATE = (
);

my %VARIATION_ALLELE_PREDICATE = (
    11 => \&mismatched_allele_string,
    12 => \&multiple_alleles
);

my %VARIATION_FEATURE_PREDICATE = (
    1 => \&multiple_mappings,
    2 => \&reference_mismatch,
    3 => \&multiple_alleles,
    5 => \&no_mapping,
    13 => \&illegal_character_alleles,
    14 => \&ambiguous_alleles,
    15 => \&inconsistent_coords
);

# Accepted alleles
my @ACCEPTED_ALLELE = (
    'HGMD_MUTATION'
);

my %AMBIG_REGEXP_HASH = (
    'M' =>  '[AC]',
    'R' =>  '[AG]',
    'W' =>  '[AT]',
    'S' =>  '[CG]',
    'Y' =>  '[CT]',
    'K' =>  '[GT]',
    'V' =>  '[ACG]',
    'H' =>  '[ACT]',
    'D' =>  '[AGT]',
    'B' =>  '[CGT]',
    'X' =>  '[ACGT]',
    'N' =>  '[ACGT]'
);
#ÊGet a string containing the possible ambiguity nucleotides
my $AMBIGUITIES = join("",keys(%AMBIG_REGEXP_HASH));
# Add the code for uracil in case some allele should have that
%AMBIG_REGEXP_HASH = (%AMBIG_REGEXP_HASH,('U' =>  'T'));

# The maximum number of mappings before the variation is flagged
my $MAX_MAP_WEIGHT = 3;

# The maximum number of different alleles a variation is permitted to have 
my $MAX_ALLELES = 3;

#ÊThe option definitions
my @defs = (
    'registry_file=s',
    'qc=s@',
    'output_dir=s',
    'variation_id_range=s',
    'task_management_file=s',
    'task_id=i',
    'species=s',
    'group=s',
    'scope=s',
    'parallelize=i',
    'source_id=i@',
    'help!'
);

#ÊParse the command line and store the results in the options hash
my %options;
GetOptions(\%options,@defs);

# Check that we got a registry configuration file
die ("You need to provide a registry configuration file") unless (defined($options{'registry_file'}));
# Check that a species was specified
die ("You need to provide a species") unless (defined($options{'species'}));

#ÊIf no output dir was specified, use the current working one
my $outdir = $options{'output_dir'};
$outdir ||= "";
# Append a slash if we have a directory
if (length($outdir)) {
    $outdir .= "/";
}

#ÊLoad the registry and get a DBAdaptor to the variation database we're processing (or the group specified on the command line)
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all($options{'registry_file'});
my $species = $options{'species'};
my $group = $options{'group'};
$group ||= 'variation';
my $dba = $registry->get_DBAdaptor($species,$group) or die ("Could not get a DBAdaptor for $species - $group");

#ÊIf the option to parallelize was specified, we will chunk the task into the desired sizes and create the corresponding task management file
if ($options{'parallelize'}) {
    
    # Check that a desired task_management_file was specified
    die ("You must specify a file where the task parameters will be written") unless (defined($options{'task_management_file'}));
    
    my $chunksize = $options{'parallelize'};
    
    # Get the min and max variation_ids and simply assume that the data is evenly distributed on average w.r.t. variation_id
    my $stmt = qq{
        SELECT
            MIN(variation_id),
            MAX(variation_id)
        FROM
            variation
    };
    my ($min_id,$max_id) = @{$dba->dbc->db_handle->selectall_arrayref($stmt)->[0]};
    
    # Divide the id range into chunks and write to management file
    open (TASK,">",$options{'task_management_file'}) or die ("Could not open " . $options{'task_management_file'} . " for writing");
    my $offset = $min_id;
    my $task_id = 0;
    while ($offset <= $max_id) {
        $task_id++;
        print TASK join("\t",($task_id,$offset,($offset+$chunksize-1))) . "\n";
        $offset += $chunksize;
    }
    close(TASK);
    
    print STDOUT "The task has been divided into chunks of $chunksize. The parameters have been written to " . $options{'task_management_file'} . ". You should submit this as a job array over the indexes 1-$task_id\n";
    exit(0);
} 

# We will probably need a core dbadaptor as well so create one
my $dba_core = $registry->get_DBAdaptor($species,'core') or warn ("Could not get a DBAdaptor for $species - core");
 
#ÊGet the range of variations we should work on. This can either be specified by:
# 1. A variation_id range specified on the command line
# 2. Provided in a task management file specified on the command line. This overrides a specified range. 
#    If this is the case then a job index corresponding to a row in the task management file must be specified.
#    This can either be done on the command line or through the LSB_JOBINDEX environment variable (which gets set by LSF in a jobarray submission).
#    The latter overrides the former.
# 3. None of the above, in which case all variations will be processed
my ($lower_id,$upper_id);
if (defined($options{'task_management_file'})) {
    
    my $job_index = $ENV{'LSB_JOBINDEX'};
    $job_index ||= $options{'task_id'};
    
    # Check that we have a job index
    die ("A task management file was specified but not a task index, can not proceed") unless (defined($job_index));
    
    # Get the variation_id range for this job index
    open(TASK,"<",$options{'task_management_file'}) or die ("Could not open task management file " . $options{'task_management_file'} . " for parsing");
    while (<TASK>) {
        chomp;
        my @arr = split(/\s+/,$_); 
        ($lower_id,$upper_id) = ($arr[1],$arr[2]) if ($arr[0] == $job_index);
    }
    close(TASK);
    
    # Check that we could find the range
    die ("Could not find the corresponding variation_id range for task index $job_index") unless (defined($lower_id) && defined($upper_id));
    
    # Print the job assignment to STDERR
    print STDERR "Job $job_index works on range $lower_id - $upper_id ";
} 
#ÊElse, we check for a comma-separated range
elsif (defined($options{'variation_id_range'})) {
    ($lower_id,$upper_id) = split(",",$options{'variation_id_range'});
}

my $failed_variation_file = $outdir . "failed_variation.txt";
my $failed_allele_file = $outdir . "failed_allele.txt";
my $loadfile = {
    'variation' => $failed_variation_file,
    'allele' => $failed_allele_file
};

### Now, get the data from the database

# Get the haplotype seq region ids
our $HAPLOTYPE_IDS = get_haplotype_seq_region_ids($dba_core);
# Get the failed description ids
my %failed_description = %{get_failed_description($dba,$options{'qc'})};
my @failed_description_ids = keys(%failed_description);


# A hash to hold the variation_ids and the tests that it failed
my %failed_variation;

# A hash to hold the allele_ids and the tests that it failed
my %failed_allele;
        
#ÊCheck if we should do the checking for variations
my $scope = lc($options{'scope'});
$scope ||= 'variation';
if ($scope eq 'variation') {
    
    #ÊLoop over the variation features and flag them as appropriate
    
    #ÊIf a variation_id range was specified, create a condition on it
    my $condition = get_range_condition($lower_id,$upper_id,"v");
    
    # If a source_id condition was specified, append this to the condition
    $condition .= " AND " . get_source_condition($options{'source_id'},"v");
    
    my $stmt = qq{
        SELECT
            v.variation_id,
            v.name,
            vf.variation_feature_id,
            vf.seq_region_id,
            vf.seq_region_start,
            vf.seq_region_end,
            vf.seq_region_strand,
            vf.allele_string,
            ras.ref_allele,
            ra.seq_region_strand,
            'variation'
        FROM
            variation v LEFT JOIN 
            variation_feature vf ON (
                vf.variation_id = v.variation_id
            ) LEFT JOIN 
            (
                tmp_ref_allele ra JOIN
                tmp_ref_allele_seq ras ON (
                    ras.ref_allele_seq_id = ra.ref_allele_seq_id
                )
            ) ON (
                ra.variation_feature_id = vf.variation_feature_id
            )
        WHERE
            $condition
        ORDER BY
            v.variation_id; 
    };
    my $sth = $dba->dbc->prepare($stmt);
    
    # Execute the query
    $sth->execute();
    
    # Loop over the variation features
    my @vf_arr;
    my @row = $sth->fetchrow_array();
    while (@row) {
        
        # Add the row to the array grouping the same variation_ids into an array
        push(@vf_arr,[@row]);
        
        # Get the next row
        my @nextrow = $sth->fetchrow_array();
        
        #ÊIf we are switching variation or we have no more rows, do the checks
        if (!scalar(@nextrow) || $nextrow[0] != $row[0]) {
            
            #ÊExecute the predicates
            if (scalar(@vf_arr)) {
                my @failed;
                # Cache the results in a hash
                my $cache = {};
                map {
                    push(@failed,$_) if (exists($VARIATION_FEATURE_PREDICATE{$_}) && $VARIATION_FEATURE_PREDICATE{$_}->(\@vf_arr,$cache));
                } @failed_description_ids;
                
                $failed_variation{$row[0]} = \@failed if (scalar(@failed));
            }            
            
            # Empty the variation array
            splice(@vf_arr);
        }
        
        @row = @nextrow;
    }
}


if ($scope eq 'allele') {
    
    #ÊLoop over the variation features and flag them as appropriate
    
    #ÊIf a variation_id range was specified, create a condition on it
    my $condition = get_range_condition($lower_id,$upper_id,"a");
    
    my $stmt = qq{
        SELECT
            a.allele_id,
            a.subsnp_id,
            a.variation_id,
            vf.seq_region_id,
            vf.allele_string,
            vf.seq_region_end,
            vf.seq_region_strand,
            a.allele,
            NULL,
            NULL,
            'allele'
        FROM
            allele a LEFT JOIN 
            variation_feature vf ON (
                vf.variation_id = a.variation_id
            )
        WHERE
            $condition
        ORDER BY
            a.variation_id,
            a.subsnp_id; 
    };
    my $sth = $dba->dbc->prepare($stmt);

    # Execute the query
    $sth->execute();
    
    # Loop over the joined rows. We'll send off checks both for individual alleles, subsnps and variations
    my @variation;
    my @subsnp;
    my @allele;
    my @row = $sth->fetchrow_array();
    while (@row) {
        
        # Variation array
        push(@variation,[@row]);
        push(@subsnp,[@row]);
        push(@allele,[@row]);
        
        # Get the next row
        my @nextrow = $sth->fetchrow_array();
        
        #ÊIf we are switching allele or we have no more rows, do the checks for alleles
        if (!scalar(@nextrow) || $nextrow[0] != $row[0]) {
            
            #ÊExecute the predicates
            if (scalar(@allele)) {
                my @failed;
                # Cache the results in a hash
                my $cache = {};
                map {
                    push(@failed,$_) if (exists($ALLELE_PREDICATE{$_}) && $ALLELE_PREDICATE{$_}->(\@allele,$cache));
                } @failed_description_ids;
                
                if (scalar(@failed)) {
                    map {$failed_allele{$_->[0]} = \@failed} @allele;
                }
            }            
            
            # Empty the array
            splice(@allele);
        }
        
        #ÊIf we are switching subsnp or we have no more rows, do the checks for subsnp
        if (!scalar(@nextrow) || $nextrow[1] != $row[1]) {
            
            #ÊExecute the predicates
            if (scalar(@subsnp)) {
                my @failed;
                # Cache the results in a hash
                my $cache = {};
                map {
                    push(@failed,$_) if (exists($SUBSNP_PREDICATE{$_}) && $SUBSNP_PREDICATE{$_}->(\@subsnp,$cache));
                } @failed_description_ids;
                
                if (scalar(@failed)) {
                    map {$failed_allele{$_->[0]} = \@failed} @subsnp;
                }
            }            
            
            # Empty the array
            splice(@subsnp);
        }
        
        #ÊIf we are switching variation or we have no more rows, do the checks for variations
        if (!scalar(@nextrow) || $nextrow[2] != $row[2]) {

            #ÊExecute the predicates
            if (scalar(@variation)) {
                my @failed;
                # Cache the results in a hash
                my $cache = {};
                map {
                    push(@failed,$_) if (exists($VARIATION_ALLELE_PREDICATE{$_}) && $VARIATION_ALLELE_PREDICATE{$_}->(\@variation,$cache));
                } @failed_description_ids;
                
                if (scalar(@failed)) {
                    $failed_variation{$row[2]} = \@failed;
                }
            }            
            
            # Empty the variation feature array
            splice(@variation);
        }
        
        @row = @nextrow;
    }

    
}

foreach my $scope (('variation','allele')) {
    
    my %h;
    if ($scope eq 'variation') {
        %h = %failed_variation;
    }
    else {
        %h = %failed_allele;
    }
    
    # Only dump to file if we have any results
    next unless (scalar(keys(%h)));
    
    # Open the loadfile (append) and get a lock on it 
    open(LOAD,">>",$loadfile->{$scope}) or die ("Could not open loadfile " . $loadfile->{$scope} . " for writing");
    flock(LOAD,LOCK_EX);
    
    #ÊWrite the ids and the failed_description_id to the load file
    foreach my $id (keys(%h)) {
        map {print LOAD "$id\t$_\n"} @{$h{$id}};
    }
    
    close(LOAD);
}

#ÊIf we finished successfully, print that to STDERR
print STDERR " Finished ok!\n";

#ÊCheck if a variation is mapped to more than the maximum allowed number of (non-haplotype) genomic locations 
sub multiple_mappings {
    my $variation_features = shift;
    my $cache = shift;
    
    # If the result of this test has been cached return it
    my $failed_description_id = 1;
    unless (exists($cache->{$failed_description_id})) {
        $cache->{$failed_description_id} = _multiple_mappings($variation_features,$cache);
    }
    return $cache->{$failed_description_id};
}
sub _multiple_mappings {
    my $variation_features = shift;
    my $cache = shift;
    
    my $count = 0;
    foreach my $vf (@{$variation_features}) {
        
        next unless (defined($vf->[3]));
        next if (grep {$vf->[3] == $_} @{$HAPLOTYPE_IDS});
        
        $count++;
        
        return 1 if ($count > $MAX_MAP_WEIGHT);
    }
    
    return 0;
}

#ÊCheck if the allele string provided by dbSNP is in agreement with the alleles of all subsnps belonging to the variation 
sub mismatched_allele_string {
    my $rows = shift;
    my $cache = shift;
    
    # If the result of this test has been cached return it
    my $failed_description_id = 11;
    unless (exists($cache->{$failed_description_id})) {
        $cache->{$failed_description_id} = _mismatched_allele_string($rows,$cache);
    }
    return $cache->{$failed_description_id};
}
sub _mismatched_allele_string {
    my $rows = shift;
    my $cache = shift;
    
    # If this variation has no mapping, it won't have any allele string associated
    return 0 if (no_mapping($rows,$cache));
    
    # Get the unique alleles from the subsnps
    my %ss = map {$_->[7] => 1} @{$rows};
    
    # Get the unique alleles from the variation feature allele string
    my %vf = map {map {$_ => 1} split(/\//,$_->[4])} @{$rows};
    
    # Check that all subsnp alleles are present in the allele_string
    map {return 1 unless (exists($vf{$_}))} keys(%ss);
    
    # Check that all allele_string alleles are present in the subsnp alleles
    map {return 1 unless (exists($ss{$_}))} keys(%vf);
    
    return 0;
}

# Check if a variation has no mappings
sub no_mapping {
    my $rows = shift;
    my $cache = shift;
    
    # If the result of this test has been cached return it
    my $failed_description_id = 5;
    unless (exists($cache->{$failed_description_id})) {
        $cache->{$failed_description_id} = _no_mapping($rows,$cache);
    }
    return $cache->{$failed_description_id};
}
sub _no_mapping {
    my $rows = shift;
    my $cache = shift;
    
    return (defined($rows->[0][3]) ? 0 : 1);
}    

# Check if the coordinates given for a variation is not compatible with its allele string
sub inconsistent_coords {
    my $rows = shift;
    my $cache = shift;
    
    # If the result of this test has been cached return it
    my $failed_description_id = 15;
    unless (exists($cache->{$failed_description_id})) {
        $cache->{$failed_description_id} = _inconsistent_coords($rows,$cache);
    }
    return $cache->{$failed_description_id};
}
sub _inconsistent_coords {
    my $rows = shift;
    my $cache = shift;
    
    # If this variation has no mappings, it shouldn't be classified as inconsistent
    return 0 if (no_mapping($rows,$cache));
    # If this variation contains illegal characters, there's no point in checking for inconsistent coordinates
    return 0 if (illegal_character_alleles($rows,$cache));
    
    #ÊThe only things we accept is if the position is a deletion or if at least one of the alleles are of the same length as the position
    foreach my $variation_feature (@{$rows}) {
        
        expand(\$variation_feature->[7]);
        my $ref_len = ($variation_feature->[5] - $variation_feature->[4] + 1);
        
        #ÊMatching lengths or deletion and insertion in allele string?
        next if (grep {($_ eq '-' && $ref_len == 0) || (length($_) == $ref_len)} split(/\//,$variation_feature->[7]));
        
        # Else, this is inconsistent coordinates
        return 1;
    }
    
    return 0;
}

#ÊCheck if the allele string alleles does not agree with the reference sequence
sub reference_mismatch {
    my $rows = shift;
    my $cache = shift;
    
    # If the result of this test has been cached return it
    my $failed_description_id = 2;
    unless (exists($cache->{$failed_description_id})) {
        $cache->{$failed_description_id} = _reference_mismatch($rows,$cache);
    }
    return $cache->{$failed_description_id};
}
sub _reference_mismatch {
    my $rows = shift;
    my $cache = shift;
    
    # If this variation has no mappings, it shouldn't be classified as a mismatch
    return 0 if (no_mapping($rows,$cache));
    
    # Get the unique reference alleles
    my $ref_allele = _unique_reference_allele($rows);
    
    # Get the unique allele strings
    my $allele_string = _unique_allele_string($rows);
    
    # Loop over the allele strings and match them to the reference alleles
    foreach my $as (@{$allele_string}) {
        
        expand(\$as);
        map {
            my $allele = $_;
            return 0 if (grep {mismatch($allele,$_) == 0} @{$ref_allele});
        } split(/\//,$as);
        
    }
    
    # Nothing matched
    return 1;
}

# Check if a sequence (possibly) ambiguous mismatches another
sub mismatch {
    my $allele = shift;
    my $reference = shift;
     
    # If they match
    return 0 if ($allele eq $reference);
    #ÊReturn mismatch if allele doesn't contains ambig codes
    return 1 unless (ambiguous(\$allele));
    
    # Turn the sequence into regexps if necessary
    ambiguity_to_regexp(\$allele);
    
    # By now, the allele should only contain nucleotide characters and brackets. 
    
    # Do a regexp matching
    return 0 if ($reference =~ m/^$allele$/);
    return 1;
}

# Check if the allele string contains too many single nucleotide alleles
sub multiple_alleles {
    my $rows = shift;
    my $cache = shift;
    
    # If the result of this test has been cached return it
    my $failed_description_id = ($rows->[0][10] eq 'variation' ? 3 : 12);
    unless (exists($cache->{$failed_description_id})) {
        $cache->{$failed_description_id} = _multiple_alleles($rows,$cache);
    }
    return $cache->{$failed_description_id};
}
sub _multiple_alleles {
    my $rows = shift;
    my $cache = shift;
    
    # If this variation has no mappings, it won't have any allele strings
    #return 0 if (no_mapping($rows,$cache) && $rows->[0][10] eq 'variation');
    
    # Get the unique allele strings
    my $allele_string = _unique_allele_string($rows);
    
    foreach my $a_string (@{$allele_string}) {
        expand(\$a_string);
        my $count = grep {$_ =~ m/^[ACGT]$/i} split(/\//,$a_string);
        return 1 if ($count > $MAX_ALLELES);
    }
    
    return 0; 
}

#ÊCheck if a variation's allele strings contain ambiguity codes
sub ambiguous_alleles {
    my $rows = shift;
    my $cache = shift;
    
    # If the result of this test has been cached return it
    my $failed_description_id = 14;
    unless (exists($cache->{$failed_description_id})) {
        $cache->{$failed_description_id} = _ambiguous_alleles($rows,$cache);
    }
    return $cache->{$failed_description_id};
}
sub _ambiguous_alleles {
    my $rows = shift;
    my $cache = shift;
    
    my @alleles;
    #ÊCheck if we are dealing with a variation feature or alleles
    if ($rows->[0][10] eq 'variation') {
    
        # If this variation has no mappings, it won't have any illegal characters in the allele_string
        #return 0 if (no_mapping($rows,$cache));
    
        # Get the unique allele strings
        my $allele_string = _unique_allele_string($rows);
    
        map {push(@alleles,split(/\//,$_))} @{$allele_string};
    }
    else {
        push(@alleles,$rows->[0][7]);
    }
    
    foreach my $allele (@alleles) {
            
        #ÊExpand the allele
        expand(\$allele);
        #ÊReport the allele if it contains 'illegal' characters
        return 1 if (ambiguous(\$allele));
        
    }
    
    return 0;
}

# Check if an allele contains ambiguity codes, but make sure that it doesn't contain 'illegal' characters
sub ambiguous {
    my $allele_ref = shift;
    return (${$allele_ref} =~ m/[$AMBIGUITIES]/i && !illegal_characters($allele_ref));
}

#ÊCheck if a variation's allele strings contain illegal characters
sub illegal_character_alleles {
    my $rows = shift;
    my $cache = shift;
    
    # If the result of this test has been cached return it
    my $failed_description_id = 13;
    unless (exists($cache->{$failed_description_id})) {
        $cache->{$failed_description_id} = _illegal_character_alleles($rows,$cache);
    }
    return $cache->{$failed_description_id};
}
sub _illegal_character_alleles {
    my $rows = shift;
    my $cache = shift;
    
    my @alleles;
    #ÊCheck if we are dealing with a variation feature or alleles
    if ($rows->[0][10] eq 'variation') {
    
        # If this variation has no mappings, it won't have any illegal characters in the allele_string
        #return 0 if (no_mapping($rows,$cache));
    
        # Get the unique allele strings
        my $allele_string = _unique_allele_string($rows);
    
        map {push(@alleles,split(/\//,$_))} @{$allele_string};
    }
    else {
        push(@alleles,$rows->[0][7]);
    }
    
    foreach my $allele (@alleles) {
            
        #ÊExpand the allele
        expand(\$allele);
        #ÊReport the allele if it contains 'illegal' characters
        return 1 if (illegal_characters(\$allele));
        
    }
    
    return 0;
}

#ÊCheck if an allele is a 'NOVARIATION'
sub novariation_alleles {
    my $rows = shift;
    my $cache = shift;
    
    # If the result of this test has been cached return it
    my $failed_description_id = 4;
    unless (exists($cache->{$failed_description_id})) {
        $cache->{$failed_description_id} = _novariation_alleles($rows,$cache);
    }
    return $cache->{$failed_description_id};
}
sub _novariation_alleles {
    my $rows = shift;
    my $cache = shift;
    
    return 1 if (grep {novariation(\$_->[7])} @{$rows});
    return 0;
}

#ÊKeep a list of accepted alleles that won't be flagged as containing illegal characters. Check if an allele is in this list
sub accepted {
    my $allele_ref = shift;
    map {return 1 if ($_ eq ${$allele_ref})} @ACCEPTED_ALLELE;
    return 0;
}

# Check if an allele is 'NOVARIATION'
sub novariation {
    my $allele_ref = shift;
    return (${$allele_ref} eq 'NOVARIATION');
}

# Check if an allele contains 'illegal' characters
sub illegal_characters {
    my $allele_ref = shift;
    return (${$allele_ref} =~ m/[^ACGTU\-$AMBIGUITIES]/i && !accepted($allele_ref) && !novariation($allele_ref));
}

# Replace ambiguity codes in a sequence with a suitable regular expression
sub ambiguity_to_regexp {
    my $seq_ref = shift;
    
    ${$seq_ref} =~ s/([U$AMBIGUITIES])/$AMBIG_REGEXP_HASH{$1}/ig;
};

#ÊPrivate method to get the unique allele strings from variation features
sub _unique_allele_string {
    my $variation_features = shift;
    
    # Check first if this is just a single row
    return [$variation_features->[0][7]] if (scalar(@{$variation_features}) == 1);
    
    # Get the unique allele strings
    my %allele_string;
    map {
        $allele_string{$_->[7]}++;
    } @{$variation_features};
    
    my @unique = keys(%allele_string);
    
    # If it is alleles rather than a variation we're looking at, create an allele string from the alleles
    if ($variation_features->[0][10] eq 'allele') {
        my $as = join("/",@unique);
        @unique = ($as);
    }
    return \@unique;
}

#ÊPrivate method to get the reference alleles from variation features
sub _unique_reference_allele {
    my $variation_features = shift;
    
    # Check first if this is just a single row
    if (scalar(@{$variation_features}) == 1) {
        # Flip the reference allele if necessary
        reverse_comp($variation_features->[0][8]) unless ($variation_features->[0][9] == $variation_features->[0][6]);
        return [$variation_features->[0][8]];
    }
    
    # Get the unique reference alleles
    my %ref_allele;
    map {
        # Flip the reference allele if necessary
        reverse_comp(\$_->[8]) unless ($_->[9] == $_->[6]);
        $ref_allele{$_->[8]}++;
    } @{$variation_features};
    
    my @unique = keys(%ref_allele);
    return \@unique;
}

sub get_haplotype_seq_region_ids {
    my $dba_core = shift;
    
    #ÊThe haplotype regions have attribs 'non reference'. So do the LRGs however, so filter by name to exclude these
    my $stmt = qq{
        SELECT
            sr.seq_region_id
        FROM
            seq_region sr JOIN
            seq_region_attrib sra ON (
                sra.seq_region_id = sr.seq_region_id
            ) JOIN 
            attrib_type at ON (
                at.attrib_type_id = sra.attrib_type_id
            )
        WHERE
            sr.name NOT LIKE 'lrg%' AND
            at.name LIKE 'non reference'
    };
    my $haplotype_ids = $dba_core->dbc->db_handle->selectcol_arrayref($stmt);
    
    return $haplotype_ids;
}

sub get_range_condition {
    my $lower_id = shift;
    my $upper_id = shift;
    my $alias = shift;
    
    return " 1 " unless (defined($lower_id) && defined($upper_id));
    
    return (defined($alias) ? " $alias\." : " ") . qq{variation_id BETWEEN $lower_id AND $upper_id };
}

sub get_source_condition {
    my $ids = shift;
    my $alias = shift;
    
    return " 1 " unless (defined($ids) && scalar(@{$ids}));
    
    my $condition = " (" . (defined($alias) ? "$alias\." : "") . "source_id = " . join(" OR " . (defined($alias) ? "$alias\." : "") . "source_id = ",@{$ids}) . ") ";
    return $condition;
}

sub get_failed_description {
    my $dba = shift;
    my $ids = shift;
    
    my $condition = " 1 ";
    if (defined($ids) && scalar(@{$ids})) {
        $condition = " failed_description_id IN (" . join(",",@{$ids}) . ") ";
    }
    
    my $stmt = qq{
        SELECT
            failed_description_id,
            description
        FROM
            failed_description
        WHERE
            $condition
    };
    #ÊGet a hashref of the descriptions with the failed_description_id as key
    my $description = $dba->dbc->db_handle->selectall_hashref($stmt,'failed_description_id');
    
    return $description;
}



=head
#ÊLoop over the failed_description_ids and for each, call the corresponding subroutine. Each check will return a hashref with arrayrefs of failed variation_ids and allele_ids, respectively and we write these to the corresponding dump file.
foreach my $failed_description_id (keys(%{$failed_description})) {
    
    # Print some progress information to stdout
    print STDOUT Progress::location() . "\tFlagging variations/alleles for '" . $failed_description->{$failed_description_id}{'description'} . "' (failed_description_id = $failed_description_id)\n";  
    
    # Warn and skip if we don't know how to perform this check
    unless (exists($PREDICATE{$failed_description_id})) {
        warn ("Can not determine the corresponding subroutine to use for consistency check '" . $failed_description->{$failed_description_id}{'description'} . "' (failed_description_id = $failed_description_id). Skipping");
        next;
    }
    
    # Call the checking subroutine
    my $routine = $PREDICATE{$failed_description_id};
    my $flagged = $routine->($dba,$lower_id,$upper_id,$dba_core);
    
    # Loop over the flagged variations and alleles and write them to the dump files
    foreach my $type (('variation','allele')) {
        
        #ÊGet the ids that were returned
        my $ids = $flagged->{$type} || [];
        
        #ÊIf no ids were flagged, skip
        next unless (scalar(@{$ids}));
        
        # Print some progress information to stdout
        print STDOUT Progress::location() . "\tDumping flagged " . $type . "s to loadfile\n";
    
        # Open the loadfile (append) and get a lock on it 
        open(LOAD,">>",$loadfile->{$type}) or die ("Could not open loadfile " . $loadfile->{$type} . " for writing");
        flock(LOAD,LOCK_EX);
        
        #ÊWrite the ids and the failed_description_id to the load file
        while (my $id = shift(@{$ids})) {
            print LOAD join("\t",($id,$failed_description_id)) . "\n";
        }
        
        close(LOAD);
    }
}

sub get_haplotype_condition {
    my $dba_core = shift;
    
    my $haplotype_seq_region_ids = get_haplotype_seq_region_ids($dba_core);
    return " 1 " unless (defined($haplotype_seq_region_ids) && scalar(@{$haplotype_seq_region_ids}));
    return " seq_region_id NOT IN (" . join(",",@{$haplotype_seq_region_ids}) . ") ";
}

#ÊCheck if a variation is mapped to more than the maximum allowed number of (non-haplotype) genomic locations 
sub multiple_mappings {
    my $dba = shift;
    my $lower_id = shift;
    my $upper_id = shift;
    my $dba_core = shift;
    
    #ÊIf a range was specified, create a condition on it
    my $condition = get_range_condition($lower_id,$upper_id);
    $condition .= " AND " . get_haplotype_condition($dba_core);

    my $stmt = qq{
        SELECT 
            variation_id
        FROM
            variation_feature
        WHERE
            $condition 
    };
    
    # Add the group and condition on maximum mappings
    $stmt .= qq{ 
        GROUP BY
            variation_id
        HAVING
            COUNT(*) > $MAX_MAP_WEIGHT
    };
    
    # Execute the query and get the result
    my $flagged_variation_ids = $dba->dbc->db_handle->selectcol_arrayref($stmt);
    
    # Return a hashref with the result
    return {'variation' => $flagged_variation_ids};
}

#ÊCheck whether the variation has at least one allele that matches the reference
sub reference_mismatch {
    my $dba = shift;
    my $lower_id = shift;
    my $upper_id = shift;
    
    #ÊIf a range was specified, create a condition on it
    my $condition = get_range_condition($lower_id,$upper_id);
    
    # Statement to get the variation alleles
    my $stmt = qq{
        SELECT
            allele_id,
            subsnp_id,
            variation_id,
            allele
        FROM
            allele
        WHERE
            $condition
        ORDER BY
            variation_id
    };
    my $sth = $dba->dbc->db_handle($stmt);
    
    #ÊStatement to get the reference sequence for each variation_feature
    $stmt = qq{
        SELECT
            ras.ref_allele,
            ra.seq_region_strand
        FROM
            variation_feature vf JOIN
            tmp_ref_allele ra ON (
                ra.variation_feature_id = vf.variation_feature_id
            ) JOIN
            tmp_ref_allele_seq ON (
                ras.ref_allele_seq_id = ra.ref_allele_seq_id
            )
        WHERE
            vf.variation_id = ?
    };
    my $seq_sth = $dba->dbc->prepare($stmt);
    
    # Get the alleles
    $sth->execute();
    my ($allele_id,$subsnp_id,$variation_id,$allele,$refseq,$refstrand,$last_variation_id);
    $sth->bind_columns(\$allele_id,\$subsnp_id,\$variation_id,\$allele);
    
    $last_variation_id = -1;
    while ($sth->fetch()) {
        
        #ÊIf we switched variation, get the possible reference sequences
        if ($variation_id != $last_variation_id) {
            
            $seq_sth->execute($variation_id);
            $seq_sth->bind_columns(\$refseq,\$refstrand);
            
        }
        
    }
    
}

#ÊCheck that a variation does not have more than the maximum allowed number of single-nucleotide alleles (based on subsnps) 
sub multiple_alleles {
    my $dba = shift;
    my $lower_id = shift;
    my $upper_id = shift;
    
    #ÊIf a range was specified, create a condition on it
    my $condition = get_range_condition($lower_id,$upper_id);
    
    #ÊStatement to get the alleles for the variation_id range
    my $stmt = qq{
        SELECT
            variation_id,
            allele 
        FROM
            allele
        WHERE
            $condition
        ORDER BY
            variation_id
    };
    my $sth = $dba->dbc->prepare($stmt);
    
    # Execute the statement and bind the result columns
    $sth->execute();
    my ($variation_id,$allele,$last_variation_id,$last_flagged);
    $sth->bind_columns(\$variation_id,\$allele);
    my %alleles;
    $last_variation_id = -1;
    
    # An array to hold the variation_id for flagged variations
    my @flagged;
    
    #ÊLoop over the alleles
    while ($sth->fetch()) {
        
        #ÊReset the allele hash and the flagged status if we are moving to a new variation_id
        if ($variation_id != $last_variation_id) {
            %alleles = ();
            $last_flagged = 0;
            $last_variation_id = $variation_id;
        }
        
        # Skip if we have already flagged this variation
        next if ($last_flagged);
        
        # If this is a single bp allele and it's not a deletion, add it to the hash
        if (length($allele) == 1 && $allele ne '-') {
            $alleles{$allele}++;
            
            # Check the size of the hash and flag the variation if it is greater than the maximum number of allowed alleles
            if (scalar(keys(%alleles)) > $MAX_ALLELES) {
                push(@flagged,$variation_id);
                $last_flagged = 1;
            }
        }
    }
    
    # Return the flagged variations
    return {'variation' => \@flagged}; 
}

# Check that the variation has a mapping to the genome
sub no_mapping {
    my $dba = shift;
    my $lower_id = shift;
    my $upper_id = shift;
    
    #ÊIf a range was specified, create a condition on it
    my $condition = get_range_condition($lower_id,$upper_id,"v");
    
    #ÊStatement to check for unmapped variations
    my $stmt = qq{
        SELECT
            v.variation_id
        FROM
            variation v LEFT JOIN
            variation_feature vf ON (
                vf.variation_id = v.variation_id
            )
        WHERE
            $condition AND
            vf.variation_feature_id IS NULL
    };
    # Execute the query and get the result
    my $flagged_variation_ids = $dba->dbc->db_handle->selectcol_arrayref($stmt);
    
    return {'variation' => $flagged_variation_ids};
}

#ÊCheck if this is a 'NoVariation'
sub no_variation {
    my $dba = shift;
    my $lower_id = shift;
    my $upper_id = shift;
    
    #ÊIf a range was specified, create a condition on it
    my $condition = get_range_condition($lower_id,$upper_id);
    
    #ÊStatement to get the alleles that are 'NoVariation'
    my $stmt = qq{
        SELECT
            allele_id,
            variation_id
        FROM
            allele
        WHERE
            $condition AND
            allele LIKE 'novariation'
    };
    
    return _check_allele_variation($dba,$stmt);
}

# Check that there are no disallowed (e.g. 'N') alleles
sub disallowed_alleles {
    my $dba = shift;
    my $lower_id = shift;
    my $upper_id = shift;
    
    #ÊIf a range was specified, create a condition on it
    my $condition = get_range_condition($lower_id,$upper_id);
    
    # Define a number of regexps for things that we do allow and catch the rest
    my $normal_variation = '^[-ACGT]+$';
    my $microsatellite = '^\\([ACGTMRWSYKVHDBXN]+\\)[0-9]+';
    my $novariation = '^NOVARIATION$';
    my $hgmd = '^HGMD_MUTATION$';
    
    #ÊStatement to catch non-accepted alleles
    my $stmt = qq{
        SELECT
            allele_id,
            variation_id
        FROM
            allele
        WHERE
            $condition AND
            allele NOT REGEXP '$normal_variation' AND
            allele NOT REGEXP '$microsatellite' AND
            allele NOT REGEXP '$novariation' AND
            allele NOT REGEXP '$hgmd'
    };
    
    return _check_allele_variation($dba,$stmt);
}  
  
# 'internal' function that checks alleles and whether all alleles for the corresponding variation have failed
sub _check_allele_variation {
    my $dba = shift;
    my $stmt = shift;
    
    my $sth = $dba->dbc->prepare($stmt);
    $sth->execute();
    my ($allele_id,$variation_id);
    $sth->bind_columns(\$allele_id,\$variation_id);
    my %variation_ids;
    my @flagged_alleles;
    my @flagged_variations;
    
    # Loop over the alleles and flag them. At the same time, count the number of alleles for each variation_id that has this allele string
    while ($sth->fetch()) {
        push(@flagged_alleles,$allele_id);
        $variation_ids{$variation_id}++;
    }
    
    # In order to determine if the variation should be flagged as well as the allele, count the number of alleles for each variation and see if it corresponds to the number of failed alleles
    $stmt = qq{
        SELECT
            COUNT(*) 
        FROM
            allele
        WHERE
            variation_id = ?
    };
    $sth = $dba->dbc->prepare($stmt);
    
    # Loop over the variaiton_ids concerned
    while (my ($variation_id,$count) = each(%variation_ids)) {
        
        $sth->execute($variation_id);
        
        # If the count matches the number of alleles, we should flag the variation as well 
        if ($count == $sth->fetchrow_arrayref()->[0]) {
            push(@flagged_variations,$variation_id);
        } 
    }
    
    # Return the flagged variations and alleles
    return {'variation' => \@flagged_variations, 'allele' => \@flagged_alleles};
}
=cut
  
