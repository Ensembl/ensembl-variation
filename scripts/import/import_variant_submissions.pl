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


use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Sequence qw( reverse_comp );
use Getopt::Long;
use DBI qw(:sql_types);

our $DEFAULT_FLANKING_SIZE = 400;

my $submission_report;
my $registry_file;
my $verbose;
my $help;
my $source;
my $species;

usage() if (!scalar(@ARGV));
 
GetOptions(
    'submission_report=s' => \$submission_report,
    'registry_file=s' => \$registry_file,
    'source=s' => \$source,
    'species=s' => \$species,
    'verbose!' => \$verbose,
    'help!' => \$help
);

# Use human by default
$species ||= 'Human';

usage() if ($help);

die ("The source name of the submission is required") unless (defined($source));
die ("A submission report file is required") unless (defined($submission_report));
die ("A registry configuration file is required") unless (defined($registry_file));

# Load the registry and get a variation feature adaptor
Bio::EnsEMBL::Registry->load_all($registry_file);
my $vf_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species,'Variation','VariationFeature');
# Get a db_handle to the variation database
my $dbh = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'Variation')->dbc->db_handle();

# Check that the source exists in the source table, otherwise add it
import_source($source,$dbh);

# Open the submission report for parsing
open(FH,'<',$submission_report) or die ("Could not open $submission_report for reading");
while (<FH>) {
    chomp;
    
    # Skip comments
    next if (m/^#/);
    
    # Split the line
    my ($local_id,$hgvs,$ssid,$rsid,$condition,$omim) = split(/\t/);
    
    # Skip if local id or hgvs notation could not be found
    next unless (defined($local_id) && defined($hgvs));
    
    # Set the ssID to be undefined if it's just an empty string
    $ssid = undef unless (defined($ssid) && length($ssid));
    
    # If the rsId is not all digits (after stripping an initial rs), warn and skip to the next record
    $rsid =~ s/^rs//;
    warn ("Can not import variant with rsId '$rsid', skipping row beginning with '" . substr($_,0,35) . "...'") unless ($rsid =~ m/^\d+$/);
    
    # Get a VariationFeature from the HGVS notation
    my $vf = $vf_adaptor->fetch_by_hgvs_notation($hgvs);
    
    # Warn if we didn't get a variation feature back
    if (!defined($vf)) {
        warn ("Could not parse $hgvs into a variation feature");
        next;
    }
    
    # Set the name, source and synonym of the variation
    $vf->variation->source($source);
    
    # In case an rsId was supplied, use that
    if (defined($rsid)) {
        $vf->variation->name("rs$rsid");
        $vf->variation->source("dbSNP");
        $vf->variation->add_synonym($source,$local_id);
    }
    # Else, use the ssID if necessary
    elsif (defined($ssid)) {
        $vf->variation->name("ss$ssid");
        $vf->variation->add_synonym($source,$local_id);
    }
    # Else, use the local id
    else {
        $vf->variation->name($local_id);
    }
    $vf->source($vf->variation->source);
    
    # If defined, add the subsnp_id to the alleles
    if (defined($ssid)) {
        map {$_->subsnp($ssid)} @{$vf->variation->get_all_Alleles()};
    }
    
    # Set the validation_status of the variation_feature and variation to 'precious'
    $vf->add_validation_state('precious');
    $vf->variation->add_validation_state('precious');
    
    # Transform the variation feature to the chromosome coordinate system
    my $chr_vf = $vf->transform('chromosome');
    
    # Warn if the transform was unsuccessful but store the variation feature using LRG coordinates
    warn ("Could not transform $hgvs to the chromosome coordinate system") unless (defined($chr_vf));
    $chr_vf ||= $vf;
    
    # Check if there already exists a variation feature in the same location, in which case we should only add synonyms to this one
    my $slice = $chr_vf->feature_Slice();
    my $existing_vfs = $vf_adaptor->fetch_all_by_Slice($slice);
    
    # If we already have a variation here and it has the same rsID as the one we are importing, add synonyms to the submitter and add the alleles if necessary 
    if (scalar(@{$existing_vfs}) && grep {$_->variation_name() eq $chr_vf->variation->name()} @{$existing_vfs}) {
        
        my ($ext_vf) = grep {$_->variation_name() eq $chr_vf->variation->name()} @{$existing_vfs};
        
        # Set the variation_id of our vf to match the one in the database
        $chr_vf->dbID($ext_vf->dbID());
        $chr_vf->variation->dbID($ext_vf->variation->dbID());
        
        # Flip our vf if it's on a different strand than the existing vf
        if ($chr_vf->strand() != $ext_vf->strand()) {
          flip_variation_feature(\$chr_vf);
        }
        
        # Add the alleles of our variation
        import_alleles($chr_vf->variation,$dbh);
        
        # Add the synonyms
        import_variation_synonyms($chr_vf->variation,$dbh);
        
        # Add 'precious' to the validation_state of the existing variation if it's not set
        unless (grep {$_ =~ m/precious/} @{$ext_vf->get_all_validation_states()}) {
            
            my $vid = $ext_vf->variation->dbID();
            my $stmt = qq{
                UPDATE
                    variation_feature
                SET 
                    validation_status = CONCAT(validation_status,',precious')
                WHERE
                    variation_id = $vid
            };
            $dbh->do($stmt);
            
        }
        unless (grep {$_ =~ m/precious/} @{$ext_vf->variation->get_all_validation_states()}) {
            
            my $vid = $ext_vf->variation->dbID();
            my $stmt = qq{
                UPDATE
                    variation
                SET 
                    validation_status = CONCAT(validation_status,',precious')
                WHERE
                    variation_id = $vid
            };
            $dbh->do($stmt);
            
        }
        
    }
    # Else, add everything to the database
    else {
        
        # If the variation feature is on the negative strand on the chromsome, flip everything around to be on the positive strand
        if ($chr_vf->strand() < 0) {
            flip_variation_feature(\$chr_vf);
        }
    
        # 1) import the variation object
        import_variation($chr_vf->variation,$dbh);
        # 2) import the alleles
        import_alleles($chr_vf->variation,$dbh);
        # 3) import the flanking sequence
        import_flanking_sequence($chr_vf,$dbh);
        # 4) import the variation_synonyms
        import_variation_synonyms($chr_vf->variation,$dbh);
        # 5) import the variation_feature
        import_variation_feature($chr_vf,$dbh);
    }
    
    # Verify that the vf does not fail
    # check_failed($chr_vf,$dbh);
=head    
    # If the input data contains phenotypes, attach a variation annotation object to the variation
    if (defined($condition) && $condition ne 'NULL' && defined($omim) && $omim ne 'NULL') {
        
        my $study = Bio::EnsEMBL::Variation::Study-new(
            -name => 
        );
        my $v_annotation = Bio::EnsEMBL::Variation::VariationAnnotation->new(
            '-source_name' => $source,
            '-variation' => $chr_vf->variation(),
            '-phenotype_description' => $condition,
            '-study' => 'MIM:' . $omim,
            '-adaptor' => $vf_adaptor->db->get_VariationAnnotationAdaptor()
        );
        
        import_variation_annotation($v_annotation,$dbh);
    }
=cut    
}
close(FH);

sub flip_variation_feature {
  my $vf = shift;
  
  # Flip the allele string
  my @alleles = split(/\//,${$vf}->allele_string());
  map {reverse_comp(\$_)} @alleles;
  ${$vf}->allele_string(join("/",@alleles));
  
  # Flip the alleles
  foreach my $allele (@{${$vf}->variation()->get_all_Alleles()}) {
      my $seq = $allele->allele();
      reverse_comp(\$seq);
      $allele->allele($seq);
  }
  
  # Change the strand
  ${$vf}->strand(-1 * ${$vf}->strand());
}

sub import_variation_annotation {
    my $va = shift;
    my $dbh = shift;
    
    my $stmt = qq{
        SELECT
            va.variation_annotation_id
        FROM
            variation_annotation va JOIN
            phenotype p ON (
                p.phenotype_id = va.phenotype_id
            ) JOIN
            source s ON (
                s.source_id = va.source_id
            )
        WHERE
            va.variation_id = ? AND
            p.description = ? AND
            s.name = ?
        LIMIT 1
    };
    my $check_sth = $dbh->prepare($stmt);
    
    $stmt = qq{
        INSERT IGNORE INTO
            variation_annotation (
                variation_id,
                phenotype_id,
                source_id,
                study
            )
        VALUES (
            ?,
            ?,
            (
                SELECT
                    source_id
                FROM
                    source
                WHERE
                    name = ?
                LIMIT 1
            ),
            ?
        )
    };
    my $sth = $dbh->prepare($stmt);
    
    $check_sth->execute(
        $va->variation->dbID(),
        $va->phenotype_description(),
        $va->source_name()
    );
    my $row = $check_sth->fetchrow_arrayref();
    my $dbid;
    $dbid = $row->[0] if (defined($row));
    $va->dbID($dbid);
    
    unless (defined($va->dbID())) {
        
        my $phenotype_id = import_phenotype($va->phenotype_description(),$dbh);
        
        $sth->execute(
            $va->variation->dbID(),
            $phenotype_id,
            $va->source_name(),
            $va->study()
        );
        
        $va->dbID($dbh->last_insert_id(undef,undef,undef,undef));
    }
    
    return $va->dbID();
}

sub import_variation_synonyms {
    my $variation = shift;
    my $dbh = shift;
    
    my @dbIDs;
    
    my $check_stmt = qq{
        SELECT
            vs.variation_synonym_id
        FROM
            variation_synonym vs JOIN
            source s ON (
                s.source_id = vs.source_id
            )
        WHERE
            vs.variation_id = %d AND
            vs.subsnp_id %s AND
            s.name = '%s' AND
            vs.name = '%s'
        LIMIT 1
    };
    
    my $stmt = qq{
        INSERT IGNORE INTO
            variation_synonym (
                variation_id,
                subsnp_id,
                source_id,
                name
            )
        VALUES (
            ?,?,
            (
                SELECT
                    source_id
                FROM
                    source
                WHERE
                    name = ?
                LIMIT 1
            ),?
        )
    };
    my $sth = $dbh->prepare($stmt);
    
    # If there is a subsnp_id attached to an allele object on the variation, get it
    my $subsnp_id = $variation->get_all_Alleles->[0]->subsnp() if (scalar(@{$variation->get_all_Alleles()}));
    
    # Loop over all synonym sources
    foreach my $source (@{$variation->get_all_synonym_sources()}) {
        
        # Loop over all synonyms from this source
        foreach my $synonym (@{$variation->get_all_synonyms($source)}) {
            
            # Check if the synonym exists
            $stmt = sprintf($check_stmt,$variation->dbID(),(defined($subsnp_id) ? "= $subsnp_id" : "IS NULL"),$source,$synonym);
            my $dbid = $dbh->selectall_arrayref($stmt)->[0][0];
            
            unless (defined($dbid)) {
            
                $sth->execute(
                    $variation->dbID(),
                    $subsnp_id,
                    $source,
                    $synonym
                );
                
                $dbid = $dbh->last_insert_id(undef,undef,undef,undef);
            }
            
            push(@dbIDs,$dbid);
        }
    }
    
    return \@dbIDs;
}

sub import_variation_feature {
    my $vf = shift;
    my $dbh = shift;
    
    # Check if we already have a variation feature with this variation and position and in that case, return that dbId and do nothing
    my $stmt = qq{
        SELECT
            vf.variation_feature_id
        FROM
            variation_feature vf JOIN
            seq_region sr ON (
                vf.seq_region_id = sr.seq_region_id
            )
        WHERE
            vf.variation_id = ? AND
            sr.name = ? AND
            vf.seq_region_start = ? AND
            vf.seq_region_end = ? AND
            vf.seq_region_strand = ?
        LIMIT 1
    };
    my $sth = $dbh->prepare($stmt);
    
    $stmt = qq{
        INSERT IGNORE INTO
            variation_feature (
                seq_region_id,
                seq_region_start,
                seq_region_end,
                seq_region_strand,
                variation_id,
                allele_string,
                variation_name,
                map_weight,
                source_id,
                validation_status
            )
        VALUES (
            (
                SELECT
                    seq_region_id
                FROM
                    seq_region
                WHERE
                    name = ?
                LIMIT 1
            ),
            ?,?,?,?,?,?,?,
            (
                SELECT
                    source_id
                FROM
                    source
                WHERE
                    name = ?
                LIMIT 1
            ),
            ?
        )
    };
    
    $sth->execute(
        $vf->variation->dbID(),
        $vf->seq_region_name(),
        $vf->seq_region_start(),
        $vf->seq_region_end(),
        $vf->seq_region_strand()
    );
    my $row = $sth->fetchrow_arrayref;
    my $dbid;
    $dbid = $row->[0] if (defined($row));
    $vf->dbID($dbid);
    
    unless (defined($vf->dbID())) {
        
        $sth = $dbh->prepare($stmt);
        $sth->execute(
            $vf->seq_region_name(),
            $vf->seq_region_start(),
            $vf->seq_region_end(),
            $vf->seq_region_strand(),
            $vf->variation->dbID(),
            $vf->allele_string(),
            $vf->variation->name(),
            $vf->map_weight(),
            $vf->source(),
            (defined($vf->get_all_validation_states()) ? join(",",@{$vf->get_all_validation_states()}) : undef)
        );
        
        $vf->dbID($dbh->last_insert_id(undef,undef,undef,undef));
    }
    
    return $vf->dbID();
}

sub import_alleles {
    my $variation = shift;
    my $dbh = shift;
    
    # For each allele, check whether we already have that allele for the variation, subsnp and sample, in which case we do nothing
    my $check_stmt = qq{
        SELECT
            allele_id
        FROM
            allele
        WHERE
            variation_id = %d AND
            allele = '%s'
# Skip subsnp_id for now            AND subsnp_id %s
        LIMIT 1
    };
    
    my @dbIDs;
    
    my $stmt = qq{
        INSERT IGNORE INTO
            allele (
                variation_id,
                subsnp_id,
                sample_id,
                allele,
                frequency,
                count
            )
        VALUES (
            ?,?,?,?,?,?
        )
    };
    my $sth = $dbh->prepare($stmt);
    
    # Loop over all alleles for this variation object
    foreach my $allele (@{$variation->get_all_Alleles()}) {
        
        # Check if the allele already exists
        $stmt = sprintf($check_stmt,$variation->dbID(),$allele->allele(),(defined($allele->subsnp()) ? "= " . $allele->subsnp() : "IS NULL "));
        my $dbid = $dbh->selectall_arrayref($stmt)->[0][0];
        $allele->dbID($dbid);
        
        # Insert the allele if it wasn't already in the db
        unless (defined($allele->dbID())) {
            $sth->execute(
                $variation->dbID(),
                $allele->subsnp(),
                (defined($allele->population()) ? $allele->population->dbID() : undef),
                $allele->allele(),
                $allele->frequency(),
                $allele->count()
            );
            $allele->dbID($dbh->last_insert_id(undef,undef,undef,undef));
        }
        
        push(@dbIDs,$allele->dbID());
    }
    
    return \@dbIDs;
}

sub import_variation {
    my $variation = shift;
    my $dbh = shift;
    
    # First, check if the variation name already exists. In that case, just set the dbID of our object and return
    my $stmt = qq{
        SELECT
            variation_id
        FROM
            variation
        WHERE
            name = ?
        LIMIT 1
    };
    my $sth = $dbh->prepare($stmt);
    
    $stmt = qq{
        INSERT IGNORE INTO
            variation (
                source_id,
                name,
                validation_status,
                ancestral_allele,
                flipped
            )
        VALUES (
            (
                SELECT
                    source_id
                FROM
                    source
                WHERE
                    name = ?
                LIMIT 1
            ),?,?,?,?
        )
    };
    
    $sth->execute($variation->name());
    my $row = $sth->fetchrow_arrayref;
    my $dbid;
    $dbid = $row->[0] if (defined($row));
    $variation->dbID($dbid);
    
    unless (defined($variation->dbID())) {
    
        $sth = $dbh->prepare($stmt);
        
        $sth->execute(
            $variation->source(),
            $variation->name(),
            (defined($variation->get_all_validation_states()) ? join(",",@{$variation->get_all_validation_states()}) : undef),
            $variation->ancestral_allele(),
            undef
        );
        $variation->dbID($dbh->last_insert_id(undef,undef,undef,undef));
    
    }
    
    return $variation->dbID();
}

sub import_flanking_sequence {
    my $vf = shift;
    my $dbh = shift;
    
    my $stmt = qq{
        INSERT IGNORE INTO
            flanking_sequence (
                variation_id,
                up_seq_region_start,
                up_seq_region_end,
                down_seq_region_start,
                down_seq_region_end,
                seq_region_id,
                seq_region_strand
            )
        VALUES (
            ?,?,?,?,?,
            (
                SELECT
                    seq_region_id
                FROM
                    seq_region
                WHERE
                    name = ?
                LIMIT 1
            ),?
        )
    };
    my $sth = $dbh->prepare($stmt);
    
    $sth->execute(
        $vf->variation->dbID(),
        ($vf->seq_region_start - $DEFAULT_FLANKING_SIZE),
        ($vf->seq_region_start - 1),
        ($vf->seq_region_end + 1),
        ($vf->seq_region_end + $DEFAULT_FLANKING_SIZE),
        $vf->seq_region_name,
        $vf->seq_region_strand
    );
    
    return $vf->variation->dbID();
}

sub import_phenotype {
    my $description = shift;
    my $dbh = shift;
    
    my $stmt = qq{
        SELECT
            phenotype_id
        FROM
            phenotype
        WHERE
            description LIKE ?
        LIMIT 1
    };
    my $check_sth = $dbh->prepare($stmt);
    
    $stmt = qq{
        INSERT IGNORE INTO
            phenotype (
                description
            )
        VALUES (
            ?
        )
    };
    my $sth = $dbh->prepare($stmt);
    
    $check_sth->execute(
        $description
    );
    my $row = $check_sth->fetchrow_arrayref();
    my $dbid;
    $dbid = $row->[0] if (defined($row));
    
    unless (defined($dbid)) {
        $sth->execute(
            $description
        );
        $dbid = $dbh->last_insert_id(undef,undef,undef,undef);
    }
    
    return $dbid;
}

sub check_failed {
    my $vf = shift;
    my $dbh = shift;
    
    my $stmt = qq{
        INSERT IGNORE INTO
            failed_variation (
                variation_id,
                subsnp_id,
                failed_description_id
            )
        VALUES (
            ?,?,?
        )
    };
    my $sth = $dbh->prepare($stmt);
    
    $stmt = qq{
        UPDATE
            variation
        SET
            validation_status = CONCAT(validation_status,',failed');
        WHERE
            variation_id = ?
    };
    my $fail_sth = $dbh->prepare($stmt);
    
    # Get the subsnp id if available
    my $subsnp_id;
    $subsnp_id = $vf->variation->get_all_Alleles->[0]->subsnp() if (defined($vf->variation->get_all_Alleles));
    
    # Check if the reference allele matches the reference sequence
    unless (grep {$_->allele() eq $vf->feature_Slice->seq()} @{$vf->variation->get_all_Alleles()}) {
        warn ("None of the variant alleles for variation_feature_id " . $vf->dbID() . " (" . $vf->variation->name() . ", " . $vf->allele_string . ") match the reference sequence (" . $vf->seq_region_name . ":" . $vf->seq_region_start . "-" . $vf->seq_region_end . ":" . $vf->seq_region_strand . ")");
        $sth->execute($vf->variation->dbID(),$subsnp_id,2);
        # $fail_sth->execute($vf->variation->dbID());
    }
}

sub import_source {
    my $source = shift;
    my $dbh = shift;
    
    my $stmt = qq{
        SELECT 
            source_id
        FROM
            source
        WHERE
            name = '$source'
        LIMIT 1
    };
    my $dbid = $dbh->selectall_arrayref($stmt)->[0][0];
    
    unless (defined($dbid)) {
        $stmt = qq{
            INSERT INTO
                source (
                    name,
                    type
                )
            VALUES (
                '$source',
                'lsdb'
            )
        };
        $dbh->do($stmt);
        $dbid = $dbh->last_insert_id(undef,undef,undef,undef);
    }
    
    return $dbid;
}

sub usage {
    
    print STDOUT qq{
              
    };
    
    exit(0);
}