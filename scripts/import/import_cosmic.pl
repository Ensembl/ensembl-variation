# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

# Import COSMIC data into an Ensembl variation schema database

use strict;
use warnings;

use Getopt::Long;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(SO_variation_class);

my $USE_DB = 1;
my $ADD_SETS = 1;
my $VERBOSE = 1;

my $import_file;
my $registry_file;
my $version;
my $help;

GetOptions(
    "import|i=s"    => \$import_file,
    "registry|r=s"  => \$registry_file,
    "verbose|v"     => \$VERBOSE,
    "version=s"     => \$version,
    "test|t"        => \$USE_DB,
    "help|h"        => \$help,
);

unless (defined($registry_file) && defined($import_file) && defined($version)) {
    print "Must supply an import file, a registry file and a version ...\n" unless $help;
    $help = 1;
}
if ($import_file =~ /\.gz$/) {
  print "Please unzip the input file before running the script!\n";
  exit(0);
}
if ($help) {
    print "Usage: $0 --import <import_file> --registry <reg_file> --version <cosmic_version>\n";
    exit(0);
}

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_all($registry_file);

my $sa = $registry->get_adaptor(
    'human', 'core', 'slice'
);

my $dbh = $registry->get_adaptor(
    'human', 'variation', 'variation'
)->dbc->db_handle;

open my $INPUT, "<$import_file" or die "Can't open '$import_file'";    

my $set_version_sth = $dbh->prepare(qq{UPDATE source SET version = ? WHERE source_id = ?});
    
# check to see if we already have the COSMIC source
my $src_sth = $dbh->prepare(qq{ SELECT source_id FROM source WHERE name LIKE "COSMIC%" });
$src_sth->execute;

my ($source_id) = $src_sth->fetchrow_array;

if ($source_id) {
    print "Found existing source_id: $source_id\n";
    $set_version_sth->execute($version,$source_id);
}
else {
    # if not, add it
    my $sth = $dbh->prepare(qq{
        INSERT INTO source (name, description, url, somatic_status, version) 
        VALUES (
            'COSMIC', 
            'Somatic mutations found in human cancers from the COSMIC project', 
            'http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/',
            'somatic',
            $version
        );
    });

    $sth->execute;

    $source_id = $dbh->last_insert_id(undef, undef, undef, undef);

    print "New source_id: $source_id\n";
}

# set up the various queries and inserts we'll need

my %failed_desc_list = ( 1 => 'Variation has no associated sequence',
                         2 => 'Mapped position is not compatible with reported alleles'
                       );

my $vf_table = 'variation_feature';
my $cosmic_phe_prefix = 'COSMIC:tumour_site:';

# Variations
my $find_existing_var_sth = $dbh->prepare(qq{
    SELECT variation_id FROM variation WHERE name = ?
});

# Variation set COSMIC
my $get_cosmic_set_id_sth = $dbh->prepare(qq{
    SELECT variation_set_id FROM variation_set WHERE name = "COSMIC phenotype variants"
});

$get_cosmic_set_id_sth->execute;

my ($cosmic_set_id) = $get_cosmic_set_id_sth->fetchrow_array;

die "Didn't find COSMIC set id?" unless defined $cosmic_set_id;


# Variation set phenotypes
my $get_phenotype_set_id_sth = $dbh->prepare(qq{
    SELECT variation_set_id FROM variation_set WHERE name like '%phenotype/disease-associated variants%'
});
$get_phenotype_set_id_sth->execute;

my ($phenotype_set_id) = $get_phenotype_set_id_sth->fetchrow_array;

die "Didn't find phenotype/disease set id?" unless defined $phenotype_set_id;


my $add_var_sth = $dbh->prepare(qq{
    INSERT INTO variation (source_id, name, flipped, class_attrib_id, somatic) VALUES (?,?,?,?,1)
});

my $add_var_set_sth = $dbh->prepare(qq{
    INSERT INTO variation_set_variation (variation_id, variation_set_id) VALUES (?,?)
});

my $add_failed_var_sth = $dbh->prepare(qq{
    INSERT INTO failed_variation (variation_id, failed_description_id) VALUES (?,?)
});

my $add_vf_sth = $dbh->prepare(qq{
    INSERT INTO $vf_table (variation_id, seq_region_id, seq_region_start,
        seq_region_end, seq_region_strand, variation_name, allele_string, map_weight, 
        source_id, class_attrib_id, somatic)
    VALUES (?,?,?,?,?,?,?,?,?,?,1)
});

my $select_populations_sth = $dbh->prepare(qq{
    SELECT population_id, name FROM population WHERE name like "COSMIC%"
});

my $add_population_sth = $dbh->prepare(qq{
    INSERT INTO population (name, size, description) 
    VALUES (?,?,?)
});

my $add_allele_with_code_sth = $dbh->prepare(qq{
    INSERT INTO allele (variation_id, population_id, allele_code_id, frequency)
    VALUES (?,?,?,?)
});

my $select_phenotypes_sth = $dbh->prepare(qq{
    SELECT phenotype_id, description FROM phenotype WHERE description LIKE "$cosmic_phe_prefix%"
});

my $add_phenotype_sth = $dbh->prepare(qq{
    INSERT INTO phenotype (description) VALUE (?)
});

my $add_phe_feature_sth = $dbh->prepare(qq{
    INSERT INTO phenotype_feature (phenotype_id, source_id, type, object_id, seq_region_id,
        seq_region_start, seq_region_end, seq_region_strand)
    VALUES (?,?,'Variation',?,?,?,?,?)
});

my $add_phe_feature_attrib_sth = $dbh->prepare(qq{
    INSERT INTO phenotype_feature_attrib (phenotype_feature_id, attrib_type_id, value)
    SELECT ?, attrib_type_id, ? FROM attrib_type WHERE code=?    
});

# allele code id
my $select_allele_codes_sth = $dbh->prepare(qq{
    SELECT allele_code_id, allele FROM allele_code WHERE allele IN ('A','T','G','C','-');
});

my $get_allele_code_sth = $dbh->prepare(qq{
    SELECT allele_code_id
    FROM allele_code
    WHERE allele=?
});
my $count_like_allele_code_sth = $dbh->prepare(qq{
    SELECT count(allele_code_id)
    FROM allele_code
    WHERE allele like ?
});
my $insert_allele_code_sth = $dbh->prepare(qq{
    INSERT INTO allele_code (allele) VALUES (?)
});

# get the failed description ID we will use

my $get_failed_desc_sth = $dbh->prepare(qq{
    SELECT failed_description_id 
    FROM   failed_description
    WHERE  description = ?
});

my $get_class_attrib_ids_sth = $dbh->prepare(qq{
    SELECT a.value, a.attrib_id
    FROM   attrib a, attrib_type at
    WHERE  a.attrib_type_id = at.attrib_type_id
    AND    at.code = 'SO_term'
});

my $patched_version = 0;

my %class_attrib_ids;
my %cosmic_genes_list;
my %cosmic_populations_list;
my %cosmic_phenotypes_list;
my %allele_codes_list;

#### Pre load some of the data ####
# Variant classes list
$get_class_attrib_ids_sth->execute;
while (my ($value, $attrib_id) = $get_class_attrib_ids_sth->fetchrow_array) {
  $class_attrib_ids{$value} = $attrib_id;
}

# Existing COSMIC populations (empty if you used the "remove_cosmic.pl" script before running this script)
$select_populations_sth->execute;
while (my ($p_id, $p_name) = $select_populations_sth->fetchrow_array) {
  $cosmic_populations_list{$p_name} = $p_id;
}

# Existing COSMIC phenotypes
$select_phenotypes_sth->execute;
while (my ($phe_id, $phe_name) = $select_phenotypes_sth->fetchrow_array) {
  $cosmic_phenotypes_list{$phe_name} = $phe_id;
}

# Allele code IDs (A,T,G,C,')
$select_allele_codes_sth->execute;
while (my ($ac_id, $ac_name) = $select_allele_codes_sth->fetchrow_array) {
  $allele_codes_list{$ac_name} = $ac_id;
}


# loop over the input file

MAIN_LOOP : while(<$INPUT>) {

    next unless /COSM/;
  
    my $line = $_;
    
    my (
        $cosmic_version,
        $cosmic_id, 
        $cosmic_gene,
        $accession,
        $class,
        $cds, 
        $aa, 
        $chr, 
        $start, 
        $stop, 
        $mut_nt, 
        $mut_aa,
        $is_snp,
        $tumour_site,
        $mutated_samples,
        $total_samples,
        $freq
    ) = split /,/, $_;


    next unless $cosmic_id =~ /COSM\d+/;
    
    while (!defined($freq)) {
      my $new_line = <$INPUT>;
      chomp $line;
      $line = "$line$new_line";
      
      (
        $cosmic_version,
        $cosmic_id, 
        $cosmic_gene,
        $accession,
        $class,
        $cds, 
        $aa, 
        $chr, 
        $start, 
        $stop, 
        $mut_nt, 
        $mut_aa,
        $is_snp,
        $tumour_site,
        $mutated_samples,
        $total_samples,
        $freq
      ) = split /,/, $line;
    }

    unless ($patched_version) {

        my ($curr_version) = $cosmic_version =~ /COSMIC v(\d+)/;

        die "No version information in file" unless $curr_version;

        $set_version_sth->execute($curr_version, $source_id);

        print "Set COSMIC version to $curr_version\n";
        
        $patched_version = 1;
    }

    my %cosmic_attrib = ('associated_gene' => $cosmic_gene);

    print "$cosmic_id:";
    
    if ($chr == 23) {
        $chr = 'X';
    }
    elsif ($chr == 24) {
        $chr = 'Y';
    }

    # check if we should swap the start & end

    if ($start > $stop + 1) {
        ($start, $stop) = ($stop, $start);
        print "Swapped start & end incorrect coords: start: $start stop: $stop\n";
    }

#    if ($class =~ /insertion/i && $mut_nt =~ /(^-)|^\s*$/ && $cds !~ /dup/i && $cds !~ /del/) {
#        if ($stop == $start+1) {
#            ($start, $stop) = ($stop, $start);
#            print "Swapped start & end for insertion: start: $start stop: $stop\n";
#        }
#        else {
#            print "Wierd start stop for insertion?\n";
#        }
#    }

    my $fail_variation = 0;

    my $dont_rev_comp = 0;

    if ($mut_nt =~ /^\s+$/) {
        #print "$cosmic_id: $cds $class\n";
        $mut_nt = '';
        
        if ($cds =~ /del([A-Z]+)/) {
            $mut_nt = "$1>-";
        }
        elsif ($cds =~ /del([0-9]+)/) {
            $mut_nt = "($1 BP DELETION)>-";
            $dont_rev_comp = 1;
        }
        elsif ($cds =~ /del?/) {
            $mut_nt = "DELETION>-";
            $dont_rev_comp = 1;
        }
        elsif ($cds =~ /ins([A-Z]+)/) {
            $mut_nt = "->$1";
        }
        elsif ($cds =~ /ins([0-9]+)/) {
            $mut_nt = "->($1 BP INSERTION)";
            $dont_rev_comp = 1;
        }
        elsif ($cds =~ /ins?/) {
            $mut_nt = "->INSERTION";
            $dont_rev_comp = 1;
        }
        elsif ($cds =~ /([A-Z]+)>([A-Z]+)/) {
            $mut_nt = "$1>$2";
        }
        elsif ($cds =~ />([A-Z]+)/) {
            $mut_nt = "->$1";
        }
        else {
            $fail_variation = 1;
            print "buggy: failing variation\n";
        }
    }

    $mut_nt = uc($mut_nt);

    # cosmic doesn't prefix insertions with a '-', so add one if necessary

    $mut_nt = '-'.$mut_nt if $mut_nt =~ /^>/;
    
    if ($mut_nt =~ /^-/) {
        if ($stop == $start+1) {
            ($start, $stop) = ($stop, $start);
            print "Swapped start & end for insertion: start: $start stop: $stop\n";
        }
        else {
            if ($start == $stop) {
                $start++;
                print "Insertion start == stop, patched to: start: $start stop: $stop\n";
            }
            else {
                $fail_variation = 2;
                print "Wierd start stop for insertion?\n";
            }
        }
    }

    if ($mut_nt =~ /([A-Z]+)>-$/) {
        my $seq = $1;
        my $seq_len = length($seq);
        if (($seq_len) > 1 && ($stop == $start)) {
            $stop = $start + $seq_len - 1;
            print "Patched stop coordinate for deletion\n";
        }
    }

    my $slice = $sa->fetch_by_region('chromosome', $chr, $start, $stop);

    if ($slice) {

        # try to find a matching ensembl gene

        my $ens_gene;

        my @genes = @{ $slice->get_all_Genes };

        my $check_strand = 0;

        my $strand;

        GENE : for my $gene (@genes) {
        
            # Use existing gene found
            if ($cosmic_genes_list{$cosmic_gene}{$accession}) {
              if ($cosmic_genes_list{$cosmic_gene}{$accession} eq $gene) {
                $ens_gene = $gene;
                last;
              }
            }
        
            if ($gene->external_name eq $cosmic_gene) {
                $ens_gene = $gene;
                last;
            }
            else {
                # try synonyms of this gene
                for my $db_entry (@{ $gene->get_all_DBEntries } ) {
                    for my $syn (@{ $db_entry->get_all_synonyms }) {
                        if ($syn eq $cosmic_gene) {
                            $ens_gene = $gene;
                            last GENE;
                        }
                    }
                }
            }
            
            if ($gene->display_xref && $gene->display_xref->display_id eq $cosmic_gene) {
                $ens_gene = $gene;
                last GENE;
            }

            if ($accession =~ /ENST/) {
                for my $tran (@{ $gene->get_all_Transcripts }) {
                    if ($tran->stable_id eq $accession) {
                        $ens_gene = $gene;
                        last GENE;
                    }
                }
                 
                # Extract the gene name (e.g. CDKN2A_ENST00000361570 => extracts CDKN2A)
                my @g_names = split('_',$cosmic_gene);
                if (scalar @g_names > 1 and $g_names[1] =~ /ENST/) {
                    if ($gene->external_name eq $g_names[0]) {
                        $ens_gene = $gene;
                        last GENE;
                    }
                    else {
                        # try synonyms of this gene
                        for my $db_entry (@{ $gene->get_all_DBEntries } ) {
                            for my $syn (@{ $db_entry->get_all_synonyms }) {
                                if ($syn eq $g_names[0]) {
                                    $ens_gene = $gene;
                                    last GENE;
                                }
                            }
                        }
                    }
                }  
            }            
            
            if ($gene->description && $gene->description =~ /$cosmic_gene/) {
                $ens_gene = $gene;
                last GENE;
            }
            
        }

        if (!$ens_gene && @genes == 1) {
            $ens_gene = $genes[0];
            print "Guess gene is ".$ens_gene->external_name."\n";
            $check_strand = 1;
        }

        if ($ens_gene) {
            $strand = $ens_gene->strand;
            $cosmic_genes_list{$cosmic_gene}{$accession} = $ens_gene if (!$cosmic_genes_list{$cosmic_gene}{$accession});
        }
        else {
            my $num = @genes;

            my $strands_match = 0;

            my $poss_strand;

            if ($num > 0) {
                $strands_match = 1;
                $poss_strand = $genes[0]->strand;
                for my $gene (@genes) {
                    if ($poss_strand != $gene->strand) {
                        $strands_match = 0;
                        $poss_strand = undef;
                        last;
                    }
                }
            }

            $strand = $poss_strand;

            print "Can't find a matching gene for $cosmic_gene or $accession ";
            print "($num genes lie here - strands match: $strands_match)\n";
            
            $check_strand = 1;
        }
        
        my $ens_ref = $slice->seq || '-';

        my $ens_ref_comp = $ens_ref;

        reverse_comp(\$ens_ref_comp);
       
        my ($cs_ref_nt, $cs_mut_nt) = split />/, $mut_nt;
        my ($cs_ref_aa, $cs_mut_aa) = split />/, $mut_aa;

        if ($check_strand && $fail_variation != 1) {
            if ($cs_ref_nt eq $ens_ref) {
                
                if (defined $strand && $strand == 1) {
                    print "Guessed strand 1 looks OK\n";
                }
                else {
                    print "cosmic reference matches ensembl reference, using 1\n";
                    $strand = 1;
                }
            }
            elsif ($cs_ref_nt eq $ens_ref_comp) {
                if (defined $strand && $strand == -1) {
                    print "Guessed strand -1 looks OK\n";
                }
                else {
                    print "cosmic reference matches reverse complemented ensembl reference, using -1\n";
                    $strand = -1;
                }
            }

            unless (defined $strand) {
                print "Reference bases mismatch and no evidence for any strand - leaving\n";
                next;
            }
        }
        
        unless (defined $strand) {
          print "Reference bases mismatch and no evidence for any strand - leaving\n";
          next;
        }

        # these two variables are used when adding this feature to the DB because we 
        # add everything as if it were on the forward strand, regardless of the cosmic 
        # strand

        my $flipped = 0;
        my $var_strand = $strand;
        
        my $forward_strand_ref_allele = $cs_ref_nt;
        my $forward_strand_mut_allele = $cs_mut_nt;

        my $mismatching_reference;

        unless ($fail_variation == 1) {
            if ($cs_ref_nt ne $ens_ref) {
              if ($cs_ref_nt =~ /^[^ATGC]$/ && $class =~ /^Substitution/) {
                if ($cds =~ /c\.\d+(\w)>\w$/) {
                  $cs_ref_nt = uc($1) if ($1 !~ /$cs_ref_nt/i);
                }  
              } elsif ($cs_ref_nt =~ /^\d+$/ && $class =~ /^Deletion/) {
                 $cs_ref_nt = $ens_ref;
              }
            }
        
            my $flag_update_allele = 0;
            
            # Reverse complement if alleles match the reverse strand sequence
            if ($strand == 1) {
              if ($cs_ref_nt ne $ens_ref) {
                unless ($dont_rev_comp) {
                  reverse_comp(\$cs_ref_nt);
                  reverse_comp(\$cs_mut_nt);
                }    
                if ($cs_ref_nt ne $ens_ref) {
                  $mismatching_reference = "cosmic: $cs_ref_nt vs. ensembl: $ens_ref";
                }  
                else {
                  $flipped = 1;
                  $flag_update_allele = 1;
                }
              }
            }
            # Reverse complement and change strand if alleles match the reverse strand sequence
            else {
              if ($cs_ref_nt eq $ens_ref_comp) {
                unless ($dont_rev_comp) {
                  reverse_comp(\$cs_ref_nt);
                  reverse_comp(\$cs_mut_nt);
                  $flipped = 1;
                  $flag_update_allele = 1;
                  $var_strand = 1;
                }    
              } elsif ($cs_ref_nt eq $ens_ref) {
                $var_strand = 1;
              }  else {  
                $mismatching_reference = "cosmic: $cs_ref_nt vs. ensembl: $ens_ref";
              }  
            }
            
            if ($flag_update_allele == 1 && !$mismatching_reference) {
              $forward_strand_ref_allele = $cs_ref_nt;
              $forward_strand_mut_allele = $cs_mut_nt;
            }
        }
       
        if ($mismatching_reference) {
            my $ens_gene_id = ($ens_gene) ? $ens_gene->stable_id : 'Ensembl gene not found';
            print "Reference bases don't match for $cosmic_gene ($ens_gene_id): $mismatching_reference\n";

#            print "Other variations here:\n";
#
#            my $vfs = $slice->get_all_VariationFeatures;
#
#            for my $vf (@$vfs) {
#                print $vf->allele_string, "\n";
#            }
#
#            next;
        }
        
        my $allele_string = ($fail_variation == 1) ? '' : "$forward_strand_ref_allele/$forward_strand_mut_allele";
        
        my $class_id = $class_attrib_ids{
            SO_variation_class($allele_string, defined $mismatching_reference ? 0 : 1)
        };

        die "No class attrib id for allele string: $allele_string" unless defined $class_id;

        # handle large deletions, insertions and indels
        if (defined($forward_strand_ref_allele)) {
          if (length($forward_strand_ref_allele) > 50 || ($stop-$start+1) > 50) {
            $allele_string = ($forward_strand_mut_allele eq '-') ? "LARGE_DELETION" : "LARGE_INDEL";
          }
        }
        if (defined($forward_strand_mut_allele)) {
          if (length($forward_strand_mut_allele) > 50 || ($stop-$start+1) > 50) { 
            $allele_string = ($forward_strand_ref_allele eq '-') ? "LARGE_INSERTION" : "LARGE_INDEL";
          }
        }

        unless ($USE_DB) {
            print " all OK\n";
            next;
        }

        # all OK so far so let's (try to) add stuff to the database

        my $seq_region_id = $sa->get_seq_region_id($slice);

        # first check to see if we already have this variation

        $find_existing_var_sth->execute($cosmic_id);

        my ($variation_id) = $find_existing_var_sth->fetchrow_array;

        if (!$variation_id) {

            # it's new, so add the main variation object,

            $add_var_sth->execute($source_id, $cosmic_id, $flipped, $class_id);

            $variation_id = $dbh->last_insert_id(undef, undef, undef, undef);

            if ($fail_variation) {
                my $failed_description_id = get_failed_description_id($fail_variation);
                $add_failed_var_sth->execute($variation_id, $failed_description_id);
                print "Failed variation\n";
                next MAIN_LOOP if ($fail_variation == 1);
            }
            
            if ($ADD_SETS) {
                $add_var_set_sth->execute($variation_id, $cosmic_set_id);
            }

            # the variation feature object.

            $add_vf_sth->execute(
                $variation_id,
                $seq_region_id,
                $start,
                $stop,
                $var_strand,
                $cosmic_id,
                $allele_string,
                1,
                $source_id,
                $class_id,
            ); 
        }
        
        # add sample and allele information

        # first check if there is an existing population
        
        my $population_id;

        my $sample_name = "COSMIC:gene:$cosmic_gene:tumour_site:$tumour_site";
        
        if ($cosmic_populations_list{$sample_name}) {
            $population_id = $cosmic_populations_list{$sample_name};
        }
        else {
            # there is not, so add it

            $add_population_sth->execute(
                $sample_name, 
                $total_samples,
                "$total_samples $tumour_site tumours examined through the $cosmic_gene gene"
            );

            $population_id = $dbh->last_insert_id(undef, undef, undef, undef);
            $cosmic_populations_list{$sample_name} = $population_id;
        }

        my $mut_freq = $mutated_samples / $total_samples;
          
        # add the mutant allele
        my $mut_allele_code_id = get_allele_code($forward_strand_mut_allele);
        $add_allele_with_code_sth->execute(
            $variation_id,
            $population_id,
            $mut_allele_code_id,
            $mut_freq,
        );
        
        # and the reference allele
        my $ref_allele_code_id = get_allele_code($forward_strand_ref_allele);
        $add_allele_with_code_sth->execute(
            $variation_id,
            $population_id,
            $ref_allele_code_id,
            (1 - $mut_freq)
        );

        # try to find an existing phenotype, or create a new phenotype entry

        my $phenotype_id;
        
        my $phenotype_name = "$cosmic_phe_prefix$tumour_site";

        if ($cosmic_phenotypes_list{$phenotype_name}) {
            $phenotype_id = $cosmic_phenotypes_list{$phenotype_name};
        }
        else {
            $add_phenotype_sth->execute($phenotype_name);
            $phenotype_id = $dbh->last_insert_id(undef, undef, undef, undef);
            $cosmic_phenotypes_list{$phenotype_name} = $phenotype_id;
        }

        # add variation annotation (phenotype_feature)

        $add_phe_feature_sth->execute(
          $phenotype_id,
          $source_id,
          $cosmic_id,
          $seq_region_id,
          $start,
          $stop,
          1
        );
        my $pf_id = $dbh->last_insert_id(undef, undef, undef, undef);
        
        # add phenotype feature attributes
        foreach my $attrib (keys(%cosmic_attrib)) {
          $add_phe_feature_attrib_sth->execute($pf_id,$cosmic_attrib{$attrib},$attrib);
        }
        print " added to DB\n";
    }
    else {
        print "Failed to get slice for chr $chr $start - $stop\n";
        next;
    }
}


sub get_failed_description_id {
  my $failed_id = shift;
  my $description = $failed_desc_list{$failed_id};
  
  $get_failed_desc_sth->execute($description);
  my ($failed_description_id) = $get_failed_desc_sth->fetchrow_array;
  
  return $failed_description_id;
}


sub get_allele_code {
    my $allele = shift;
    
    my $allele_code_id;

    # Check if an entry exists in the allele_code hash
    if ($allele_codes_list{$allele}) { # A,T,G,C,-
      $allele_code_id = $allele_codes_list{$allele};
    }
    # Check if an entry exists in the allele_code table
    else {
      $get_allele_code_sth->execute($allele);
      ($allele_code_id) = $get_allele_code_sth->fetchrow_array;
    }
    
    # Check with large allele, because of issue with the allele index 
    # (limited to unique 1000 first bases)
    if (!defined($allele_code_id) && length($allele) >= 1000) {
      my $indexed_string = substr($allele,0,1000).'%';
      $count_like_allele_code_sth->execute($indexed_string);
      my ($count_allele) = $count_like_allele_code_sth->fetchrow_array;
      if ($count_allele > 0) {
        $allele = '(LARGE ALLELE)';
        $get_allele_code_sth->execute($allele);
        ($allele_code_id) = $get_allele_code_sth->fetchrow_array;
      }
    }
    
    # Insert a new allele in the allele_code table
    if (!defined($allele_code_id)) {
      $insert_allele_code_sth->execute($allele);
      $allele_code_id = $dbh->last_insert_id(undef, undef, undef, undef);
    }
    
    return $allele_code_id;
}
