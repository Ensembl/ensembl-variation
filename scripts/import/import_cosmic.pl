use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

my $USE_DB = 1;
my $CHECK_CONS = 0;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_all;

my $vfa = $registry->get_adaptor(
    'human',
    'variation',
    'variationfeature'
);

my $sa = $registry->get_adaptor(
    'human', 'core', 'slice'
);

my $dbh = DBI->connect(
    'DBI:mysql:database=grsr_cosmic_import;host=ens-variation;port=3306',
    #'DBI:mysql:database=grsr_cosmic_test;host=127.0.0.1;port=13306',
    'ensadmin',
    'ensembl',
);

# check to see if we already have the COSMIC source

my $src_sth = $dbh->prepare(qq{
    SELECT source_id
    FROM source
    WHERE name LIKE "COSMIC%"
});

$src_sth->execute;

my $existing_src = $src_sth->fetchrow_arrayref;

my $source_id;

if ($existing_src) {
    $source_id = $existing_src->[0];
    
    print "Found existing source_id: $source_id\n";
}
else {

    # if not, add it

    my $sth = $dbh->prepare(qq{
        INSERT INTO source (name, description, somatic) 
        VALUES ('COSMIC', 'Somatic mutations found in human cancers from the COSMIC project', 1);
    });

    $sth->execute;

    $source_id = $dbh->last_insert_id(undef, undef, undef, undef);

    print "New source_id: $source_id\n";
}

# set up the various queries and inserts we'll need

my $find_existing_var_sth = $dbh->prepare(qq{
    SELECT variation_id FROM variation WHERE name = ?
});

my $add_var_sth = $dbh->prepare(qq{
    INSERT INTO variation (source_id, name) VALUES (?,?)
});

my $add_flank_sth = $dbh->prepare(qq{
    INSERT INTO flanking_sequence (variation_id, seq_region_id, seq_region_strand, 
        up_seq_region_start, up_seq_region_end, down_seq_region_start, down_seq_region_end)
    VALUES (?,?,?,?,?,?,?)
});

my $add_vf_sth = $dbh->prepare(qq{
    INSERT INTO variation_feature (variation_id, seq_region_id, seq_region_start,
        seq_region_end, seq_region_strand, variation_name, allele_string, map_weight, 
        source_id)
    VALUES (?,?,?,?,?,?,?,?,?)
});

my $find_existing_sample_sth = $dbh->prepare(qq{
    SELECT sample_id FROM sample WHERE name = ? AND size = ? 
});

my $add_sample_sth = $dbh->prepare(qq{
    INSERT INTO sample (name, size, description) 
    VALUES (?,?,?)
});

my $add_allele_sth = $dbh->prepare(qq{
    INSERT INTO allele (variation_id, sample_id, allele, frequency)
    VALUES (?,?,?,?)
});

my $find_matching_phenotype_sth = $dbh->prepare(qq{
    SELECT * FROM phenotype WHERE description LIKE ? ORDER BY phenotype_id LIMIT 1
});

my $add_phenotype_sth = $dbh->prepare(qq{
    INSERT INTO phenotype (description) VALUE (?)
});

my $add_annotation_sth = $dbh->prepare(qq{
    INSERT INTO variation_annotation (variation_id, phenotype_id, source_id, 
        associated_gene, variation_names)
    VALUES (?,?,?,?,?)
});

# loop over the input file

while(<>) {
    
    my (
        $cosmic_id, 
        $cosmic_gene, 
        $_cds, 
        $_aa, 
        $chr, 
        $start, 
        $stop, 
        $mut_nt, 
        $mut_aa,
        $tumour_site,
        $mutated_samples,
        $total_samples
    ) = split /,/;

    next unless $cosmic_id =~ /COSM\d+/;

    print "$cosmic_id: ";
    
    if ($stop != $start) {
        print "Mutation covers more than 1 bp\n";
        next;
    }

    $chr = 'X' if $chr == 23;

    my $slice = $sa->fetch_by_region('chromosome', $chr, $start, $stop);

    if ($slice) {

        # try to find a matching ensembl gene

        my $ens_gene;

        GENE : for my $gene (@{ $slice->get_all_Genes } ) {
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
        }
        
        unless ($ens_gene) {
            print "Can't find a matching gene for $cosmic_gene\n";
            next;
        }
        
        my $ens_ref = $slice->seq;

        my $ens_ref_comp = $ens_ref;

        reverse_comp(\$ens_ref_comp);
       
        my ($cs_ref_nt, $cs_mut_nt) = split />/, $mut_nt;
        my ($cs_ref_aa, $cs_mut_aa) = split />/, $mut_aa;

        my $strand = $ens_gene->strand;

        # these two variables are used when adding this feature to the DB because we 
        # add everything as if it were on the forward strand, regardless of the cosmic 
        # strand

        my $forward_strand_ref_allele = $cs_ref_nt;
        my $forward_strand_mut_allele = $cs_mut_nt;

        my $mismatching_reference;

        if ($strand == 1) {
            if ($cs_ref_nt ne $ens_ref) {
                $mismatching_reference = "cosmic: $cs_ref_nt vs. ensembl: $ens_ref";
            }
        }
        else {
            reverse_comp(\$forward_strand_ref_allele);
            reverse_comp(\$forward_strand_mut_allele);
            
            if ($cs_ref_nt ne $ens_ref_comp) {
                $mismatching_reference = "cosmic: $cs_ref_nt vs. ensembl: $ens_ref_comp";
            }
        }
       
        if ($mismatching_reference) {
            print "Reference bases don't match for $cosmic_gene (".
                $ens_gene->stable_id."): $mismatching_reference";

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

        my %trans_gene = map { $_->dbID => 1 } @{ $ens_gene->get_all_Transcripts };

        if ($CHECK_CONS) {

            my $vf = Bio::EnsEMBL::Variation::VariationFeature->new(
                -start          => 1,
                -end            => 1,
                -slice          => $slice,        
                -allele_string  => "$cs_ref_nt/$cs_mut_nt",
                -strand         => $strand,
                -map_weight     => 1,
                -adaptor        => $vfa,       
                -variation_name => $cosmic_id,
            );

            my $found_equiv = 0;

            my $cs_pep = "$cs_ref_aa/$cs_mut_aa";

            for my $tv (@{ $vf->get_all_TranscriptVariations }) {
                if ($trans_gene{$tv->transcript->dbID}) {
                    my $ens_pep = $tv->pep_allele_string;
                    if ($ens_pep =~ /^([A-Z])$/) {
                        $ens_pep = "$1/$1"
                    }
                    #print "$ens_pep vs. $cs_pep\n";
                    if ($ens_pep eq $cs_pep) {
                        $found_equiv++;
                        last;
                    }
                }
            }

            unless ($found_equiv) {
                print "No consequences match\n";
                next;
            }

            print "consequence looks good - ";
        }

        unless ($USE_DB) {
            print " all OK\n";
            next;
        }

        # all OK so far so let's (try to) add stuff to the database

        eval {

            # first check to see if we already have this variation
            
            $find_existing_var_sth->execute($cosmic_id);

            my $existing_var = $find_existing_var_sth->fetchrow_arrayref;

            my $variation_id;

            if ($existing_var) {
                $variation_id = $existing_var->[0];
            }
            else {

                # it's new, so add the main variation object,

                $add_var_sth->execute($source_id, $cosmic_id);

                $variation_id = $dbh->last_insert_id(undef, undef, undef, undef);

                # the flanking region (just using the reference sequence),

                my $seq_region_id = $sa->get_seq_region_id($slice);

                $add_flank_sth->execute(
                    $variation_id, 
                    $seq_region_id,
                    1,
                    $start - 100,
                    $start - 1,
                    $stop + 1,
                    $stop + 100
                );

                # and the variation feature object.

                $add_vf_sth->execute(
                    $variation_id,
                    $seq_region_id,
                    $start,
                    $stop,
                    1,
                    $cosmic_id,
                    "$forward_strand_ref_allele/$forward_strand_mut_allele",
                    1,
                    $source_id
                ); 
            }

            # add sample and allele information

            # first check if there is an existing sample
            
            my $sample_id;

            my $sample_name = "COSMIC:gene:$cosmic_gene:tumour_site:$tumour_site";
            
            $find_existing_sample_sth->execute($sample_name, $total_samples);

            my $existing_sample = $find_existing_sample_sth->fetchrow_arrayref;
            
            if ($existing_sample) {
                $sample_id = $existing_sample->[0];
            }
            else {
                # there is not, so add it

                $add_sample_sth->execute(
                    $sample_name, 
                    $total_samples,
                    "$total_samples $tumour_site tumours examined through the $cosmic_gene gene"
                );

                $sample_id = $dbh->last_insert_id(undef, undef, undef, undef);
            }

            my $mut_freq = $mutated_samples / $total_samples;
            
            # add the mutant allele

            $add_allele_sth->execute(
                $variation_id,
                $sample_id,
                $forward_strand_mut_allele,
                $mut_freq,
            );

            # and the reference allele

            $add_allele_sth->execute(
                $variation_id,
                $sample_id,
                $forward_strand_ref_allele,
                (1 - $mut_freq),
            );

            # try to find an existing phenotype, or create a new phenotype entry

            my $phenotype_id;

            #$find_matching_phenotype_sth->execute('%'.$tumour_site.'%cancer%');

            #my $matching_phen = $find_matching_phenotype_sth->fetchrow_arrayref;

            #if ($matching_phen) {
            #    $phenotype_id = $matching_phen->[0];
            #}
            #else {
            
            my $phenotype_name = "COSMIC:tumour_site:$tumour_site";

            $find_matching_phenotype_sth->execute($phenotype_name);
            
            my $existing_cosmic_phen = $find_matching_phenotype_sth->fetchrow_arrayref;
            
            if ($existing_cosmic_phen) {
                $phenotype_id = $existing_cosmic_phen->[0];
            }
            else {
                $add_phenotype_sth->execute($phenotype_name);
                $phenotype_id = $dbh->last_insert_id(undef, undef, undef, undef);
            }
            
            #}

            # add variation annotation

            $add_annotation_sth->execute(
                $variation_id,
                $phenotype_id,
                $source_id,
                $cosmic_gene,
                $cosmic_id
            );

        }; # end eval around all DB operations
        
        if ($@) {
            print "failed to add to DB: $@\n";
            next;
        }

        print " added to DB\n";
    }
    else {
        print "Failed to get slice for chr $chr $start - $stop\n";
        next;
    }
}


