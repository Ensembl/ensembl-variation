#!/usr/bin/env perl

=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/legal/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut


=head

    This script will post-process the variation_feature table to update the variation_set_id column to contain the
    primary keys of the variation_sets which the corresponding variation is part of.
    
=cut


use strict;
use warnings;
use Getopt::Long;
use DBI;
use ImportUtils qw(update_table);

use Bio::EnsEMBL::Registry;

our $MAX_VARIATION_SETS = 64;

my @option_defs = (
  'species=s',
  'group=s',
  'registry_file=s',
	'sv!',
  'clean!',
  'quiet!',
  'help!'
);

my %options;
GetOptions(\%options,@option_defs);

my $species = $options{'species'};
my $group = $options{'group'} || q{variation};
my $registry_file = $options{'registry_file'};
my $clean = $options{'clean'};
my $quiet = $options{'quiet'};
my $help = $options{'help'};
my $sv_prefix = $options{'sv'} ? 'structural_' : '';

usage() if ($help);

# This is just for convenience in case we have a custom name for the variation_feature table
our $VARIATION_FEATURE_TABLE = $sv_prefix.'variation_feature';
our $VAR_COL       = $sv_prefix.'variation_id';
our $VAR_SET_TABLE = 'variation_set_'.$sv_prefix.'variation';

#ÊCheck that required parameters were passed
die ("Required argument '-species' was not specified") unless (defined($species));
die ("Required argument '-registry_file' was not specified") unless (defined($registry_file));

# Load the registry from the supplied file
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all($registry_file);

#ÊGet a dbadaptor to the variation database
my $dbVar = $registry->get_DBAdaptor($species,$group) or die ("Could not get variation DBAdaptor for $species and $group");

#ÊCall the post-processing subroutine
post_process($dbVar,$clean,$quiet);

sub post_process {
    my $dbVar = shift;
    my $clean = shift;
    my $quiet = shift;
    
    my $stmt; 
    
    #ÊFirst of all, make sure that no primary key for variation_sets is too large to fit inside the set
    $stmt = qq{
        SELECT
            MAX(variation_set_id)
        FROM
            variation_set
    };
    my $max_id = $dbVar->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
    die ("There are variation_sets with primary keys greater than $MAX_VARIATION_SETS. This will not be possible to represent with the MySQL SET data type in the $VARIATION_FEATURE_TABLE!") if ($max_id > $MAX_VARIATION_SETS);

    ###Ê
    ## Define the statements we will be using
    my $tmp_table = ($sv_prefix ne '') ? q{tmp_vs_svf_upd} : q{tmp_vs_vf_upd};    
    
    # Create a temporary table to hold the variation_set_id string until we're ready to update variation_feature
    $stmt = qq{
        CREATE TABLE
            $tmp_table (
                $VAR_COL INT NOT NULL, 
                variation_set_id SET(
                    '1','2','3','4','5','6','7','8',
                    '9','10','11','12','13','14','15','16',
                    '17','18','19','20','21','22','23','24',
                    '25','26','27','28','29','30','31','32',
                    '33','34','35','36','37','38','39','40',
                    '41','42','43','44','45','46','47','48',
                    '49','50','51','52','53','54','55','56',
                    '57','58','59','60','61','62','63','64'
                ) NOT NULL DEFAULT '', 
                PRIMARY KEY ($VAR_COL) 
            );
    };
    my $tmp_tbl_sth = $dbVar->dbc->prepare($stmt);
    
    # Insert the variation_id and a comma-separated list of sets it explicitly belongs to into the temp_table
    $stmt = qq{
        INSERT INTO
            $tmp_table (
                $VAR_COL,
                variation_set_id
            )
        SELECT
            vsv.$VAR_COL,
            GROUP_CONCAT(vsv.variation_set_id)
        FROM
            $VAR_SET_TABLE vsv
        GROUP BY
            vsv.$VAR_COL
    };
    my $ins_expl_sth = $dbVar->dbc->prepare($stmt);
    
    #ÊAdd the implicit parent sets to the list of variation_sets
    $stmt = qq{
        UPDATE
            $tmp_table t,
            variation_set_structure vss
        SET
            t.variation_set_id = CONCAT(
                t.variation_set_id,
                ',',
                vss.variation_set_super
            )
        WHERE
            FIND_IN_SET(
                vss.variation_set_sub,
                t.variation_set_id
            ) AND 
            NOT FIND_IN_SET(
                vss.variation_set_super,
                t.variation_set_id
            )
    };
    my $upd_impl_sth = $dbVar->dbc->prepare($stmt);
    
    ####
    ## Post-process
    
    #ÊFirst, create the temporary table (drop it if it already exists)
    $stmt = qq{
        DROP TABLE IF EXISTS
            $tmp_table
    };
    print STDOUT localtime() . "\tDropping temp table $tmp_table if it exists..." unless ($quiet);
    $dbVar->dbc->do($stmt);
    print STDOUT "done!\n" unless ($quiet);
    
    print STDOUT localtime() . "\tCreating temp table $tmp_table..." unless ($quiet);
    $tmp_tbl_sth->execute();
    print STDOUT "done!\n" unless ($quiet);
    
    # Next, insert the explicitly assigned variation_set_ids
    print STDOUT localtime() . "\tAdding explicitly assigned variation_set primary keys to $tmp_table..." unless ($quiet);
    $ins_expl_sth->execute();
    print STDOUT "done!\n" unless ($quiet);
    
    # Now, insert the implicit variation_set relationships. Because the sets can be nested to any level, we need to iterate until no more rows are updated
    my $update_count = 1;
    while ($update_count) {
        
        print STDOUT localtime() . "\tAdding implicitly assigned variation_set primary keys to $tmp_table..." unless ($quiet);
        $upd_impl_sth->execute();
        print STDOUT "done!\n" unless ($quiet);
        $update_count = $upd_impl_sth->rows();
        print STDOUT "\t$update_count rows were updated, will " . ($update_count ? "" : "NOT ") . "repeat the statement\n" unless ($quiet);
        
    }
    
    # Finally, update the variation_set_id column in variation_feature
    print STDOUT localtime() . "\tUpdating $VARIATION_FEATURE_TABLE with the variation_set_id column from $tmp_table..." unless ($quiet);
    
		update_table($dbVar->dbc, $tmp_table, ($VARIATION_FEATURE_TABLE, $VAR_COL, $VAR_COL, 'variation_set_id', 'variation_set_id'), $clean);
    
    print STDOUT "done!\n" unless ($quiet);
    
    # ...and lastly, drop the temporary table
    print STDOUT localtime() . "\tDropping the temp table $tmp_table..." unless ($quiet);
    $stmt = qq{
        DROP TABLE
            $tmp_table
    };
    $dbVar->dbc->do($stmt);
    print STDOUT "done!\n" unless ($quiet);
    
    # And that's it
    print STDOUT localtime() . "\tPost-processing complete!\n" unless ($quiet);
}

sub usage {
    
    print STDOUT qq{
Usage:

  $0 
    -registry_file f 
    -species s 
    -group g 
    -quiet 
    
Description:

  Post-process the variation_feature table to update the variation_set_id column to be a
  list of the primary keys of all variation_sets which the corresponding variation belongs to.

Command line switches:

  -registry_file f  (Required)
                    An Ensembl registry configuration file with connection details to the
                    relevant database.

  -species s        (Required)
                    The species of the database for which to do the post-processing.

  -group g          (Optional)
                    The group identifier for the database to post-process as it is specified
                    in the registry configuration file. Defaults to 'variation'.

  -sv               (Optional)
                    If specified, will run the script for the structural_variation_feature table
	                
  -clean            (Optional)
                    If specified, all pre-existing variation_set_id entries in variation_feature are truncated. 
                    
  -quiet            (Optional)
                    Suppresses output.
                    
  -help             Displays this message.
        
};
    
    exit(0);
}
