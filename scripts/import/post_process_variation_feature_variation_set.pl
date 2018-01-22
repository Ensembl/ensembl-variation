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


=head

    This script will post-process the variation_feature table to update the variation_set_id column to contain the
    primary keys of the variation_sets which the corresponding variation is part of.
    
=cut


use strict;
use warnings;
use Getopt::Long;
use DBI;
use ImportUtils qw(update_table load);

use Bio::EnsEMBL::Registry;

our $MAX_VARIATION_SETS = 64;

my @option_defs = (
  'species=s',
  'group=s',
  'registry_file=s',
	'sv!',
  'clean!',
  'quiet!',
  'help!',
	'tmpdir=s',
	'tmpfile=s',
  'no_mtmp',
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
my $no_mtmp = $options{'no_mtmp'};

my $TMP_DIR = $options{tmpdir} || '/tmp/';
my $TMP_FILE = $options{tmpfile} || $$.'_vfvs.txt';

$ImportUtils::TMP_DIR = $TMP_DIR;
$ImportUtils::TMP_FILE = $TMP_FILE;

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

    unless( defined $no_mtmp){
	
	# create MTMP table
	my $mtmp_table_name = 'MTMP_variation_set_'.$sv_prefix.'variation';
	
	$stmt = qq{
		CREATE TABLE IF NOT EXISTS $mtmp_table_name (
			$VAR_COL int(11) unsigned NOT NULL,
			variation_set_id int(11) unsigned NOT NULL,
			KEY $VAR_COL ($VAR_COL),
			KEY variation_set_id (variation_set_id)
		) engine=MyISAM
	};
    $dbVar->dbc->do($stmt);
	
	# truncate it
    $dbVar->dbc->do(qq{TRUNCATE $mtmp_table_name});
	
	# now populate it
	$stmt = qq{
		SELECT $VAR_COL, variation_set_id
		FROM $tmp_table
	};

	print STDOUT localtime() . "\tDumping data for $mtmp_table_name\n";
	
	my $retrieve_sth = $dbVar->dbc->prepare($stmt, {mysql_use_result => 1});
	$retrieve_sth->execute;
	
	my ($var_id, $set_ids);
	$retrieve_sth->bind_columns(\$var_id, \$set_ids);
	
	open TMP, ">$TMP_DIR/$TMP_FILE" or die "Could not write to tmp file $TMP_DIR/$TMP_FILE";
	
	while($retrieve_sth->fetch) {
		print TMP "$var_id\t$_\n" for split(',', $set_ids);
	}
	
	close TMP;
	
	print STDOUT localtime() . "\tLoading table $mtmp_table_name from dumped data\n";
	
	# Some other post-processing "cleaning" queries (variation and variation_feature)
	if (!$options{'sv'}) {
	  foreach my $col ('evidence_attribs', 'clinical_significance') {
            $dbVar->dbc->do(qq[update variation set $col = NULL where $col = '';]);
            $dbVar->dbc->do(qq[update $VARIATION_FEATURE_TABLE set $col = NULL where $col = '';]);
            print STDOUT localtime() . "\tCleaning column $col from the variation and $VARIATION_FEATURE_TABLE tables ..." unless ($quiet);
          }
        }

	load($dbVar->dbc, $mtmp_table_name);
	
    }
    
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

  -no_mtmp          (Optional)
                    Don't create the MTMP_variation_set_variation table
        
};
    
    exit(0);
}
