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

# Build a test database based on a variation database

use Bio::EnsEMBL::Registry;

# The default number of variations and structural variations to build the test database from
my $VARIATION_SIZE  = 300000;
my $STRUCT_VAR_SIZE = 100000;

# Get the registry configuration file from the command line
my $registry_file = shift;
die ("You must specify a registry configuration file") unless (defined($registry_file));

# Get the species to use
my $species = shift;
die ("You must specify the species") unless (defined($species));

# Get the variation database foreign key relationships file
my $foreign_key_file = shift;
die ("You must supply the variation database foreign key file") unless (defined($foreign_key_file));

# If a file of variation ids is specified, use that
my $variation_id_file = shift;
warn ("Will use variation ids from $variation_id_file") if (defined($variation_id_file));

# Parse the foreign keys into a hash. Keep the tables that point to variation_id in variation and 
# structural_variation_id in structural_variation
my %foreign_keys;
my @variation_tables;
my @struct_var_tables;
open (KEY,"<",$foreign_key_file) or die ("Could not open foreign key file $foreign_key_file for parsing");
while (<KEY>) {
    chomp;
    
    my ($table,$col,$ref_table,$ref_col) = $_ =~ m/^ALTER\s+TABLE\s+(\S+)\s+ADD\s+FOREIGN\s+KEY\s+\((\S+)\)\s+REFERENCES\s+([^\s\(]+)\s*\(([^\s\)]+)\)/i;
    next unless (defined($table) && defined($col) && defined($ref_table) && defined($ref_col));
 
    $foreign_keys{$table} = [] unless (exists($foreign_keys{$table}));
    push(@{$foreign_keys{$table}},[$col,$ref_table,$ref_col]);
    
    push(@variation_tables, [$table, $col, $ref_col]) if ($ref_table eq 'variation' && $ref_col eq 'variation_id');
    push(@variation_tables, [$table, $col, $ref_col]) if ($table eq 'phenotype_feature' && $ref_table eq 'variation' );

    push(@struct_var_tables,[$table, $col, $ref_col]) if ($ref_table eq 'structural_variation' && $ref_col eq 'structural_variation_id');
    push(@struct_var_tables,[$table, $col, $ref_col]) if ($table eq 'phenotype_feature' && $ref_table eq 'structural_variation' );

}
close (KEY);

# Load the registry
my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_all($registry_file);

# Get db adaptors to the source and destination databases
my $source = $reg->get_DBAdaptor($species,'source_db') or die ("Could not get a DBAdaptor for the source db");
my $dest   = $reg->get_DBAdaptor($species,'dest_db') or die ("Could not get a DBAdaptor for the destination db");


#### Variation ####

# If a variation id file was specified, parse the variation_ids to use from there.
my @variation_ids;
if (defined($variation_id_file)) {
    
    open(IDS,"<",$variation_id_file) or die ("Could not open $variation_id_file for reading");
    while (<IDS>) {
        chomp;
        my ($id) = $_ =~ m/^(\d+)/; 
        push(@variation_ids,$id);
    }
    close(IDS);
} 
# Otherwise, get some random variation ids from the source database
else {
        
        # Get random variation ids              
        @variation_ids = get_random_ids('variation', $VARIATION_SIZE , 'variation');
}

# Populate the destination database
populate_data_table('variation', \@variation_tables, \@variation_ids);

# Create and fill the subsnp_map table
$dest->dbc->do(qq{ CREATE TABLE IF NOT EXISTS subsnp_map (
                   variation_id int(11) unsigned NOT NULL,
                   subsnp_id int(11) unsigned NOT NULL,
                   PRIMARY KEY (variation_id,subsnp_id),
                   KEY variation_idx (variation_id)
                 )});

my $stmt = qq{
     INSERT IGNORE INTO subsnp_map ( variation_id, subsnp_id )
     SELECT DISTINCT variation_id, subsnp_id FROM allele
   };

$dest->dbc->do($stmt);



#### Structural variation ####

# Get random structural variation ids
my @struct_var_ids = get_random_ids('structural_variation', $STRUCT_VAR_SIZE, 'structural variation');

# Populate the destination database
populate_data_table('structural_variation', \@struct_var_tables, \@struct_var_ids);

# Fill the study table (for structural variation data)
# Get the columns of the foreign table
my $cols = get_table_columns($dest,'study');
my $source_db = $source->dbc->dbname();    
# Create the insert statement
my $ins_col_str = join(',',@{$cols});
my $sel_col_str = 'src.' . join(', src.',@{$cols});
my $ins_stmt = qq{
        INSERT IGNORE INTO
        study ($ins_col_str)
    SELECT
        $sel_col_str
        FROM
        structural_variation dst,
                        $source_db\.study src
                WHERE
                        src.study_id = dst.study_id
        };
my $ins_sth = $dest->dbc->prepare($ins_stmt);
$ins_sth->execute();
$ins_sth->finish();


# Lastly, process all the foreign key relationships without recursion
my $srcdb = $source->dbc->dbname();
foreach my $table (keys(%foreign_keys)) {
  foreach my $row (@{$foreign_keys{$table}}) {
    add_foreign_data($dest,$table,$row->[0],$row->[1],$row->[2],$srcdb,{}); 
  }     
}


# Fill the phenotype_ontology_accession table
# Get the columns of the foreign table
my $poa_table = 'phenotype_ontology_accession';
my $poa_cols = get_table_columns($dest,$poa_table);
my $source_db = $source->dbc->dbname();    
# Create the insert statement
my $poa_ins_col_str = join(',',@{$poa_cols});
my $poa_sel_col_str = 'src.' . join(', src.',@{$poa_cols});
my $poa_ins_stmt = qq{
  INSERT IGNORE INTO
    $poa_table ($poa_ins_col_str)
    SELECT
      $poa_sel_col_str
    FROM
      phenotype dst,
      $source_db\.$poa_table src
    WHERE
      src.phenotype_id = dst.phenotype_id
};
my $poa_ins_sth = $dest->dbc->prepare($poa_ins_stmt);
$poa_ins_sth->execute();
$poa_ins_sth->finish();


###########
# Methods #
###########

sub get_random_ids {
        my $table = shift;
        my $size  = shift;
        my $type  = shift;
        
        my $column = "$table\_id";
        my @ids_list;
        
        warn (localtime() . "\tWill extract $size random $type"."s from source database");
    
        # Get the min and max variation_ids
        my $stmt = qq{
                        SELECT
                                        MIN($column),
                                        MAX($column)
                        FROM
                                        $table
        };
        my ($min,$max) = @{$source->dbc->db_handle->selectall_arrayref($stmt)->[0]};
    
        # Randomly pull ids from the database until we've reached the desired number of rows
        $stmt = qq{
                        SELECT $column
                        FROM   $table
                        WHERE  $column = ?
                        LIMIT  1
        };
        my $sth = $source->dbc->prepare($stmt);
        my %ids;
        while (scalar(keys(%ids)) < $size) {
                my $id = $min + int(rand($max-$min+1));
                $sth->execute($id);
                $sth->bind_columns(\$id);
                $sth->fetch();
                next unless defined($id);
                $ids{$id}++;
        }
        return keys(%ids);
}


sub populate_data_table {
        my $table       = shift;
        my $asso_tables = shift;
        my $table_ids   = shift;
        my $ref_column  = "$table\_id";
        
        my $srcdb = $source->dbc->dbname();
        my $cols = get_table_columns($source,$table);
        my $sel_col_str = 'tbl.' . join(', tbl.',@{$cols});
        my $ins_col_str = join(',',@{$cols});

        my $ins_stmt = qq{
    INSERT INTO
        $table ($ins_col_str)
    SELECT
        $sel_col_str
    FROM
        $srcdb\.$table tbl
    WHERE
        tbl.$ref_column = ? 
        };
        print "$ins_stmt\n";
        my $ins_sth = $dest->dbc->prepare($ins_stmt);

        map {$ins_sth->execute($_)} @$table_ids;
    
        # First we'll import the "core" data, i.e. data from the tables having $ref_column as a foreign key
        foreach my $row (@$asso_tables) {
    
    my $a_table = $row->[0];
    my $column = $row->[1];
    my $ref_column2 = $row->[2];
    # Recursively add foreign relationships to this table
    add_foreign_data($dest,$table,$ref_column2,$a_table,$column,$srcdb,\%foreign_keys);
    
        }
}


sub add_foreign_data {
    my $dba = shift;
    my $table = shift;
    my $column = shift;
    my $foreign_table = shift;
    my $foreign_column = shift;
    my $source_db = shift;
    my $foreign_keys = shift;
    
    warn (localtime() . "\tProcessing foreign relationship between $table ($column) and $foreign_table ($foreign_column)");
    
    # Get the columns of the foreign table
    my $cols = get_table_columns($dba,$foreign_table);
    
    # Get the join condition
    my $constraint = qq{src.$foreign_column = dst.$column};
    
    # Create the insert statement
    my $ins_col_str = join(',',@{$cols});
    my $sel_col_str = 'src.' . join(', src.',@{$cols});
    my $ins_stmt = qq{
        INSERT IGNORE INTO
            $foreign_table ($ins_col_str)
        SELECT
            $sel_col_str
        FROM
            $table dst,
            $source_db\.$foreign_table src
        WHERE
            $constraint
    };
    my $ins_sth = $dba->dbc->prepare($ins_stmt);
    my $rc = $ins_sth->execute();
    
    # Recursively add foreign data but if no new rows were added, we'll skip the recursion
    if ($rc > 0 && exists($foreign_keys->{$foreign_table})) {
        foreach my $row (@{$foreign_keys->{$foreign_table}}) {
        
            add_foreign_data($dba,$foreign_table,$row->[0],$row->[1],$row->[2],$source_db,$foreign_keys);
        
        }
    }
}


sub get_table_columns {
    my $dba = shift;
    my $table = shift;
    
    my $stmt = qq{
        SHOW COLUMNS FROM 
            $table
    };
    my $cols = $dba->dbc->db_handle->selectcol_arrayref($stmt);
    return $cols;
}

