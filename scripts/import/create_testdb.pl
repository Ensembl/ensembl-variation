#! perl

use strict;
use warnings;

#ÊBuild a test database based on a variation database

use Bio::EnsEMBL::Registry;

# The default number of variations to build the test database from
my $VARIATION_SIZE = 300000;

# Get the registry configuration file from the command line
my $registry_file = shift;
die ("You must specify a registry configuration file") unless (defined($registry_file));

# Get the species to use
my $species = shift;
die ("You must specify the species") unless (defined($species));

# Get the variation database foreign key relationships file
my $foreign_key_file = shift;
die ("You must supply the variation database foreign key file") unless (defined($foreign_key_file));

#ÊIf a file of variation ids is specified, use that
my $variation_id_file = shift;
warn ("Will use variation ids from $variation_id_file") if (defined($variation_id_file));

# Parse the foreign keys into a hash. Keep the tables that point to variation_id in variation
my %foreign_keys;
my @variation_tables;
open (KEY,"<",$foreign_key_file) or die ("Could not open foreign key file $foreign_key_file for parsing");
while (<KEY>) {
    chomp;
    
    my ($table,$col,$ref_table,$ref_col) = $_ =~ m/^ALTER\s+TABLE\s+(\S+)\s+ADD\s+FOREIGN\s+KEY\s+\((\S+)\)\s+REFERENCES\s+([^\s\(]+)\s*\(([^\s\)]+)\)/i;
    next unless (defined($table) && defined($col) && defined($ref_table) && defined($ref_col));
    
    $foreign_keys{$table} = [] unless (exists($foreign_keys{$table}));
    push(@{$foreign_keys{$table}},[$col,$ref_table,$ref_col]);
    
    push(@variation_tables,[$table,$col]) if ($ref_table eq 'variation' && $ref_col eq 'variation_id');
}
close (KEY);

# Load the registry
my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_all($registry_file);

# Get a db adaptor to the source database
my $source = $reg->get_DBAdaptor($species,'source_db') or die ("Could not get a DBAdaptor for the source db");

#ÊIf a variation id file was specified, parse the variation_ids to use from there.
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
    
    warn (localtime() . "\tWill extract $VARIATION_SIZE random variations from source database");
    
    # Get the min and max variation_ids
    my $stmt = qq{
        SELECT
            MIN(variation_id),
            MAX(variation_id)
        FROM
            variation
    };
    my ($min,$max) = @{$source->dbc->db_handle->selectall_arrayref($stmt)->[0]};
    
    # Randomly pull variation ids from the database until we've reached the desired number of variations
    $stmt = qq{
        SELECT
            variation_id
        FROM
            variation
        WHERE
            variation_id = ?
        LIMIT 1
    };
    my $sth = $source->dbc->prepare($stmt);
    my %ids;
    while (scalar(keys(%ids)) < $VARIATION_SIZE) {
        my $id = $min + int(rand($max-$min+1));
        $sth->execute($id);
        $sth->bind_columns(\$id);
        $sth->fetch();
        next unless defined($id);
        $ids{$id}++;
    }
    @variation_ids = keys(%ids);
     
}

# Populate the destination database

# Get a DBAdaptor to the destination database
my $dest = $reg->get_DBAdaptor($species,'dest_db') or die ("Could not get a DBAdaptor for the destination db");

#ÊFirst, populate the variation table by looping over the variation ids
my $table = 'variation';
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
        tbl.variation_id = ? 
};
my $ins_sth = $dest->dbc->prepare($ins_stmt);

map {$ins_sth->execute($_)} @variation_ids;
    
# First we'll import the "core" data, i.e. data from the tables having variation_id as a foreign key
foreach my $row (@variation_tables) {
    
    my $table = $row->[0];
    my $column = $row->[1];
    
    # Recursively add foreign relationships to this table
    add_foreign_data($dest,'variation','variation_id',$table,$column,$srcdb,\%foreign_keys);
    
}

# Lastly, process all the foreign key relationships without recursion
foreach my $table (keys(%foreign_keys)) {
    
    foreach my $row (@{$foreign_keys{$table}}) {
        
        add_foreign_data($dest,$table,$row->[0],$row->[1],$row->[2],$srcdb,{});
        
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
    
    #ÊGet the join condition
    my $constraint = qq{src.$foreign_column = dst.$column};
    
    #ÊCreate the insert statement
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

