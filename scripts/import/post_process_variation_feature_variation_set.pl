#! perl -w

=head

    This script will post-process the variation_feature table to update the variation_set_id column to contain the
    primary keys of the variation_sets which the corresponding variation is part of.
    It does so by fetching all variations for either a specified variation set (useful if only updating) or for
    all available variation sets. For now, we'll assume that the script and server can handle ALL variations
    belonging to a set simultaneously. If not, it may be necessary to split the query by e.g. seq_region or
    variation_id range. Note also that the script disables the keys. This may not be entirely necessary,
    it will probably suffice to just disable the index on the variation_set_id column...
    
=cut


use strict;
use warnings;
use Getopt::Long;

use Bio::EnsEMBL::Registry;

my $MAX_VARIATION_SETS = 64;
# This is just for convenience in case we have a custom name for the variation_feature table
my $VARIATION_FEATURE_TABLE = 'variation_feature';

my @option_defs = (
  'species=s',
  'registry_file=s',
  'chromosome=s@',
  'verbose!'
);

my %options;
GetOptions(\%options,@option_defs);

my $species = $options{'species'};
my $registry_file = $options{'registry_file'};

#ÊCheck that required parameters were passed
die ("Required argument '-species' was not specified") unless (defined($species));
die ("Required argument '-registry_file' was not specified") unless (defined($registry_file));

# Load the registry from the supplied file
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all($registry_file);

#ÊGet a dbadaptor to the variation database
my $dbVar = $registry->get_DBAdaptor($species,"variation") or die ("Could not get variation DBAdaptor for $species");

#ÊIf no chromosomes were specified, get all available in the database
my @chromosomes;
if (defined($options{'chromosome'})) {
    @chromosomes = @{$options{'chromosome'}};
}
else {
    my $stmt = qq{
        SELECT DISTINCT
            name
        FROM
            seq_region;
    };
    @chromosomes = @{$dbVar->dbc->db_handle->selectcol_arrayref($stmt)};
}

# This statement will set the variation_set_id column for the sets that are explicitly set in the variation_set_variation table
my $stmt = qq{
    UPDATE
        variation_feature vf USE INDEX (pos_idx) JOIN
        seq_region sr ON (
            sr.seq_region_id = vf.seq_region_id
        )
    SET
        vf.variation_set_id = (
            SELECT
                GROUP_CONCAT(vsv.variation_set_id)
            FROM
                variation_set_variation vsv USE INDEX (PRIMARY)
            WHERE
                vsv.variation_id = vf.variation_id
            GROUP BY
                vsv.variation_id
        )
    WHERE
        sr.name = ?
};
my $set_explicit_sets = $dbVar->dbc->prepare($stmt);

#ÊThis statment will update the variation_set_id column to also contain the parent set in case there is a subset structure
$stmt = qq{
    UPDATE
        variation_feature vf USE INDEX (pos_idx) JOIN
        seq_region sr ON (
            sr.seq_region_id = vf.seq_region_id
        ),
        variation_set_structure vss
    SET
        vf.variation_set_id = (vf.variation_set_id | POW(2,(vss.variation_set_super - 1)))
    WHERE
        sr.name = ? AND
        FIND_IN_SET(vss.variation_set_sub,vf.variation_set_id)
};
my $set_inherited_sets = $dbVar->dbc->prepare($stmt);

# Loop over the chromosomes
foreach my $chromosome (@chromosomes) {
    
    #ÊFirst, get the explicit sets
    my $updated = $set_explicit_sets->execute($chromosome);
    
    #ÊThen set the implicit sets by looping until no more rows are updated
    while ($updated) {
        
        $updated = $set_inherited_sets->execute($chromosome);
        
    }
}

#ÊGet a VariationSetAdaptor
my $vs_adaptor = $dbVar->get_VariationSetAdaptor() or die ("Could not get a VariationSetAdaptor from variation db");

#ÊGet a VariationFeatureAdaptor
my $vf_adaptor = $dbVar->get_VariationFeatureAdaptor() or die ("Could not get VariationFeatureAdaptor from variation db");

# If a variation_set was specified, get an object for it
my @variation_sets;
if (defined($options{'variation_set_id'})) {
    foreach my $variation_set_id (@{$options{'variation_set_id'}}) {
        my $variation_set = $vs_adaptor->fetch_by_dbID($variation_set_id) or warn ("Could not get a VariationSet for variation_set_id '$variation_set_id'");
        push(@variation_sets,$variation_set);
    }
}
elsif (defined($options{'variation_set_name'})) {
    foreach my $variation_set_name (@{$options{'variation_set_name'}}) {
        my $variation_set = $vs_adaptor->fetch_by_name($variation_set_name) or warn ("Could not get a VariationSet for variation_set_name '$variation_set_name'");
        push(@variation_sets,$variation_set);
    }
}
# Else, get all VariationSets
else {
    push(@variation_sets,@{$vs_adaptor->fetch_all()});
}

#ÊSince we just need the variation_feature_ids for variations belonging to a variation set, do a direct SQL query instead of going through the API
my $stmt = qq{
    SELECT
        vf.variation_feature_id
    FROM
        variation_set_variation vsv JOIN
        `$VARIATION_FEATURE_TABLE` vf ON (
            vf.variation_id = vsv.variation_id
        )
    WHERE
        vsv.variation_set_id = ?
};
my $sel_sth = $dbVar->dbc->prepare($stmt);

# Prepare a statement for updating the variation_feature table
$stmt = qq{
    UPDATE
        `$VARIATION_FEATURE_TABLE`
    SET 
        variation_set_id = (variation_set_id | ?)
    WHERE
        variation_feature_id = ?
    LIMIT 1
};
my $up_sth = $dbVar->dbc->prepare($stmt);

#print STDOUT localtime() . qq{\tDisabling keys on the variation_feature table...} if (defined($options{'verbose'}));

#ÊTemporarily disable the keys on the variation_feature to speed up the updates
$stmt = qq{
    ALTER TABLE
        `$VARIATION_FEATURE_TABLE`
    DISABLE KEYS
};
#$dbVar->dbc->do($stmt);

#print STDOUT qq{ Done!\n} if (defined($options{'verbose'}));

# Wrap the actual processing in an eval statement since if we fall over, we need to enable the keys on the table again
eval {
    
    #ÊLoop over the variation sets
    while (my $variation_set = shift(@variation_sets)) {
        
        # Skip any undefs
        next unless (defined($variation_set));
        
        # If the variation_set_id is above 64, we cannot store it in the set. Either we have hit the limit for what we can actually store or the variation_set primary key column is sparsely populated. Anyway, can't add this set to the VariationFeature so warn and skip
        if ($variation_set->dbID() > $MAX_VARIATION_SETS) {
            warn ("The VariationSet '" . $variation_set->name() . "' has primary key " . $variation_set->dbID() . "! This number is too large to store in the VariationFeature table." .
                  "If possible, review the primary keys in the variation_set table (and related tables) and try to assign a lower key value." .
                  "If not possible, we cannot store any more VariationSets with this construct!");
            next;
        }
        
        # Calculate the bit value for this VariationFeature's set membership
        my $bitvalue = (2 ** ($variation_set->dbID() - 1));
        
        print STDOUT localtime() . "\tProcessing variation features belonging to VariationSet '" . $variation_set->name() . "', bitvalue for this set (dbID = " . $variation_set->dbID() . ") is $bitvalue...\n" if (defined($options{'verbose'}));
        
        # Execute the SQL to get the variation features
        $sel_sth->execute($variation_set->dbID());
        
        # Loop over the results
        while (my $vf = $sel_sth->fetchrow_arrayref()) {
        
            # Update the variation_set_id set
            $up_sth->execute($bitvalue,$vf->[0]);
                
        }
        
        print STDOUT qq{Done!\n} if (defined($options{'verbose'}));
    }
}; warn ($@) if ($@);

#ÊEnable the keys again
#print STDOUT localtime() . qq{\tEnabling keys on the variation_feature table...} if (defined($options{'verbose'}));
$stmt = qq{
    ALTER TABLE
        `$VARIATION_FEATURE_TABLE`
    ENABLE KEYS
};
#$dbVar->dbc->do($stmt);

#print STDOUT qq{ Done!\n} if (defined($options{'verbose'}));

