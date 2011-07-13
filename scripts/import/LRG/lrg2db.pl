#!perl -w

use strict;

use Getopt::Long;
use List::Util qw (min max);
use LRG::LRG;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use DBI qw(:sql_types);

my $xmlfile;
my $host;
my $port;
my $user;
my $pass;
my $dbname;
my $verbose;
my $hgnc_symbol;
my $purge;
my $lrg_id;
my $gene_id;
my $keep_fixed;
my $keep_mapping;
my $keep_updatable;
my $keep_request;
my @update_annotation_set;

GetOptions(
  'host=s'		=> \$host,
  'port=i'		=> \$port,
  'dbname=s'		=> \$dbname,
  'user=s'		=> \$user,
  'pass=s'		=> \$pass,
  'xmlfile=s'		=> \$xmlfile,
  'verbose!'		=> \$verbose,
  'hgnc_symbol=s'       => \$hgnc_symbol,
  'lrg_id=s'            => \$lrg_id,
  'purge!'              => \$purge,
  'keep_fixed!'         => \$keep_fixed,
  'keep_mapping!'       => \$keep_mapping,
  'keep_updatable!'     => \$keep_updatable,
  'keep_request!'       => \$keep_request,
  'replace_updatable=s'   => \@update_annotation_set
);

die("Database credentials (-host, -port, -dbname, -user) need to be specified!") unless (defined($host) && defined($port) && defined($dbname) && defined($user));
die("An input LRG XML file must be specified") unless (defined($xmlfile) || defined($purge));
die("A correspondiong HGNC symbol must be specified") unless (defined($hgnc_symbol) || (defined($purge) && defined($lrg_id)));
die("An LRG id must be specified") unless (!defined($purge) || defined($lrg_id) || defined($hgnc_symbol));

# Get a database connection
print STDOUT localtime() . "\tConnecting to database $dbname\n" if ($verbose);
my $db_adaptor = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host => $host,
  -user => $user,
  -pass => $pass,
  -port => $port,
  -dbname => $dbname
) or die("Could not get a database adaptor for $dbname on $host:$port");
print STDOUT localtime() . "\tConnected to $dbname on $host:$port\n" if ($verbose);

# If the database should be purged of all LRG XML-related data, do that
if ($purge) {
    
    # Get the database gene_id for the HGNC symbol and LRG id
    my $stmt = qq{
        SELECT
            gene_id,
            symbol,
            lrg_id
        FROM
            gene
        WHERE
    };
    if (defined($lrg_id)) {
        $stmt .= qq{ lrg_id = '$lrg_id'};
    }
    else {
        $stmt .= qq{ symbol = '$hgnc_symbol'};
    }
    $stmt .= qq{ LIMIT 1};
    ($gene_id,$hgnc_symbol,$lrg_id) = @{$db_adaptor->dbc->db_handle->selectall_arrayref($stmt)->[0]} or die ("No gene could be found with " . (defined($hgnc_symbol) ? "HGNC symbol $hgnc_symbol" : "LRG id $lrg_id"));
    
    # Delete the mapping data
    $stmt = qq{
        DELETE FROM
            lasm,
            lm,
            lms,
            lmd
        USING
            lrg_annotation_set las JOIN (
                lrg_annotation_set_mapping lasm LEFT JOIN (
                    lrg_mapping lm LEFT JOIN
                    lrg_mapping_span lms USING (mapping_id) LEFT JOIN
                    lrg_mapping_diff lmd USING (mapping_span_id)
                ) USING (mapping_id) 
            ) USING (annotation_set_id)
        WHERE
            las.gene_id = $gene_id
    };
    if (!$keep_mapping) {
        print STDOUT localtime() . "\tRemoving mapping from updatable annotation for $lrg_id\n" if ($verbose);
        $db_adaptor->dbc->do($stmt);
    }
    else {
        print STDOUT localtime() . "\tKeeping mapping in updatable annotation for $lrg_id\n" if ($verbose);
    }
    
    # Delete the updatable annotation data associated with this LRG. First, just delete the XML column. If the annotation set is not linked to any mapping, delete it altogether.
    if (!$keep_updatable) {
        $stmt = qq{
            UPDATE
                lrg_annotation_set las
            SET
                las.xml = NULL
            WHERE
                gene_id = $gene_id
        };
        print STDOUT localtime() . "\tRemoving annotation sets from updatable annotation for $lrg_id\n" if ($verbose);
        $db_adaptor->dbc->do($stmt);
        
        $stmt = qq{
        DELETE FROM
            las
        USING
            lrg_annotation_set las
        WHERE
            las.gene_id = $gene_id AND
            NOT EXISTS (
                SELECT
                    *
                FROM
                    lrg_annotation_set_mapping lasm
                WHERE
                    lasm.annotation_set_id = las.annotation_set_id
            )
        };
        $db_adaptor->dbc->do($stmt);
    }
    else {
        print STDOUT localtime() . "\tKeeping annotation sets in updatable annotation for $lrg_id\n" if ($verbose);
    }
    
    # Delete the LRG request entries
    $stmt = qq{
        DELETE FROM
            lrg_request
        WHERE
            gene_id = $gene_id
    };
    if (!$keep_request) {
        print STDOUT localtime() . "\tRemoving requester information for $lrg_id\n" if ($verbose);
        $db_adaptor->dbc->do($stmt);
    }
    else {
        print STDOUT localtime() . "\tKeeping LRG request data for $lrg_id\n" if ($verbose);
    }
    
    # Delete the fixed annotation data
    $stmt = qq{
        DELETE FROM
            ld,
            lt,
            lc,
            lp,
            lgs,
            lcs,
            lps,
            ls,
            le,
            li,
            lcc
        USING
            lrg_data ld LEFT JOIN
            lrg_transcript lt USING (gene_id) LEFT JOIN
            lrg_cdna lc USING (transcript_id) LEFT JOIN
            lrg_cds lp USING (transcript_id) LEFT JOIN
            lrg_genomic_sequence lgs USING (gene_id) LEFT JOIN
            lrg_cdna_sequence lcs USING (cdna_id) LEFT JOIN
            lrg_peptide_sequence lps USING (cds_id) LEFT JOIN
            lrg_sequence ls ON (lps.sequence_id = ls.sequence_id OR lgs.sequence_id = ls.sequence_id OR lcs.sequence_id = ls.sequence_id) LEFT JOIN
            lrg_exon le USING (transcript_id) LEFT JOIN
            lrg_intron li ON li.exon_5 = le.exon_id LEFT JOIN
            lrg_cds_codon lcc USING (cds_id)
        WHERE
            ld.gene_id = $gene_id
    };
    if (!$keep_fixed) {
        print STDOUT localtime() . "\tRemoving entries from fixed annotation for $lrg_id\n" if ($verbose);
        $db_adaptor->dbc->do($stmt);
    }
    else {
        print STDOUT localtime() . "\tKeeping fixed annotation for $lrg_id\n" if ($verbose);
    }
    
    print STDOUT localtime() . "\tDone removing $lrg_id from database\n";
    exit;
}

print STDOUT localtime() . "\tCreating LRG object from input XML file $xmlfile\n" if ($verbose);
my $lrg = LRG::LRG::newFromFile($xmlfile) or die("ERROR: Could not create LRG object from XML file!");

# Get the fixed section node
my $fixed = $lrg->findNode('fixed_annotation') or die ("ERROR: Could not find fixed annotation section in LRG file");

# Get the LRG id node
my $node = $fixed->findNode('id') or die ("ERROR: Could not find LRG identifier");
$lrg_id = $node->content();

# Get the organism and taxon_id
$node = $fixed->findNode('organism') or die ("ERROR: Could not find organism tag");
my $taxon_id = $node->data()->{'taxon'};
my $organism = $node->content();

# Get the moltype
$node = $fixed->findNode('mol_type') or die ("ERROR: Could not find moltype tag");
my $moltype = $node->content();

# Get the creation date
$node = $fixed->findNode('creation_date') or die ("ERROR: Could not find creation date tag");
my $creation_date = $node->content();

# Get the LRG sequence
$node = $fixed->findNode('sequence') or die ("ERROR: Could not find LRG sequence tag");
my $lrg_seq = $node->content();

# Get the database gene_id for the HGNC symbol and LRG id
my $stmt = qq{
    SELECT
        gene_id
    FROM
        gene
    WHERE
        symbol = '$hgnc_symbol' AND
        (
            lrg_id IS NULL OR
            lrg_id = '$lrg_id'
        )
    LIMIT 1
};
$gene_id = $db_adaptor->dbc->db_handle->selectall_arrayref($stmt)->[0][0] or die ("No gene could be find with HGNC symbol $hgnc_symbol and LRG id $lrg_id");

# If we are updating the annotation sets, do that here
if (scalar(@update_annotation_set)) {
        
    # Parse the updatable section to get the annotation sets
    my $updatable = $lrg->findNode('updatable_annotation') or die ("ERROR: Could not find updatable annotation section in LRG file $xmlfile");
    my $annotation_sets = $updatable->findNodeArray('annotation_set');
    $annotation_sets ||= [];
    
    while (my $annotation_set = shift(@{$annotation_sets})) {
    
        parse_annotation_set($annotation_set,$gene_id,$db_adaptor,\@update_annotation_set);
        
    }
    
    exit(0);
}

# Set the LRG id for the gene
lrg_id($gene_id,$db_adaptor,$lrg_id);

# Insert the data into the db
$stmt = qq{
    INSERT INTO
        lrg_data (
            gene_id,
            organism,
            taxon_id,
            moltype,
            creation_date
        )
    VALUES (
        '$gene_id',
        '$organism',
        $taxon_id,
        '$moltype',
        '$creation_date'
    )
};
print STDOUT localtime() . "\tAdding LRG data for $lrg_id to database\n" if ($verbose);
$db_adaptor->dbc->do($stmt);

# Insert the sequence
add_sequence($gene_id,'genomic',$db_adaptor,$lrg_seq);

# Some usefule prepared statements
my $tr_ins_stmt = qq{
    INSERT INTO
        lrg_transcript (
            gene_id,
            transcript_name
        )
    VALUES (
        '$gene_id',
        ?
    )
};
my $cdna_ins_stmt = qq{
    INSERT INTO
        lrg_cdna (
            transcript_id,
            lrg_start,
            lrg_end
        )
    VALUES (
        ?,
        ?,
        ?
    )
};
my $cds_ins_stmt = qq{
    INSERT INTO
        lrg_cds (
            transcript_id,
            lrg_start,
            lrg_end,
            codon_start
        )
    VALUES (
        ?,
        ?,
        ?,
        ?
    )
};
my $cc_ins_stmt = qq{
    INSERT INTO
        lrg_cds_codon (
            cds_id,
            type,
            codon
        )
    VALUES (
        ?,
        ?,
        ?
    )
};
my $exon_ins_stmt = qq{
    INSERT INTO
        lrg_exon (
            transcript_id,
            lrg_start,
            lrg_end,
            cdna_start,
            cdna_end,
            peptide_start,
            peptide_end
        )
    VALUES (
        ?,
        ?,
        ?,
        ?,
        ?,
        ?,
        ?
    )
};
my $intron_ins_stmt = qq{
    INSERT INTO
        lrg_intron (
            exon_5,
            exon_3,
            phase
        )
    VALUES (
        ?,
        ?,
        ?
    )
};
my $tr_ins_sth = $db_adaptor->dbc->prepare($tr_ins_stmt);
my $cdna_ins_sth = $db_adaptor->dbc->prepare($cdna_ins_stmt);
my $cc_ins_sth = $db_adaptor->dbc->prepare($cc_ins_stmt);
my $cds_ins_sth = $db_adaptor->dbc->prepare($cds_ins_stmt);
my $exon_ins_sth = $db_adaptor->dbc->prepare($exon_ins_stmt);
my $intron_ins_sth = $db_adaptor->dbc->prepare($intron_ins_stmt);

# Get the requester data
$node = $fixed->findNodeArray('source') or warn ("Could not find requester data");
if (defined($node)) {
    while (my $source = shift(@{$node})) {
        parse_source($source,$gene_id,$db_adaptor);
    }
}

# Get the transcript nodes
my $transcripts = $fixed->findNodeArray('transcript') or die ("ERROR: Could not find transcript tags");

# Parse and add each transcript to the database separately
while (my $transcript = shift(@{$transcripts})) {
    
    # Transcript name and LRG coords
    my $name = $transcript->data()->{'name'};
    my $lrg_start = $transcript->data()->{'start'};
    my $lrg_end = $transcript->data()->{'end'};
    
    # Insert the transcript into db
    $tr_ins_sth->bind_param(1,$name,SQL_VARCHAR);
    $tr_ins_sth->execute();
    my $transcript_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};
    
    # Get the cdna
    $node = $transcript->findNode('cdna/sequence') or warn("Could not get cdna sequence for transcript $name\, skipping this transcript");
    next unless ($node);
    my $seq = $node->content();
    
    # Insert the cdna into db
    $cdna_ins_sth->bind_param(1,$transcript_id,SQL_INTEGER);
    $cdna_ins_sth->bind_param(2,$lrg_start,SQL_INTEGER);
    $cdna_ins_sth->bind_param(3,$lrg_end,SQL_INTEGER);
    $cdna_ins_sth->execute();
    my $cdna_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};
    
    # Insert the cdna sequence
    add_sequence($cdna_id,'cdna',$db_adaptor,$seq);
    
    # Get the cds
    my $cds_node = $transcript->findNode('coding_region') or warn("Could not get coding region for transcript $name\, skipping this transcript");
    next unless ($cds_node);
    $lrg_start = $cds_node->data()->{'start'};
    $lrg_end = $cds_node->data()->{'end'};
    my $codon_start = $cds_node->data()->{'codon_start'};
    
    $node = $cds_node->findNode('translation/sequence') or warn("Could not get translated sequence for transcript $name\, skipping this transcript");
    next unless ($node);
    $seq = $node->content();
    
    # Insert the cds into db
    $cds_ins_sth->bind_param(1,$transcript_id,SQL_INTEGER);
    $cds_ins_sth->bind_param(2,$lrg_start,SQL_INTEGER);
    $cds_ins_sth->bind_param(3,$lrg_end,SQL_INTEGER);
    $cds_ins_sth->bind_param(4,$codon_start,SQL_INTEGER);
    $cds_ins_sth->execute();
    my $cds_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};
    
    # Insert the peptide sequence
    add_sequence($cds_id,'peptide',$db_adaptor,$seq);
    
    # Get any non-standard codons
    $node = $cds_node->findNodeArray('selenocysteine');
    $node ||= [];
    push(@{$node},@{($cds_node->findNodeArray('pyrrolysine') or [])});
    while (my $c = shift(@{$node})) {
        my $type = $c->name();
        my $codon = $c->data()->{'codon'};
        $cc_ins_sth->bind_param(1,$cds_id,SQL_INTEGER);
        $cc_ins_sth->bind_param(2,$type,SQL_VARCHAR);
        $cc_ins_sth->bind_param(3,$codon,SQL_INTEGER);
        $cc_ins_sth->execute();
    }
    
    # Get the children nodes of the transcript and iterate over the exons and introns
    my $children = $transcript->getAllNodes();
    my $phase;
    my $last_exon;
    while (my $child = shift(@{$children})) {
        # Skip if it's not an intron or exon
        next if ($child->name() ne 'exon' && $child->name() ne 'intron');
        
        # If we have an exon, parse out the data
        if ($child->name() eq 'exon') {
            my $exon_lrg_start;
            my $exon_lrg_end;
            
            my $cdna_start;
            my $cdna_end;
            
            my $peptide_start;
            my $peptide_end;
            
            $node = $child->findNode('lrg_coords') or warn("Could not get LRG coordinates for one or more exons in $name");
            if ($node) {
                $exon_lrg_start = $node->data()->{'start'};
                $exon_lrg_end = $node->data()->{'end'};
            }
            $node = $child->findNode('cdna_coords') or warn("Could not get cDNA coordinates for one or more exons in $name"); 
            if ($node) {
                $cdna_start = $node->data()->{'start'};
                $cdna_end = $node->data()->{'end'};
            }
            $node = $child->findNode('peptide_coords'); 
            if ($node) {
                $peptide_start = $node->data()->{'start'};
                $peptide_end = $node->data()->{'end'};
            }
            
            # Insert the exon into db
            $exon_ins_sth->bind_param(1,$transcript_id,SQL_INTEGER);
            $exon_ins_sth->bind_param(2,$exon_lrg_start,SQL_INTEGER);
            $exon_ins_sth->bind_param(3,$exon_lrg_end,SQL_INTEGER);
            $exon_ins_sth->bind_param(4,$cdna_start,SQL_INTEGER);
            $exon_ins_sth->bind_param(5,$cdna_end,SQL_INTEGER);
            $exon_ins_sth->bind_param(6,$peptide_start,SQL_INTEGER);
            $exon_ins_sth->bind_param(7,$peptide_end,SQL_INTEGER);
            $exon_ins_sth->execute();
            my $exon_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};
            
            # If an intron was preceeding this exon, we should insert that one as well
            if (defined($phase)) {
                $intron_ins_sth->bind_param(1,$last_exon,SQL_INTEGER);
                $intron_ins_sth->bind_param(2,$exon_id,SQL_INTEGER);
                $intron_ins_sth->bind_param(3,$phase,SQL_INTEGER);
                $intron_ins_sth->execute();
                # Unset phase so that it will be set for next intron
                undef($phase);
            }
            # Store the current exon_id as last_exon to be used as upstream for the next intron
            $last_exon = $exon_id;
        }
        # Else, this is an intron so get the phase
        elsif (exists($child->data()->{'phase'})) {
            $phase = $child->data()->{'phase'};
        }
    }
}



# Parse the updatable section to get the annotation sets
my $updatable = $lrg->findNode('updatable_annotation') or die ("ERROR: Could not find updatable annotation section in LRG file $xmlfile");
my $annotation_sets = $updatable->findNodeArray('annotation_set');
$annotation_sets ||= [];

while (my $annotation_set = shift(@{$annotation_sets})) {

    parse_annotation_set($annotation_set,$gene_id,$db_adaptor,\@update_annotation_set);
    
}

sub parse_annotation_set {
    my $annotation_set = shift;
    my $gene_id = shift;
    my $db_adaptor = shift;
    my $use_annotation_set = shift || [];
        
    # Use different syntax depending on whether we are inserting or updating
    my $insert_mode = (scalar(@{$use_annotation_set}) ? qq{REPLACE} : qq{INSERT IGNORE});
    my $as_ins_stmt = qq{
        $insert_mode INTO
            lrg_annotation_set (
                gene_id,
                source,
                comment,
                modification_date,
                lrg_gene_name,
                xml
            )
        VALUES (
            '$gene_id',
            ?,
            ?,
            ?,
            ?,
            ?
        )
    };
    my $asm_ins_stmt = qq{
        $insert_mode INTO
            lrg_annotation_set_mapping (
                annotation_set_id,
                mapping_id
            )
        VALUES (
            ?,
            ?
        )
    };
    my $as_ins_sth = $db_adaptor->dbc->prepare($as_ins_stmt);
    my $asm_ins_sth = $db_adaptor->dbc->prepare($asm_ins_stmt);

    # Get and parse the source
    my $source = $annotation_set->findNode('source') or die ("Could not find any source for annotation_set in $xmlfile");
    my $source_id = parse_source($source,$gene_id,$db_adaptor,$use_annotation_set) or warn ("Could not properly parse source information in annotation_set in $xmlfile");
    return $source_id if (!defined($source_id) || $source_id < 0);
    
    # Get comment and make sure it belongs to the annotation set element
    my $comment = $annotation_set->findNode('comment');
    if (defined($comment) && $comment->parent()->name() eq 'annotation_set') {
        $comment = $comment->content();
    }
    else {
        undef($comment);
    }
    
    # Get modification_date
    my $modification_date = $annotation_set->findNode('modification_date') or die ("Could not find modification date for annotation_set in $xmlfile");
    $modification_date = $modification_date->content();
    
    # Get the mapping from the updatable section
    my $mappings = $annotation_set->findNodeArray('mapping');
    my @mids;
    $mappings ||= [];
    while (my $mapping = shift(@{$mappings})) {
        push(@mids,parse_mapping($mapping,$db_adaptor,$gene_id,$xmlfile));
    }
    
    # Get the lrg_gene_name from the updatable section
    my $lrg_gene_name = $annotation_set->findNode('lrg_gene_name');
    if (defined($lrg_gene_name)) {
        $lrg_gene_name = $lrg_gene_name->content();
        warn("HGNC symbol for $lrg_id as specified in XML file is different from corresponding symbol stored in db ($lrg_gene_name in xml vs $hgnc_symbol in db)") if ($lrg_gene_name ne $hgnc_symbol);
    }
    
    # Create an XML writer to print the XML of the rest of the nodes
    my $xml_string;
    my $xml_out;
    
    # Get the xml for the rest of the annotation_set
    foreach my $name (('other_exon_naming','alternate_amino_acid_numbering','features')) {
        my $node = $annotation_set->findNode($name);
        if (defined($node)) {
            my $xml = new XML::Writer(OUTPUT => \$xml_out, DATA_INDENT => 2, DATA_MODE => 1);
            $node->xml($xml);
            $node->printNode();
            $xml_out .= "\n";
        }
    }
    
    $as_ins_sth->bind_param(1,$source_id,SQL_INTEGER);
    $as_ins_sth->bind_param(2,$comment,SQL_VARCHAR);
    $as_ins_sth->bind_param(3,$modification_date,SQL_VARCHAR);
    $as_ins_sth->bind_param(4,$lrg_gene_name,SQL_VARCHAR);
    $as_ins_sth->bind_param(5,$xml_out,SQL_VARCHAR);
    $as_ins_sth->execute();
    my $annotation_set_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};
    
    # Link the mappings to the annotation_set
    $asm_ins_sth->bind_param(1,$annotation_set_id,SQL_INTEGER);
    while (my $mid = shift(@mids)) {
        $asm_ins_sth->bind_param(2,$mid,SQL_INTEGER);
        $asm_ins_sth->execute();
    }
    
    return $annotation_set_id;
}

sub parse_mapping {
    my $mapping = shift;
    my $db_adaptor = shift;
    my $gene_id = shift;
    my $xmlfile = shift;

    my $m_ins_stmt = qq{
        INSERT INTO
            lrg_mapping (
                gene_id,
                assembly,
                chr_name,
                chr_id,
                chr_start,
                chr_end,
                most_recent
            )
        VALUES (
            '$gene_id',
            ?,
            ?,
            ?,
            ?,
            ?,
            ?
        )
    };
    my $ms_ins_stmt = qq{
        INSERT INTO
            lrg_mapping_span (
                mapping_id,
                lrg_start,
                lrg_end,
                chr_start,
                chr_end,
                strand
            )
        VALUES (
            ?,
            ?,
            ?,
            ?,
            ?,
            ?
        )
    };
    my $diff_ins_stmt = qq{
        INSERT INTO
            lrg_mapping_diff (
                mapping_span_id,
                lrg_start,
                lrg_end,
                chr_start,
                chr_end,
                type,
                lrg_sequence,
                chr_sequence
            )
        VALUES (
            ?,
            ?,
            ?,
            ?,
            ?,
            ?,
            ?,
            ?
        )
    };
    my $m_ins_sth = $db_adaptor->dbc->prepare($m_ins_stmt);
    my $ms_ins_sth = $db_adaptor->dbc->prepare($ms_ins_stmt);
    my $diff_ins_sth = $db_adaptor->dbc->prepare($diff_ins_stmt);
    
    my $assembly = $mapping->data()->{'assembly'} or die ("Could not find assembly attribute in mapping tag of $xmlfile");
    my $chr_name = $mapping->data()->{'chr_name'} or die ("Could not find chr_name attribute in mapping tag of $xmlfile");
    my $chr_start = $mapping->data()->{'chr_start'} or die ("Could not find chr_start attribute in mapping tag of $xmlfile");
    my $chr_end = $mapping->data()->{'chr_end'} or die ("Could not find chr_end attribute in mapping tag of $xmlfile");
    my $chr_id = $mapping->data()->{'chr_id'};
    my $most_recent = $mapping->data()->{'most_recent'};
    my $lrg_start;
    my $lrg_end;
    my $strand;
    
    my $spans = $mapping->findNodeArray('mapping_span') or die ("Could not find any mapping spans in mapping tag of $xmlfile");
    
    $m_ins_sth->bind_param(1,$assembly,SQL_VARCHAR);
    $m_ins_sth->bind_param(2,$chr_name,SQL_VARCHAR);
    $m_ins_sth->bind_param(3,$chr_id,SQL_VARCHAR);
    $m_ins_sth->bind_param(4,$chr_start,SQL_INTEGER);
    $m_ins_sth->bind_param(5,$chr_end,SQL_INTEGER);
    $m_ins_sth->bind_param(6,(defined($most_recent) ? 1 : 0),SQL_INTEGER);
    $m_ins_sth->execute();
    my $mapping_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};
    
    while (my $span = shift(@{$spans})) {
        $lrg_start = $span->data()->{'lrg_start'} or die ("Could not find lrg_start attribute in mapping span of $xmlfile");
        $lrg_end = $span->data()->{'lrg_end'} or die ("Could not find lrg_end attribute in mapping span of $xmlfile");
        $chr_start = $span->data()->{'start'} or die ("Could not find start attribute in mapping span of $xmlfile");
        $chr_end = $span->data()->{'end'} or die ("Could not find end attribute in mapping span of $xmlfile");
        $strand = $span->data()->{'strand'} or die ("Could not find strand attribute in mapping span of $xmlfile");
        
#        #$ms_ins_sth->bind_param(1,$mapping_id,SQL_INTEGER);
#        #$ms_ins_sth->bind_param(2,$lrg_start,SQL_INTEGER);
#        #$ms_ins_sth->bind_param(3,$lrg_end,SQL_INTEGER);
#        #$ms_ins_sth->bind_param(4,$chr_start,SQL_INTEGER);
#        #$ms_ins_sth->bind_param(5,$chr_end,SQL_INTEGER);
#        #$ms_ins_sth->bind_param(6,$strand,SQL_INTEGER);
        $ms_ins_sth->execute($mapping_id,$lrg_start,$lrg_end,$chr_start,$chr_end,$strand);
        my $span_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};
        
        my $diffs = $span->findNodeArray('diff');
        $diffs ||= [];
        while (my $diff = shift(@{$diffs})) {
            my $diff_type = $diff->data()->{'type'} or die ("Could not find type attribute in diff of $xmlfile");
            $lrg_start = $diff->data()->{'lrg_start'} or die ("Could not find lrg_start attribute in diff of $xmlfile");
            $lrg_end = $diff->data()->{'lrg_end'} or die ("Could not find lrg_end attribute in diff of $xmlfile");
            $chr_start = $diff->data()->{'start'} or die ("Could not find start attribute in diff of $xmlfile");
            $chr_end = $diff->data()->{'end'} or die ("Could not find end attribute in diff of $xmlfile");
            my $lrg_seq = $diff->data()->{'lrg_sequence'};
            my $chr_seq = $diff->data()->{'genomic_sequence'};
            
            $diff_ins_sth->bind_param(1,$span_id,SQL_INTEGER);
            $diff_ins_sth->bind_param(2,$lrg_start,SQL_INTEGER);
            $diff_ins_sth->bind_param(3,$lrg_end,SQL_INTEGER);
            $diff_ins_sth->bind_param(4,$chr_start,SQL_INTEGER);
            $diff_ins_sth->bind_param(5,$chr_end,SQL_INTEGER);
            $diff_ins_sth->bind_param(6,$diff_type,SQL_VARCHAR);
            $diff_ins_sth->bind_param(7,$lrg_seq,SQL_VARCHAR);
            $diff_ins_sth->bind_param(8,$chr_seq,SQL_VARCHAR);
            $diff_ins_sth->execute();
        }
    }
    
    return $mapping_id;
}

sub parse_source {
    my $source = shift;
    my $gene_id = shift;
    my $db_adaptor = shift;
    my $use_annotation_set = shift || [];
    
    my $lr_ins_stmt = qq{
        INSERT IGNORE INTO
            lrg_request (
                gene_id,
                lsdb_id
            )
        VALUES (
            $gene_id,
            ?
        )
    };
    my $lsdb_ins_stmt = qq{
        INSERT INTO
            lsdb (
                name,
                url
            )
        VALUES (
            ?,
            ?
        )
    };
    my $contact_ins_stmt = qq{
        INSERT INTO
            contact (
                name,
                email,
                url,
                address
            )
        VALUES (
            ?,
            ?,
            ?,
            ?
        )
    };
    my $lc_ins_stmt = qq{
        INSERT IGNORE INTO
            lsdb_contact (
                lsdb_id,
                contact_id
            )
        VALUES (
            ?,
            ?
        )
    };
    my $lg_ins_stmt = qq{
        INSERT IGNORE INTO
            lsdb_gene (
                lsdb_id,
                gene_id
            )
        VALUES (
            ?,
            '$gene_id'
        )
    };
    my $lr_ins_sth = $db_adaptor->dbc->prepare($lr_ins_stmt);
    my $lsdb_ins_sth = $db_adaptor->dbc->prepare($lsdb_ins_stmt);
    my $contact_ins_sth = $db_adaptor->dbc->prepare($contact_ins_stmt);
    my $lc_ins_sth = $db_adaptor->dbc->prepare($lc_ins_stmt);
    my $lg_ins_sth = $db_adaptor->dbc->prepare($lg_ins_stmt);

    my $lsdb_name = $source->findNode('name');
    
    # Check that the parent node is source so that we are not in the contact section
    if (defined($lsdb_name) && $lsdb_name->parent()->name() eq 'source') {
        $lsdb_name = $lsdb_name->content();
        $lsdb_name ||= "";
    }
    else {
        undef($lsdb_name);
    }
    warn ("Could not find source name") if (!defined($lsdb_name));
    # Check if a specific source was asked for and return if this is not it
    return -1 unless (scalar(@{$use_annotation_set}) == 0 || (grep {$_ =~ m/^$lsdb_name$/i} @{$use_annotation_set}));
    
    my $lsdb_url = $source->findNode('url');
    # Check that the parent node is source so that we are not in the contact section
    if (defined($lsdb_url) && $lsdb_url->parent()->name() eq 'source') {
        $lsdb_url = $lsdb_url->content();
    }
    else {
        undef($lsdb_url);
    }
    
    # If we could not get the LSDB data from the source element, or if both name and url are blank, return undef
    $lsdb_url ||= "";
    return undef if (!defined($lsdb_name) || ($lsdb_name eq "" && $lsdb_url eq ""));
    
    # Enter the LSDB info into db
    my $stmt = qq{
        SELECT
            lsdb_id
        FROM
            lsdb
        WHERE
            name = '$lsdb_name' AND
            url = '$lsdb_url'
        LIMIT 1
    };
    my $lsdb_id = $db_adaptor->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
    if (!defined($lsdb_id)) {
        $lsdb_ins_sth->bind_param(1,$lsdb_name,SQL_VARCHAR);
        $lsdb_ins_sth->bind_param(2,$lsdb_url,SQL_VARCHAR);
        $lsdb_ins_sth->execute();
        $lsdb_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};
    }
    
    # Get the contact information for this source
    my $contacts = $source->findNodeArray('contact') or warn ("Could not find contact information for source " . (defined($lsdb_name) ? $lsdb_name : ""));
    $contacts ||= [];
    while (my $contact = shift(@{$contacts})) {
        my $name = $contact->findNode('name');
        $name = $name->content() if (defined($name));
        my $email = $contact->findNode('email'); 
        $email = $email->content() if (defined($email));
        my $url = $contact->findNode('url'); 
        $url = $url->content() if (defined($url));
        my $address = $contact->findNode('address'); 
        $address = $address->content() if (defined($address));
        
        # Enter the contact information into the db
        $stmt = qq{
            SELECT
                contact_id
            FROM
                contact
            WHERE
                name = '$name'
            LIMIT 1
        };
        my $contact_id = $db_adaptor->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
        if (!defined($contact_id)) {
            $contact_ins_sth->bind_param(1,$name,SQL_VARCHAR);
            $contact_ins_sth->bind_param(2,$email,SQL_VARCHAR);
            $contact_ins_sth->bind_param(3,$url,SQL_VARCHAR);
            $contact_ins_sth->bind_param(4,$address,SQL_VARCHAR);
            $contact_ins_sth->execute();
            $contact_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};
        }
        
        # Link the contact to the LSDB
        $lc_ins_sth->bind_param(1,$lsdb_id,SQL_INTEGER);
        $lc_ins_sth->bind_param(2,$contact_id,SQL_INTEGER);
        $lc_ins_sth->execute();
    }
    
    # Set the LSDB as a requester for this LRG if necessary (that is, if this source is in the fixed section)
    if ($source->parent()->name() eq 'fixed_annotation') {
        $lr_ins_sth->bind_param(1,$lsdb_id,SQL_INTEGER);
        $lr_ins_sth->execute();
    }
        
    # If necessary, link this LSDB to the gene
    $lg_ins_sth->bind_param(1,$lsdb_id,SQL_INTEGER);
    $lg_ins_sth->execute();
    
    return $lsdb_id;
}

sub add_sequence {
    my $id = shift;
    my $type = shift;
    my $db_adaptor = shift;
    my $sequence = shift;
    
    my $FIELD_MAX_LENGTH = 65535;
    my @sids;
    
    my $stmt = qq{
        INSERT INTO
            lrg_sequence (
                sequence
            )
        VALUES (
            ?
        )
    };
    my $sth = $db_adaptor->dbc->prepare($stmt);
    
    while (length($sequence) > 0) {
        my $substr = substr($sequence,0,$FIELD_MAX_LENGTH);
        $sth->bind_param(1,$substr,SQL_VARCHAR);
        $sth->execute();
        push(@sids,$db_adaptor->dbc->db_handle->{'mysql_insertid'});
        $sequence = substr($sequence,min(length($sequence),$FIELD_MAX_LENGTH));
    }
    
    if ($type =~ m/^genomic$/i) {
        $stmt = qq{
            INSERT IGNORE INTO
                lrg_genomic_sequence (
                    gene_id,
                    sequence_id
                )
            VALUES (
                '$id',
                ?
            )
        };
    }
    elsif ($type =~ m/^cdna$/i) {
        $stmt = qq{
            INSERT IGNORE INTO
                lrg_cdna_sequence (
                    cdna_id,
                    sequence_id
                )
            VALUES (
                $id,
                ?
            )
        };
    }
    elsif ($type =~ m/^peptide$/i) {
        $stmt = qq{
            INSERT IGNORE INTO
                lrg_peptide_sequence (
                    cds_id,
                    sequence_id
                )
            VALUES (
                $id,
                ?
            )
        };
    }
    else {
        warn("Unknown sequence type '$type' specified");
        return;
    }
    $sth = $db_adaptor->dbc->prepare($stmt);
    
    while (my $sid = shift(@sids)) {
        $sth->bind_param(1,$sid,SQL_INTEGER);
        $sth->execute();
    }
}

# Get/Set the LRG id for a gene_id
sub lrg_id {
    my $gene_id = shift;
    my $db_adaptor = shift;
    my $lrg_id = shift;
    
    if (defined($lrg_id)) {
        my $stmt = qq{
            UPDATE
                gene
            SET
                lrg_id = '$lrg_id'
            WHERE
                gene_id = $gene_id AND
                lrg_id IS NULL
        };
        $db_adaptor->dbc->do($stmt);
    }
    else {
        my $stmt = qq{
            SELECT
                lrg_id
            FROM
                gene
            WHERE
                gene_id = $gene_id
            LIMIT 1
        };
        $lrg_id = $db_adaptor->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
    }
    
    return $lrg_id;
}
