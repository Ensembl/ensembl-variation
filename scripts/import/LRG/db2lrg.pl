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
my $lrg_id;
my $include_external;
my $list_lsdbs;
my @filter_list_name;
my @filter_list_lsdb;
my @filter_list_url;

my $LRG_SCHEMA_VERSION = "1.6";

GetOptions(
  'host=s'		=> \$host,
  'port=i'		=> \$port,
  'dbname=s'		=> \$dbname,
  'user=s'		=> \$user,
  'pass=s'		=> \$pass,
  'xmlfile=s'		=> \$xmlfile,
  'verbose!'		=> \$verbose,
  'hgnc_symbol=s'         => \$hgnc_symbol,
  'lrg_id=s'         => \$lrg_id,
  'include_external!'   => \$include_external,
  'list_lsdbs!'         => \$list_lsdbs,
  'filter_list_lsdb=s'  =>  \@filter_list_lsdb,
  'filter_list_name=s'  =>  \@filter_list_name,
  'filter_list_url=s'  =>  \@filter_list_url
);

die("Database credentials (-host, -port, -dbname, -user) need to be specified!") unless (defined($host) && defined($port) && defined($dbname) && defined($user));
die("An output LRG XML file must be specified") unless (defined($xmlfile) || defined($list_lsdbs));
die("Either the LRG id or an HGNC symbol must be specified") unless (defined($hgnc_symbol) || defined($lrg_id));

#ÊGet a database connection
print STDOUT localtime() . "\tConnecting to database $dbname\n" if ($verbose);
my $db_adaptor = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host => $host,
  -user => $user,
  -pass => $pass,
  -port => $port,
  -dbname => $dbname
) or die("Could not get a database adaptor for $dbname on $host:$port");
print STDOUT localtime() . "\tConnected to $dbname on $host:$port\n" if ($verbose);

# Get the gene_id
my $stmt;
my $sth;
my $gene_id;
if (defined($lrg_id)) {
    $stmt = qq{
        SELECT
            gene_id,
            symbol,
            lrg_id
        FROM
            gene
        WHERE
            lrg_id = '$lrg_id'
    };
    $stmt .= qq{ AND symbol = '$hgnc_symbol'} if (defined($hgnc_symbol));
    $stmt .= qq{ LIMIT 1};
    ($gene_id,$hgnc_symbol,$lrg_id) = @{$db_adaptor->dbc->db_handle->selectall_arrayref($stmt)->[0]} or die ("Could not find gene corresponding to LRG id $lrg_id" . (defined($hgnc_symbol) ? " and HGNC symbol $hgnc_symbol" : ""));
}
else {
    $stmt = qq{
        SELECT
            gene_id,
            symbol,
            lrg_id
        FROM
            gene
        WHERE
            symbol = '$hgnc_symbol'
        LIMIT 1
    };
    ($gene_id,$hgnc_symbol,$lrg_id) = @{$db_adaptor->dbc->db_handle->selectall_arrayref($stmt)->[0]} or die ("Could not find gene corresponding to HGNC symbol $hgnc_symbol");
}

# If a list of LSDBs were asked for, print that
if (defined($list_lsdbs)) {
    $stmt = qq{
        SELECT
            lg.lsdb_id
        FROM
            lsdb_gene lg JOIN
            lsdb l USING (lsdb_id) LEFT JOIN
            lsdb_contact lc USING (lsdb_id) LEFT JOIN
            contact c USING (contact_id)
        WHERE
            lg.gene_id = $gene_id AND
            (
                l.name NOT LIKE 'NCBI RefSeqGene' AND
                l.name NOT LIKE 'LRG' AND
                l.name NOT LIKE 'Ensembl'
            )
    };
    my $condition = join("' AND c.name NOT LIKE '",@filter_list_name);
    $stmt .= qq{ AND (c.name IS NULL OR (c.name NOT LIKE '$condition')) } if (length($condition) > 0);
    
    $condition = join("' AND l.name NOT LIKE '",@filter_list_lsdb);
    $stmt .= qq{ AND (l.name NOT LIKE '$condition') } if (length($condition) > 0);
    
    $condition = join("' AND l.url NOT LIKE '",@filter_list_url);
    $stmt .= qq{ AND (l.url IS NULL OR (l.url NOT LIKE '$condition')) } if (length($condition) > 0);
    
    $stmt .= qq{
        GROUP BY
            lg.lsdb_id
        ORDER BY
            l.name ASC,
            l.lsdb_id ASC
    };
    
    $sth = $db_adaptor->dbc->prepare($stmt);
    $sth->execute();
    my $lsdb_id;
    $sth->bind_columns(\$lsdb_id);
    
    while ($sth->fetch()) {
        my $lsdb = get_source($lsdb_id,$db_adaptor);
        next unless (defined($lsdb));
        
        foreach my $field (('name','url')) {
            my $node = $lsdb->findNode($field);
            print STDOUT ucfirst($field) . ":\t" . $node->content() . "\n" if (defined($node));
        }
        my $contacts = $lsdb->findNodeArray('contact');
        while (my $contact = shift(@{$contacts})) {
            print STDOUT "\n";
            foreach my $field (('name','address','email','url','phone','fax')) {
                my $node = $contact->findNode($field);
                print STDOUT "\t" . ucfirst($field) . ":\t" . $node->content() . "\n" if (defined($node));
            }
        }
        print STDOUT "\n";
    }
    exit;
}

#ÊStatement to get the lrg_data
$stmt = qq{
    SELECT
        organism,
        taxon_id,
        moltype,
        creation_date
    FROM
        lrg_data
    WHERE
        gene_id = '$gene_id'
    LIMIT 1
};
$sth = $db_adaptor->dbc->prepare($stmt);
$sth->execute();

my ($organism,$taxon_id,$moltype,$creation_date);
$sth->bind_columns(\$organism,\$taxon_id,\$moltype,\$creation_date);
$sth->fetch();

# Create LRG root elements and fixed_annotation element
my $lrg_root = LRG::LRG::new($xmlfile);
my $lrg = $lrg_root->addNode('lrg',{'schema_version' => $LRG_SCHEMA_VERSION});
my $fixed = $lrg->addNode('fixed_annotation');

# Set the data
$fixed->addNode('id')->content($lrg_id);
$fixed->addNode('organism',{'taxon' => $taxon_id})->content($organism);

$stmt = qq{
    SELECT DISTINCT
        lr.lsdb_id
    FROM
        lrg_request lr JOIN
        lsdb l USING (lsdb_id)
    WHERE
        lr.gene_id = $gene_id
    ORDER BY
        (l.name = 'NCBI RefSeqGene') ASC,
        l.name ASC
};
$sth = $db_adaptor->dbc->prepare($stmt);
$sth->execute();
my $lsdb_id;
$sth->bind_columns(\$lsdb_id);
while ($sth->fetch()) {
    $fixed->addExisting(get_source($lsdb_id,$db_adaptor));
}

$fixed->addNode('mol_type')->content($moltype);
$fixed->addNode('creation_date')->content($creation_date);
my $sequence = get_sequence($gene_id,'genomic',$db_adaptor);
$fixed->addNode('sequence')->content($sequence);

my $c_stmt = qq{
    SELECT
        c.codon,
        c.type
    FROM
        lrg_cds_codon c
    WHERE
        c.cds_id = ?
    ORDER BY
        c.type DESC,
        c.codon ASC
};
my $e_stmt = qq{
    SELECT
        e.lrg_start,
        e.lrg_end,
        e.cdna_start,
        e.cdna_end,
        e.peptide_start,
        e.peptide_end,
        i.phase
    FROM
        lrg_exon e LEFT JOIN
        lrg_intron i ON i.exon_5 = e.exon_id
    WHERE
        e.transcript_id = ?
    ORDER BY
        e.lrg_start ASC
};
$stmt = qq{
    SELECT
        lt.transcript_id,
        lt.transcript_name,
        cdna.cdna_id,
        cdna.lrg_start,
        cdna.lrg_end,
        cds.cds_id,
        cds.lrg_start,
        cds.lrg_end,
        cds.codon_start
    FROM
        lrg_transcript lt JOIN
        lrg_cdna cdna USING (transcript_id) JOIN
        lrg_cds cds USING (transcript_id)
    WHERE
        lt.gene_id = '$gene_id'
    ORDER BY
        lt.transcript_name ASC
};
my $c_sth = $db_adaptor->dbc->prepare($c_stmt);
my $e_sth = $db_adaptor->dbc->prepare($e_stmt);
$sth = $db_adaptor->dbc->prepare($stmt);
$sth->execute();
my ($t_id,$t_name,$cdna_id,$cdna_start,$cdna_end,$cds_id,$cds_lrg_start,$cds_lrg_end,$codon_start);
$sth->bind_columns(\$t_id,\$t_name,\$cdna_id,\$cdna_start,\$cdna_end,\$cds_id,\$cds_lrg_start,\$cds_lrg_end,\$codon_start);
while ($sth->fetch()) {
    my $transcript = $fixed->addNode('transcript',{'name' => $t_name, 'start' => $cdna_start, 'end' => $cdna_end});
    my $cdna_seq = get_sequence($cdna_id,'cdna',$db_adaptor);
    $transcript->addNode('cdna/sequence')->content($cdna_seq);
    my $cds = $transcript->addNode('coding_region',{'start' => $cds_lrg_start, 'end' => $cds_lrg_end});
    $cds->addData({'codon_start' => $codon_start}) if (defined($codon_start));
    
    # Check for non-standard codons
    $c_sth->bind_param(1,$cds_id,SQL_INTEGER);
    $c_sth->execute();
    my ($c_codon,$c_type);
    $c_sth->bind_columns(\$c_codon,\$c_type);
    while ($c_sth->fetch()) {
        $cds->addEmptyNode($c_type,{'codon' => $c_codon});
    }
    
    # Add translation element
    my $cds_seq = get_sequence($cds_id,'peptide',$db_adaptor);
    $cds->addNode('translation/sequence')->content($cds_seq);
    
    $e_sth->bind_param(1,$t_id,SQL_INTEGER);
    $e_sth->execute();
    my ($e_lrg_start,$e_lrg_end,$e_cdna_start,$e_cdna_end,$e_peptide_start,$e_peptide_end,$phase);
    $e_sth->bind_columns(\$e_lrg_start,\$e_lrg_end,\$e_cdna_start,\$e_cdna_end,\$e_peptide_start,\$e_peptide_end,\$phase);
    while ($e_sth->fetch()) {
        my $exon = $transcript->addNode('exon');
        $exon->addEmptyNode('lrg_coords',{'start' => $e_lrg_start,'end' => $e_lrg_end});
        $exon->addEmptyNode('cdna_coords',{'start' => $e_cdna_start,'end' => $e_cdna_end});
        $exon->addEmptyNode('peptide_coords',{'start' => $e_peptide_start,'end' => $e_peptide_end}) if (defined($e_peptide_start) && defined($e_peptide_end));
        $transcript->addEmptyNode('intron',{'phase' => $phase}) if (defined($phase));
    }
}

# Create the updatable section
my $updatable = $lrg->addNode('updatable_annotation');

# Get the annotation sets
my $asm_stmt = qq{
    SELECT
        lm.mapping_id
    FROM
        lrg_annotation_set_mapping lasm JOIN
        lrg_mapping lm USING (mapping_id)
    WHERE
        lasm.annotation_set_id = ?
    ORDER BY
        lm.most_recent DESC,
        lm.assembly ASC
};
$stmt = qq{
    SELECT
        annotation_set_id,
        source,
        comment,
        modification_date,
        lrg_gene_name,
        xml
    FROM
        lrg_annotation_set
    WHERE
        gene_id = '$gene_id'
    ORDER BY
        annotation_set_id ASC
};
my $asm_sth = $db_adaptor->dbc->prepare($asm_stmt);
$sth = $db_adaptor->dbc->prepare($stmt);
$sth->execute();
my ($annotation_set_id,$comment,$modification_date,$lrg_gene_name,$xml);
$sth->bind_columns(\$annotation_set_id,\$lsdb_id,\$comment,\$modification_date,\$lrg_gene_name,\$xml);
while ($sth->fetch()) {
    my $annotation_set = $updatable->addNode('annotation_set');
    # Add source information
    $annotation_set->addExisting(get_source($lsdb_id,$db_adaptor));
    # Add comment
    $annotation_set->addNode('comment')->content($comment) if (defined($comment));
    # Add modification date
    $annotation_set->addNode('modification_date')->content($modification_date);
    
    # Add mapping if necessary
    $asm_sth->bind_param(1,$annotation_set_id,SQL_INTEGER);
    $asm_sth->execute();
    my $mapping_id;
    $asm_sth->bind_columns(\$mapping_id);
    while ($asm_sth->fetch()) {
        $annotation_set->addExisting(get_mapping($mapping_id,$gene_id,$db_adaptor));
    }
    
    # Add lrg_gene_name
    $annotation_set->addNode('lrg_gene_name',{'source' => 'HGNC'})->content($lrg_gene_name) if (defined($lrg_gene_name));
    # Add the remaining XML
    my $lrg = LRG::LRG::newFromString($xml);
    while (my $node = shift(@{$lrg->{'nodes'}})) {
        $annotation_set->addExisting($node);
    }
}

# If we should add external LSDBs, add these as separate annotation sets
if (defined($include_external)) {
    $stmt = qq{
        SELECT DISTINCT
            lg.lsdb_id
        FROM
            lsdb_gene lg
        WHERE
            lg.gene_id = $gene_id AND
            NOT EXISTS (
                SELECT
                    *
                FROM
                    lrg_request lr
                WHERE
                    lr.gene_id = '$gene_id' AND
                    lr.lsdb_id = lg.lsdb_id
            ) AND
            NOT EXISTS (
                SELECT
                    *
                FROM
                    lrg_annotation_set las
                WHERE
                    las.gene_id = '$gene_id' AND
                    las.source = lg.lsdb_id
            )
        ORDER BY
            lg.lsdb_id ASC
    };
    $sth = $db_adaptor->dbc->prepare($stmt);
    $sth->execute();
    $sth->bind_columns(\$lsdb_id);
    while ($sth->fetch()) {
        my $annotation_set = $updatable->addNode('annotation_set');
        # Add source information
        $annotation_set->addExisting(get_source($lsdb_id,$db_adaptor,1));
        # Add modification date
        $annotation_set->addNode('modification_date')->content(LRG::LRG::date());
    }
}

# Dump XML to output_file
$lrg_root->printAll();


sub get_mapping {
    my $mapping_id = shift;
    my $gene_id = shift;
    my $db_adaptor = shift;
    
    # Get the mapping data
    my $m_stmt = qq{
        SELECT        
            assembly,
            chr_name,
            chr_id,
            chr_start,
            chr_end,
            most_recent
        FROM
            lrg_mapping
        WHERE
            mapping_id = $mapping_id AND
            gene_id = '$gene_id'
        LIMIT 1
    };
    my $ms_stmt = qq{
        SELECT
            mapping_span_id,
            lrg_start,
            lrg_end,
            chr_start,
            chr_end,
            strand
        FROM
            lrg_mapping_span
        WHERE
            mapping_id = $mapping_id
        ORDER BY
            lrg_start ASC
    };
    my $md_stmt = qq{
        SELECT
            type,
            chr_start,
            chr_end,
            lrg_start,
            lrg_end,
            lrg_sequence,
            chr_sequence
        FROM
            lrg_mapping_diff
        WHERE
            mapping_span_id = ?
        ORDER BY
            lrg_start ASC
    };
    my $ms_sth = $db_adaptor->dbc->prepare($ms_stmt);
    my $md_sth = $db_adaptor->dbc->prepare($md_stmt);
    
    my ($assembly,$chr_name,$chr_id,$chr_start,$chr_end,$most_recent) = @{$db_adaptor->dbc->db_handle->selectall_arrayref($m_stmt)->[0]};
    my $mapping = LRG::Node::new('mapping');
    $mapping->addData({'assembly' => $assembly,'chr_name' => $chr_name,'chr_start' => $chr_start, 'chr_end' => $chr_end});
    $mapping->addData({'chr_id' => $chr_id}) if (defined($chr_id));
    $mapping->addData({'most_recent' => $most_recent}) if ($most_recent);
    
    $ms_sth->execute();
    my ($mapping_span_id,$lrg_start,$lrg_end,$strand);
    $ms_sth->bind_columns(\$mapping_span_id,\$lrg_start,\$lrg_end,\$chr_start,\$chr_end,\$strand);
    while ($ms_sth->fetch()) {
        my $span = $mapping->addNode('mapping_span',{'lrg_start' => $lrg_start, 'lrg_end' => $lrg_end, 'start' => $chr_start, 'end' => $chr_end, 'strand' => $strand});
        
        $md_sth->bind_param(1,$mapping_span_id,SQL_INTEGER);
        $md_sth->execute();
        my ($md_type,$md_chr_start,$md_chr_end,$md_lrg_start,$md_lrg_end,$md_lrg_sequence,$md_chr_sequence);
        $md_sth->bind_columns(\$md_type,\$md_chr_start,\$md_chr_end,\$md_lrg_start,\$md_lrg_end,\$md_lrg_sequence,\$md_chr_sequence);
        while ($md_sth->fetch()) {
            my $diff = $span->addEmptyNode('diff',{'type' => $md_type, 'lrg_start' => $md_lrg_start, 'lrg_end' => $md_lrg_end, 'start' => $md_chr_start, 'end' => $md_chr_end});
            $diff->addData({'lrg_sequence' => $md_lrg_sequence}) if (defined($md_lrg_sequence));
            $diff->addData({'genomic_sequence' => $md_chr_sequence}) if (defined($md_chr_sequence));
        }
    }
    
    return $mapping;
}

sub get_source {
    my $lsdb_id = shift;
    my $db_adaptor = shift;
    my $skip_contact = shift;
    
    #ÊGet the lsdb data
    my $stmt = qq{
        SELECT
            l.name,
            l.url
        FROM
            lsdb l
        WHERE
            l.lsdb_id = $lsdb_id
    };
    my ($lsdb_name,$lsdb_url) = @{$db_adaptor->dbc->db_handle->selectall_arrayref($stmt)->[0]};
    
    # Create the source node
    my $source = LRG::Node::new('source');
    $source->addNode('name')->content($lsdb_name);
    $source->addNode('url')->content($lsdb_url) if (defined($lsdb_url) && length($lsdb_url) > 0);
    
    #ÊIf we skip contact information, return here
    return $source if ($skip_contact);
    
    # Get the source data for the requesters
    $stmt = qq{
        SELECT DISTINCT
            c.name,
            c.email,
            c.url,
            c.address
        FROM
            lsdb_contact lc JOIN
            contact c USING (contact_id)
        WHERE
            lc.lsdb_id = '$lsdb_id'
        ORDER BY
            c.name ASC
    };
    my $sth = $db_adaptor->dbc->prepare($stmt);
    
    $sth->execute();
    my ($name,$email,$url,$address);
    $sth->bind_columns(\$name,\$email,\$url,\$address);
    while ($sth->fetch()) {
        my $contact = $source->addNode('contact');
        $contact->addNode('name')->content($name) if (defined($name));
        $contact->addNode('url')->content($url) if (defined($url));
        $contact->addNode('address')->content($address) if (defined($address));
        $contact->addNode('email')->content($email) if (defined($email));
    }
    
    return $source;
}

sub get_sequence {
    my $id = shift;
    my $type = shift;
    my $db_adaptor = shift;
    
    my $stmt;
    
    if ($type =~ m/^genomic$/i) {
        $stmt = qq{
            SELECT
                ls.sequence
            FROM
                lrg_genomic_sequence lgs JOIN
                lrg_sequence ls USING (sequence_id)
            WHERE
                lgs.gene_id = '$id'
            ORDER BY
                ls.sequence_id
        };
    }
    elsif ($type =~ m/^cdna$/i) {
        $stmt = qq{
            SELECT
                ls.sequence
            FROM
                lrg_cdna_sequence lcs JOIN
                lrg_sequence ls USING (sequence_id)
            WHERE
                lcs.cdna_id = '$id'
            ORDER BY
                ls.sequence_id
        };
    }
    elsif ($type =~ m/^peptide$/i) {
        $stmt = qq{
            SELECT
                ls.sequence
            FROM
                lrg_peptide_sequence lps JOIN
                lrg_sequence ls USING (sequence_id)
            WHERE
                lps.cds_id = '$id'
            ORDER BY
                ls.sequence_id
        };
    }
    else {
        warn ("Unknown sequence type '$type' specified");
        return undef;
    }
    my $sth = $db_adaptor->dbc->prepare($stmt);
    $sth->execute();
    my ($seq,$substr);
    $sth->bind_columns(\$substr);
    while ($sth->fetch()) {
        $seq .= $substr;
    }
    
    return $seq;
}
