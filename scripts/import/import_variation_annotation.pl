#! perl -w

use strict;
use Getopt::Long;
use DBI qw(:sql_types);
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $infile;
my $host;
my $dbname;
my $user;
my $pass;
my $port;
our $verbose;
my $skip_synonyms;
my $skip_phenotypes;
my $source;
my $help;

my $UNIPROT_SOURCE_NAME = "Uniprot";
my $UNIPROT_SOURCE_DESCRIPTION = "Variants with protein annotation imported from Uniprot";
my $UNIPROT_SOURCE_URL = "http://www.uniprot.org/";

my $NHGRI_SOURCE_NAME = "NHGRI_GWAS_catalog";
my $NHGRI_SOURCE_DESCRIPTION = "Variants associated with phenotype data from the NHGRI GWAS catalog";
my $NHGRI_SOURCE_URL = "http://www.genome.gov/gwastudies/";

usage() if (!scalar(@ARGV));
 
GetOptions(
    'infile=s' => \$infile,
    'host=s' => \$host,
    'dbname=s' => \$dbname,
    'user=s' => \$user,
    'pass=s' => \$pass,
    'port=s' => \$port,
    'source=s' => \$source,
    'verbose!' => \$verbose,
    'skip_synonyms!' => \$skip_synonyms,
    'skip_phenotypes!' => \$skip_phenotypes,
    'help!' => \$help
);

usage() if ($help);

die ("An input file is required") unless (defined($infile));
die ("Database credentials are required") unless (defined($host) && defined($dbname) && defined($user) && defined($pass));

$port ||= 3306;

my $result;
my $source_name;
my $source_description;
my $source_url;

=head

    The parser subroutines parse the input file into a common data structure.
    Currently what is returned should be a reference to a hash. The hash should
    contain the key 'phenotypes' with the value being a reference to an array of
    phenotype data objects. The phenotype data object is a reference to a hash
    where the keys correspond to the column names in the variation_annotation table.
    In addition, there are the keys 'rsid' which holds the rs-id that the phenotype
    annotates and 'description' and 'name' which correspond to the columns in the
    phenotype table.
    
    In addition, the hash can also contain the key 'synonyms' and the value is a
    reference to a hash where the keys are rs-ids and each respective value is a
    reference to an array of synonyms for the rs-id.

=cut

#ÊParse the input files into a hash
if ($source =~ m/uniprot/i) {
    $result = parse_uniprot($infile);
    $source_name = $UNIPROT_SOURCE_NAME;
    $source_description = $UNIPROT_SOURCE_DESCRIPTION;
    $source_url = $UNIPROT_SOURCE_URL;
}
elsif ($source =~ m/nhgri/i) {
    $result = parse_nhgri($infile);
    $source_name = $NHGRI_SOURCE_NAME;
    $source_description = $NHGRI_SOURCE_DESCRIPTION;
    $source_url = $NHGRI_SOURCE_URL;
}
else {
    die("Source $source is not recognized");
}

my %synonym;
my @phenotypes;
if (exists($result->{'synonyms'})) {
    %synonym = %{$result->{'synonyms'}};
}
if (exists($result->{'phenotypes'})) {
    @phenotypes = @{$result->{'phenotypes'}};
}

# Connect to the variation database
print STDOUT localtime() . "\tConnecting to database $dbname\n" if ($verbose);
my $db_adaptor = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host => $host,
  -user => $user,
  -pass => $pass,
  -port => $port,
  -dbname => $dbname
) or die("Could not get a database adaptor for $dbname on $host:$port");
print STDOUT localtime() . "\tConnected to $dbname on $host:$port\n" if ($verbose);

# Get internal variation ids for the rsIds
my @rsids = map {$_->{'rsid'}} @phenotypes;
my $variation_ids = get_dbIDs(\@rsids,$db_adaptor);

# Get or add a source
my $source_id = get_or_add_source($source_name,$source_description,$source_url,$db_adaptor);
print STDOUT "$source source_id is $source_id\n" if ($verbose);

#ÊAdd the synonyms if required
add_synonyms(\%synonym,$variation_ids,$source_id,$db_adaptor) unless ($skip_synonyms);

# Now, insert phenotypes
add_phenotypes(\@phenotypes,$variation_ids,$source_id,$db_adaptor) unless ($skip_phenotypes);

# Loop over the remaining rsids (the ones that could not be find in the db) and print them out
while (my ($rs_id,$var_id) = each(%{$variation_ids})) {
    next if (defined($var_id->[0]));
    print STDOUT "$rs_id could not be found in $dbname (Synonyms: " . join(", ",@{$synonym{$rs_id}}) . ")\n";
}


sub parse_uniprot {
    my $infile = shift;
    
    my %synonym;
    my @phenotypes;

    #ÊOpen the input file for reading
    open(IN,'<',$infile) or die ("Could not open $infile for reading");
    
    # Read through the file and parse out the desired fields
    while (<IN>) {
        chomp;
        
        # A regexp to catch the meta information in the header. Just echo this to the stdout for logging purposes
        if ($_ =~ m/^\s*(Description|Name|Release)\:/) {
            print STDOUT $_ . "\n";
        }
        
        #ÊMain regexp to extract relevant variation information
        if ($_ =~ m/^(\S+)\s+(?:\w+\s+){2}(VAR\_\d+)\s+\d+\s+\w+ \-> \w+\s+(Disease|Polymorphism|Unclassified)\s+(rs\d+)\s+(.+)$/) {
            
            #ÊGet the data that was caught by the regexp
            my $gene = $1;
            my $uniprot_id = $2;
            my $rs_id = $4;
            my $phenotype = $5;
            
            push(@{$synonym{$rs_id}},$uniprot_id);
            
            # Try to further split the phenotype into a short name, description and possibly MIM id
            if ($phenotype ne '-') {
                my $description;
                my $name;
                my $mim_id;
                
                ($description,$mim_id) = $phenotype =~ m/^([^\[]+)(?:\[(MIM\:.+?)\])?$/;
                ($description,$name) = $description =~ m/^(.+?)\s*(?:\((.+?)\))?\s*$/;
                
                $mim_id &&= join(",MIM:",split(",",$mim_id));
                $mim_id =~ s/\s+//g if (defined($mim_id));
                
                push(
                    @phenotypes,
                    {
                        "rsid" => $rs_id,
                        "associated_gene" => $gene,
                        "description" => $description,
                        "name" => $name,
                        "variation_names" => $rs_id,
                        "study" => $mim_id,
                        "study_type" => 'GWAS'
                    }
                );
            }
        }
    }
    close(IN);
    
    print STDOUT "Parsed " . scalar(keys(%synonym)) . " rs-ids with Uniprot synonyms\n" if ($verbose);
    print STDOUT scalar(@phenotypes) . " phenotypes were found linked to rs-ids\n" if ($verbose);
    
    my %result = ('synonyms' => \%synonym, 'phenotypes' => \@phenotypes);
    return \%result;
}

sub parse_nhgri {
    my $infile = shift;
    
    my @phenotypes;
    
    #ÊOpen the input file for reading
    open(IN,'<',$infile) or die ("Could not open $infile for reading");
    
    # Read through the file and parse out the desired fields
    while (<IN>) {
        chomp;
        
        my (
            $catalog_date,
            $pubmed_id,
            $author,
            $pub_date,
            $journal,
            $url,
            $study,
            $phenotype,
            $initial_sample_size,
            $replication_sample_size,
            $region,
            $gene,
            $rs_risk_allele,
            $rs_id,
            $risk_frequency,
            $pvalue,
            $pval_text,
            $beta,
            $ci,
            $platform,
            $cnv
        ) = split(/\t/,$_);
        
        my %data = (
            'study_type' => 'GWAS',
            'description' => $phenotype,
            'associated_gene' => $gene,
            'associated_variant_risk_allele' => $rs_risk_allele,
            'risk_allele_freq_in_controls' => $risk_frequency,
            'p_value' => $pvalue,
        );
        
        # Parse the rsids
        my @rsids;
        $rs_id ||= "";
        while ($rs_id =~ m/(rs[0-9]+)/g) {
            push(@rsids,$1);
        }
        $data{'variation_names'} = join(',',@rsids);
        $data{'study'} = 'pubmed/' . $pubmed_id if (defined($pubmed_id));
        
        # If we didn't get any rsIds, skip this row (this will also get rid of the header)
        warn("Could not parse any rsIds from string '$rs_id'") if (!scalar(@rsids));
        next if (!scalar(@rsids));
        
        map {
            my %t_data = %{\%data};
            $t_data{'rsid'} = $_;
            push(@phenotypes,\%t_data)
        } @rsids;
    }
    close(IN);
    
    my %result = ('phenotypes' => \@phenotypes);
    
    return \%result;
}

sub parse_ega {
    my $infile;
    
    #ÊOpen the input file for reading
    open(IN,'<',$infile) or die ("Could not open $infile for reading");
    
    # Read through the file and parse out the desired fields
    while (<IN>) {
        chomp;
    }
    close(IN);
    
    my %result;
    
    return\%result;
}
sub get_dbIDs {
    my $rs_ids = shift;
    my $db_adaptor = shift;
    
    my $id_stmt = qq{
        SELECT
            v.variation_id,
            v.name
        FROM
            variation v
        WHERE
            v.name = ?
        LIMIT 1
    };
    my $syn_stmt = qq{
        SELECT
            v.variation_id,
            v.name
        FROM
            variation_synonym vs JOIN
            variation v ON vs.variation_id = v.variation_id
        WHERE
            vs.name = ?
        LIMIT 1
    };
    my $id_sth = $db_adaptor->dbc->prepare($id_stmt);
    my $syn_sth = $db_adaptor->dbc->prepare($syn_stmt);
    
    my %mapping;
    
    foreach my $rs_id (@{$rs_ids}) {
        $id_sth->bind_param(1,$rs_id,SQL_VARCHAR);
        $id_sth->execute();
        my ($var_id,$var_name);
        $id_sth->bind_columns(\$var_id,\$var_name);
        $id_sth->fetch();
        
        # If we couldn't find the rs_id, look in the synonym table
        if (!defined($var_id)) {
            $syn_sth->bind_param(1,$rs_id,SQL_VARCHAR);
            $syn_sth->execute();
            $syn_sth->bind_columns(\$var_id,\$var_name);
            $syn_sth->fetch();
        }
        
        $mapping{$rs_id} = [$var_id,$var_name];
    }
    
    return \%mapping;
}

sub get_or_add_source {
    my $source_name = shift;
    my $source_description = shift;
    my $source_url = shift;
    my $db_adaptor = shift;
    
    my $stmt = qq{
        SELECT
            source_id
        FROM
            source
        WHERE
            name = '$source_name'
        LIMIT 1
    };
    my $sth = $db_adaptor->dbc->prepare($stmt);
    $sth->execute();
    my $source_id;
    $sth->bind_columns(\$source_id);
    $sth->fetch();
    
    if (!defined($source_id)) {
        $stmt = qq{
            INSERT INTO
                source (
                    name,
                    description,
                    url
                )
            VALUES (
                '$source_name',
                '$source_description',
                '$source_url'
            )
        };
        $db_adaptor->dbc->do($stmt);
        $source_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};
        
        print STDOUT "Added source for $source_name (source_id = $source_id)\n" if ($verbose);
    }

    return $source_id;
}

sub add_phenotypes {
    my $phenotypes = shift;
    my $variation_ids = shift;
    my $source_id = shift;
    my $db_adaptor = shift;
    
    # Prepared statements
    my $phen_check_stmt = qq{
        SELECT
            phenotype_id
        FROM
            phenotype
        WHERE
            description = ?
        LIMIT 1
    };
    my $phen_ins_stmt = qq{
        INSERT INTO
            phenotype (
                name,
                description
            )
        VALUES (
            ?,
            ?
        )
    };
    my $va_check_stmt = qq{
        SELECT
            variation_annotation_id
        FROM
            variation_annotation
        WHERE
            variation_id = ? AND
            phenotype_id = ? AND
            source_id = $source_id
        LIMIT 1
    };
    my $va_ins_stmt = qq{
        INSERT INTO
            variation_annotation (
                variation_id,
                phenotype_id,
                source_id,
                study,
                study_type,
                associated_gene,
                associated_variant_risk_allele,
                variation_names,
                risk_allele_freq_in_controls,
                p_value
            )
        VALUES (
            ?,
            ?,
            $source_id,
            ?,
            ?,
            ?,
            ?,
            ?,
            ?,
            ?
        )
    };
    my $phen_check_sth = $db_adaptor->dbc->prepare($phen_check_stmt);
    my $phen_ins_sth = $db_adaptor->dbc->prepare($phen_ins_stmt);
    my $va_check_sth = $db_adaptor->dbc->prepare($va_check_stmt);
    my $va_ins_sth = $db_adaptor->dbc->prepare($va_ins_stmt);

    # First, sort the array according to the phenotype description
    my @sorted = sort {$a->{"description"} cmp $b->{"description"}} @{$phenotypes};
    my $current = "";
    my $phenotype_id;
    my $phenotype_count = 0;
    my $annotation_count = 0;
    
    while (my $phenotype = shift(@sorted)) {
        
        #ÊIf the rs could not be mapped to a variation id, skip it
        next if (!defined($variation_ids->{$phenotype->{"rsid"}}[0]));
        
        # If we have a new phenotype we need to see if it exists in the database and otherwise add it
        if ($phenotype->{"description"} ne $current) {
            undef($phenotype_id);
            $phen_check_sth->bind_param(1,$phenotype->{"description"},SQL_VARCHAR);
            $phen_check_sth->execute();
            $phen_check_sth->bind_columns(\$phenotype_id);
            $phen_check_sth->fetch();
            
            # If no phenotype was found, we need to add it
            if (!defined($phenotype_id)) {
                #$phen_ins_sth->bind_param(1,$phenotype->{"name"},SQL_VARCHAR);
                $phen_ins_sth->bind_param(1,undef,SQL_VARCHAR);
                $phen_ins_sth->bind_param(2,$phenotype->{"description"},SQL_VARCHAR);
                $phen_ins_sth->execute();
                $phenotype_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};
                $phenotype_count++;
            }
            $current = $phenotype->{"description"};
        }
        
        #ÊCheck if this phenotype already exists for this variation and source, in that case we probably want to skip it
        my $va_id;
        $va_check_sth->bind_param(1,$variation_ids->{$phenotype->{"rsid"}}[0],SQL_INTEGER);
        $va_check_sth->bind_param(2,$phenotype_id,SQL_INTEGER);
        $va_check_sth->execute();
        $va_check_sth->bind_columns(\$va_id);
        $va_check_sth->fetch();
        next if (defined($va_id));
        
        # Else, insert this phenotype.
        $va_ins_sth->bind_param(1,$variation_ids->{$phenotype->{"rsid"}}[0],SQL_INTEGER);
        $va_ins_sth->bind_param(2,$phenotype_id,SQL_INTEGER);
        $va_ins_sth->bind_param(3,$phenotype->{"study"},SQL_VARCHAR);
        $va_ins_sth->bind_param(4,$phenotype->{"study_type"},SQL_VARCHAR);
        $va_ins_sth->bind_param(5,$phenotype->{"associated_gene"},SQL_VARCHAR);
        $va_ins_sth->bind_param(6,$phenotype->{"associated_variant_risk_allele"},SQL_VARCHAR);
        $va_ins_sth->bind_param(7,$phenotype->{"variation_names"},SQL_VARCHAR);
        $va_ins_sth->bind_param(8,$phenotype->{"risk_allele_frequency_in_controls"},SQL_VARCHAR);
        $va_ins_sth->bind_param(9,$phenotype->{"p_value"},SQL_VARCHAR);
        $va_ins_sth->execute();
        $annotation_count++;
    }
    print STDOUT "$phenotype_count new phenotypes added\n" if ($verbose);
    print STDOUT "$annotation_count variations were annoteted with phenotypes\n" if ($verbose);
}

sub add_synonyms {
    my $synonyms = shift;
    my $variation_ids = shift;
    my $source_id = shift;
    my $db_adaptor = shift;
    
    #ÊIf we actually didn't get any synonyms, just return
    return if (!defined($synonyms) || !scalar(keys(%{$synonyms})));
    
    #ÊSome prepeared statements needed for inserting the synonyms into database
    my $ins_stmt = qq{
        INSERT IGNORE INTO
          variation_synonym (
            variation_id,
            source_id,
            name
          )
        VALUES (
            ?,
            $source_id,
            ?
        )
    };
    my $ins_sth = $db_adaptor->dbc->prepare($ins_stmt);
    
    my $alt_count = 0;
    my $variation_count = 0;
    
    foreach my $rs_id (keys %{$variation_ids}) {
        
        my $var_id = $variation_ids->{$rs_id}[0];
        
        #ÊIf we have a variation id, we can proceed
        if (defined($var_id)) {
            
            $variation_count++;
            
            $ins_sth->bind_param(1,$var_id,SQL_INTEGER);
            
            #ÊHandle all synonym ids for this rs_id
            while (my $alt_id = shift(@{$synonyms->{$rs_id}})) {
            
                # Add the id as synonym, if it is already present, it will just be ignored
                $ins_sth->bind_param(2,$alt_id,SQL_VARCHAR);
                $ins_sth->execute();
                $alt_count++;
            }
            
        }
    }
    
    print STDOUT "Added $alt_count synonyms for $variation_count rs-ids\n" if ($verbose);
}

sub usage {
	
  print qq{
  Usage: perl import_variation_annotation.pl [OPTION]
  
  Import variation annotation phenotype data into a Variation database
	
  Options:
    
      -verbose		Progress information is printed
      -help		Print this message
      
    Database credentials are specified on the command line
    
      -host		Variation database host name (Required)
      -dbname		Variation database name (Required)
      -user		Variation database user (Required)
      -pass		Variation database password (Required)
      -port		Variation database port (Default: 3306)
      
    An input file must be specified. This file contains the data that will be imported, typically tab-delimited
    and obtained from the UniProt or NHGRI GWAS catalog. If a new source is required, a method for parsing the
    source format into a "standard" data structure can be added.
    
      -infile	        Typically a tab-delimited file (Required)
			
    The source of the data must be specified so that the correct parser and source table entry will be used. Currently
    supported sources are 'uniprot' and 'nhgri'.
    
      -source		String indicating the source of the data (Required)
  } . "\n";
  exit(0);
}
