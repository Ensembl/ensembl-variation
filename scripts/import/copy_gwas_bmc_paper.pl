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
use Getopt::Long;
use DBI qw(:sql_types);
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use FindBin qw( $Bin );
use ImportUtils;

my $host;
my $dbname;
my $user;
my $pass;
my $port;
our $verbose;
my $skip_phenotypes;
my $skip_sets;
my $help;


my $source_name = "Open Access GWAS Database";
my $source_description = "Johnson & O\'Donnell \'An Open Access Database of Genome-wide Association Results\' PMID:19161620";
my $source_url = "http://www.biomedcentral.com/1471-2350/10/6";
my $set = 'ph_johnson_et_al';

my $study_table      = 'study';
my $phenotype_table  = 'phenotype';
my $annotation_table = 'variation_annotation';

usage() if (!scalar(@ARGV));
 
GetOptions(
    'host=s' => \$host,
    'dbname=s' => \$dbname,
    'user=s' => \$user,
    'pass=s' => \$pass,
    'port=s' => \$port,
    'verbose!' => \$verbose,
    'skip_phenotypes!' => \$skip_phenotypes,
		'skip_sets!' => \$skip_sets,
    'help!' => \$help
);

usage() if ($help);

die ("Database credentials are required") unless (defined($host) && defined($dbname) && defined($user) && defined($pass));

$port ||= 3306;

my $registry_file ||= $Bin . "/ensembl.registry";

my $result;

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

# Connect to the new variation database
print STDOUT localtime() . "\tConnecting to database $dbname\n" if ($verbose);
my $db_adaptor = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host => $host,
  -user => $user,
  -pass => $pass,
  -port => $port,
  -dbname => $dbname
) or die("Could not get a database adaptor for $dbname on $host:$port");
print STDOUT localtime() . "\tConnected to $dbname on $host:$port\n" if ($verbose);


# Connect to the current variation database
Bio::EnsEMBL::Registry->load_all( $registry_file );
my $current_db_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor('human','variation');
my $current_dbVar = $current_db_adaptor->dbc->db_handle;

my @rsids;
my $rs_id;

# Parse the input files into a hash
$result = get_gwas_data();

my @phenotypes;
if (exists($result->{'phenotypes'})) {
    @phenotypes = @{$result->{'phenotypes'}};
}

# Get internal variation ids for the rsIds
if (scalar @rsids == 0) {
	@rsids = map {$_->{'rsid'}} @phenotypes;
}
my $variation_ids = get_dbIDs(\@rsids,$db_adaptor);

# Get or add a source
my $source_id = get_or_add_source($source_name,$source_description,$source_url,$db_adaptor);
print STDOUT "The new source_id is $source_id\n" if ($verbose);

# Now, insert phenotypes
add_phenotypes(\@phenotypes,$variation_ids,$source_id,$db_adaptor) unless ($skip_phenotypes);

# Add the variation sets if required
add_set($set,$source_id,$db_adaptor) unless ($skip_sets);

# Loop over the remaining rsids (the ones that could not be find in the db) and print them out
while (my ($rs_id,$var_id) = each(%{$variation_ids})) {
    next if (defined($var_id->[0]));
    print STDOUT "$rs_id could not be found in $dbname\n";
}

sub get_gwas_data {
    
    my @phenotypes;
    my $current_source_sth = $current_dbVar->prepare(qq{SELECT source_id from source where name='Open Access GWAS Database';});
    $current_source_sth->execute();
		my $current_source_id = ($current_source_sth->fetchrow_array())[0];
		
		my $current_data_sth = $current_dbVar->prepare(qq{
			SELECT va.variation_names, va.associated_gene, va.p_value, s.external_reference, p.name, p.description, v.name 
			FROM variation v, study s, variation_annotation va, phenotype p
			WHERE v.variation_id=va.variation_id
						AND s.study_id=va.study_id
						AND va.phenotype_id=p.phenotype_id
						AND s.source_id=$current_source_id
		});
		
		$current_data_sth->execute();

		while (my @res = $current_data_sth->fetchrow_array()) {   
        push(
            @phenotypes,
            {
            		'study_type'      => 'GWAS',
								'name'            => $res[4],
            		'description'     => $res[5],
            		'associated_gene' => $res[1],
            		'p_value'         => $res[2],
								'study'           => $res[3],
								'rsid'            => $res[6],
        				'variation_names' => $res[0]
        		}
				);
    }

    my %result = ('phenotypes' => \@phenotypes);
    
    return \%result;
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
            name = "$source_name"
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
                "$source_description",
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
            $phenotype_table
        WHERE
            description = ?
        LIMIT 1
    };
    my $phen_ins_stmt = qq{
        INSERT INTO
            $phenotype_table (
                name,
                description
            )
        VALUES (
            ?,
            ?
        )
    };
		my $st_ins_stmt = qq{
        INSERT INTO
            $study_table (
                source_id,
                external_reference,
		            study_type
            )
        VALUES (
            $source_id,
            ?,
	          ?
        )
    };
	
    my $va_check_stmt = qq{
        SELECT
            variation_annotation_id
        FROM
            $annotation_table
        WHERE
            variation_id = ? AND
            phenotype_id = ? AND
            study_id = ? AND 
						variation_names = ?
        LIMIT 1
    };
		#AND variation_names = ? 
    my $va_ins_stmt = qq{
        INSERT INTO
            $annotation_table (
                variation_id,
                phenotype_id,
                study_id,
                associated_gene,
                variation_names,
                p_value
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
    my $phen_check_sth = $db_adaptor->dbc->prepare($phen_check_stmt);
    my $phen_ins_sth = $db_adaptor->dbc->prepare($phen_ins_stmt);
    my $st_ins_sth = $db_adaptor->dbc->prepare($st_ins_stmt);
		my $va_check_sth = $db_adaptor->dbc->prepare($va_check_stmt);
    my $va_ins_sth = $db_adaptor->dbc->prepare($va_ins_stmt);

    # First, sort the array according to the phenotype description
    my @sorted = sort {$a->{"description"} cmp $b->{"description"}} @{$phenotypes};
    my $current = "";
    my $phenotype_id;
		my $study_count = 0;
    my $phenotype_count = 0;
    my $annotation_count = 0;
    
    while (my $phenotype = shift(@sorted)) {
		
			# If the rs could not be mapped to a variation id, skip it
			next if (!defined($variation_ids->{$phenotype->{"rsid"}}[0]));
			
			my $sql_study = '= ?';
			my $sql_type = '= ?';
			my $sql_names = '= ?';
      
			# To avoid duplication of study entries
			if (!defined $phenotype->{"study"}) {$sql_study = 'IS NULL'; }
			if (!defined $phenotype->{"study_type"}) {$sql_type = 'IS NULL'; }	
			
			my $st_check_stmt = qq{
        	SELECT
            study_id
        	FROM
            $study_table
        	WHERE
            source_id = $source_id AND
            external_reference $sql_study AND
						study_type $sql_type
        	LIMIT 1
    	};
			my $st_check_sth = $db_adaptor->dbc->prepare($st_check_stmt);
			my $second_param_num = 2;
				
			my $study_id;
			if (defined $phenotype->{"study"}) {
       	$st_check_sth->bind_param(1,$phenotype->{"study"},SQL_VARCHAR) if (defined $phenotype->{"study"});
			}
			else { $second_param_num = 1; }
			 	$st_check_sth->bind_param($second_param_num,$phenotype->{"study_type"},SQL_VARCHAR) if (defined $phenotype->{"study_type"});
       	$st_check_sth->execute();
       	$st_check_sth->bind_columns(\$study_id);
       	$st_check_sth->fetch();
				
      if (!defined($study_id)) {
        $st_ins_sth->bind_param(1,$phenotype->{"study"},SQL_VARCHAR);
        $st_ins_sth->bind_param(2,$phenotype->{"study_type"},SQL_VARCHAR);
        $st_ins_sth->execute();
				$study_count++;
					
        $st_check_sth->bind_param(1,$phenotype->{"study"},SQL_VARCHAR) if (defined $phenotype->{"study"}); ;
				$st_check_sth->bind_param($second_param_num,$phenotype->{"study_type"},SQL_VARCHAR) if (defined $phenotype->{"study_type"});
        $st_check_sth->execute();
        $st_check_sth->bind_columns(\$study_id);
        $st_check_sth->fetch();
			}
        
      # If we have a new phenotype we need to see if it exists in the database and otherwise add it
      if (defined($phenotype->{"description"}) and $phenotype->{"description"} ne $current) {
      	undef($phenotype_id);
        $phen_check_sth->bind_param(1,$phenotype->{"description"},SQL_VARCHAR);
        $phen_check_sth->execute();
        $phen_check_sth->bind_columns(\$phenotype_id);
        $phen_check_sth->fetch();
            
        # If no phenotype was found, we need to add it
        if (!defined($phenotype_id)) {
          $phen_ins_sth->bind_param(1,undef,SQL_VARCHAR);
					$phen_ins_sth->bind_param(2,$phenotype->{"description"},SQL_VARCHAR);
          $phen_ins_sth->execute();
          $phenotype_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};
          $phenotype_count++;
        }
        $current = $phenotype->{"description"};
      }
        
				
      # Else, insert this phenotype.
			$va_ins_sth->bind_param(1,$variation_ids->{$phenotype->{"rsid"}}[0],SQL_INTEGER);
      $va_ins_sth->bind_param(2,$phenotype_id,SQL_INTEGER);
      $va_ins_sth->bind_param(3,$study_id,SQL_INTEGER);
      $va_ins_sth->bind_param(4,$phenotype->{"associated_gene"},SQL_VARCHAR);
      $va_ins_sth->bind_param(5,$phenotype->{"variation_names"},SQL_VARCHAR);
      $va_ins_sth->bind_param(6,$phenotype->{"p_value"},SQL_VARCHAR);
      $va_ins_sth->execute();
      $annotation_count++;
    }
		print STDOUT "$study_count new studies added\n" if ($verbose);
    print STDOUT "$phenotype_count new phenotypes added\n" if ($verbose);
    print STDOUT "$annotation_count variations were annoteted with phenotypes\n" if ($verbose);
}


sub add_set {
  my $set = shift;
  my $source_id = shift;
  my $db_adaptor = shift;
	
	return if (!defined($set));
	
	my $variation_set_id;
	
	# Get variation_set_id
	my $select_set_stmt = qq{
        SELECT v.variation_set_id
        FROM variation_set v, attrib a
        WHERE v.short_name_attrib_id=a.attrib_id 
				  AND a.value = ?
	};
	my $sth1 = $db_adaptor->dbc->prepare($select_set_stmt);
	$sth1->bind_param(1,$set,SQL_VARCHAR);
  $sth1->execute();
  $sth1->bind_columns(\$variation_set_id);
  $sth1->fetch();
	return if (!defined($variation_set_id));
	
	# Insert into variation_set_variation
	my $insert_set_stmt = qq{ 
		INSERT IGNORE INTO variation_set_variation (variation_id,variation_set_id)
			SELECT distinct variation_id, ? 
			FROM variation_annotation WHERE study_id IN
				(SELECT study_id FROM study WHERE source_id=?)
	};
	my $sth2 = $db_adaptor->dbc->prepare($insert_set_stmt);
	$sth2->bind_param(1,$variation_set_id,SQL_INTEGER);
  $sth2->bind_param(2,$source_id,SQL_INTEGER);
  $sth2->execute();
}


sub convert_p_value {

	my $pval = shift;
	
	my $sci_pval = '';
	# If a scientific format is not found, then ...
	if ($pval !~ /^\d+.*e.+$/i) {	
		# If a range format is found (e.g. 10^-2 > p > 10^-3)
		if ($pval =~ /^\d+\^(-\d+)/) {
			if (length("$1")==1) { $1 = "0$1"; } 
			$sci_pval = "1.00e$1"; # e.g 10^-2 > p > 10^-3 => 1.00e-2
		}
		# If a decimal format is found (e.g. 0.0023)
		elsif ($pval =~ /^\d+/){
			$sci_pval = $pval;
		#$sci_pval = sprintf("%.2e",$pval); # e.g. 0.002 => 2,30e-3
		}
		elsif ($pval =~ /^\w+/) {
			$sci_pval = "NULL";
		}
		
	}
	else {
		$pval =~ tr/E/e/;
		if ($pval =~ /^(\d+)(e-?\d+)/) {
			$pval="$1.00$2";	
		}
		if ($pval =~ /^(\d+\.\d{1})(e-?\d+)/) {
			$pval="$1"."0$2";	
		}
		if ($pval =~ /^(\d+\.\d+e-?)(\d{1})$/) {
			$pval = "$1"."0$2";
		}
		$sci_pval = $pval;
	}
	return $sci_pval;
}


sub usage {
	
  print qq{
  Usage: perl copy_gwas_bmc_paper.pl [OPTION]
  
  Import variation annotation phenotype data into a Variation database
	
  Options:
    
      -verbose		       Progress information is printed
      -help		           Print this message
			
			-skip_phenotypes   Skip the study, variation_annotation and phenotype tables insertions.
      
    Database credentials are specified on the command line
    
      -host		  Variation database host name (Required)
      -dbname		Variation database name (Required)
      -user		  Variation database user (Required)
      -pass		  Variation database password (Required)
      -port		  Variation database port (Default: 3306)
  } . "\n";
  exit(0);
}
