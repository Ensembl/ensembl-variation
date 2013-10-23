# Copyright 2013 Ensembl
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

#/usr/bin/env perl


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

use strict;
use warnings;

use Getopt::Long;
use ImportUtils qw(debug load);
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);
#use DBH;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use FindBin qw( $Bin );
use DBI qw(:sql_types);
use Data::Dumper;

my ($TMP_DIR, $TMP_FILE);
my $input_file; #file containing the raw data

my ($vhost, $vport, $vdbname, $vuser, $vpass,
    $chost, $cport, $cdbname, $cuser, $cpass);


GetOptions('tmpdir=s'  => \$ImportUtils::TMP_DIR,
	   'tmpfile=s' => \$ImportUtils::TMP_FILE,
	   'inputfile=s' => \$input_file
	   );

warn("Make sure you have a updated ensembl.registry file!\n");

my $registry_file ||= $Bin . "/ensembl.registry";
my $species = 'mouse';
Bio::EnsEMBL::Registry->load_all( $registry_file );

my $dbCore = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $dbVar = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');


#added default options
$chost = $dbCore->dbc->host; 
$cuser    ||= 'ensro';
$cport = $dbCore->dbc->port ;
$cdbname = $dbCore->dbc->dbname;

$vhost = $dbVar->dbc->host;
$vport = $dbVar->dbc->port;
$vuser    ||= 'ensadmin';
$vdbname = $dbVar->dbc->dbname;
$vpass = $dbVar->dbc->password;

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

my $region_names = {}; #reference to a hash containing the seq_region_name => seq_region_id
my $variation_id = 0; #id for the SNP
my $variation_feature_id = 0;
my $allele_id = 0;

my $source_id = create_source($dbVar,"Sanger");
$region_names = get_region_names($dbCore);


#copy the file to the /tmp in the machine
system("lsrcp ecs4a:$TMP_DIR/$input_file /tmp/$input_file");
print STDERR "Parsing file: $input_file at \t", scalar(localtime),"\n";
parse_file('/tmp/' . $input_file, $dbCore,$dbVar,$source_id, $region_names);
print STDERR "Finished parsing file $input_file \t", scalar(localtime),"\n";
system("rm /tmp/$input_file");
import_files($dbVar,"variation", "variation_feature", "allele", "flanking_sequence", "tmp_individual_genotype_single_bp");


sub parse_file{
    my $filename = shift;
    my $dbCore = shift;
    my $dbVar = shift;
    my $source_id = shift;
    my $region_names = shift;

    my $strain_names = {}; #reference to the hash containing the strain_name=>sample_id

    open IN, "<$filename" || die "Could not open file $filename: $!\n";   
#SNP: 0 C C57BL/6J 1:3427833 NT_039169.6_427833 3 2 1 1
#     1 T 129X1/SvJ 1:3427833 NT_039169.6_427833 3 2 1 1
#     0 C A/J 1:3427833 NT_039169.6_427833 3 2 1 0
#     1 T DBA/2J 1:3427833 NT_039169.6_427833 3 2 1 1
    my $previous_variation_name = '';
    my $previous_seq_region_id;
    my $previous_seq_region_start;
    my $previous_seq_region_end;
    my $ref_allele;
    my %alleles;
    while (<IN>){
	next if (/^\n/); #jump empty lines	
	chomp;
	/[SNP:]?\s\d+\s(\w)\s(\S+)\s(\S+)\s(\S+)/;
	my $allele = $1;
	my $strain_name = $2;
# 	if ($strain_name eq 'MSM'){
# 	    $strain_name = 'MSM/Ms';
# 	}
# 	elsif ($strain_name eq 'C3H'){
# 	    $strain_name = 'C3H/HeJ';
# 	}
# 	elsif ($strain_name eq 'NOD'){
# 	    $strain_name = 'NOD/DIL';
# 	}
	my $position = $3;
	my $name = $4;
	#create the different information: variation, variation_feature, flanking_sequence, allele, sample, population, individual,
	# tmp_individual_genotype_single_bp, individual_population
	#first, try to create a new entry in the individual and population table if not already present the strain
	if ($previous_variation_name eq '' || $previous_variation_name ne $name){
	    #this is a new variation, the first line contain information about the reference allele
	    if ($previous_variation_name ne $name && $previous_variation_name ne ''){
		#update variation_feature table
		if (keys %alleles != 2){
		    warn "Variation $previous_variation_name has ", scalar(keys(%alleles)),"\n";		    
		}
		delete $alleles{$ref_allele};
		my $allele_string = join("/",$ref_allele,keys %alleles);
		write_file("variation_feature",++$variation_feature_id,$previous_seq_region_id,$previous_seq_region_start,$previous_seq_region_end,1,$variation_id,$allele_string,$previous_variation_name,1,1,$source_id,'\N','INTERGENIC');
		#flush the hash with the alleles
		%alleles = ();
	    }
	    #update variation and flanking_sequence tables
	    #variation information
	    write_file("variation",++$variation_id,$source_id,$name,'\N','\N');
	    my ($seq_region_id, $seq_region_start, $seq_region_end) = get_position($position,$region_names);#get the seq_region_id and the position
	    write_file("flanking_sequence",$variation_id,'\N','\N',$seq_region_start - 100,$seq_region_start -1,$seq_region_end +1,$seq_region_end+100,$seq_region_id,1);
	    
	    $previous_variation_name = $name;
	    $previous_seq_region_id = $seq_region_id;
	    $previous_seq_region_start = $seq_region_start;
	    $previous_seq_region_end = $seq_region_end;
	    $alleles{$allele}++;
	    $ref_allele = $allele; #store the ref allele, has to be the first in the string
	}
	else{
	    if (!defined $strain_names->{$strain_name}){
		$strain_names->{$strain_name}->{population} = insert_population($dbVar,$strain_name);
		$strain_names->{$strain_name}->{individual} = insert_individual($dbVar,$strain_name);
		insert_individual_population($dbVar,$strain_names->{$strain_name}->{individual},$strain_names->{$strain_name}->{population});
	    }
	    #and insert in the allele and tmp_genotype tables
	    write_file("allele",++$allele_id,$variation_id,$allele,1,$strain_names->{$strain_name}->{population});
	    write_file("tmp_individual_genotype_single_bp",$variation_id,$allele,$allele,$strain_names->{$strain_name}->{individual});
	    $alleles{$allele}++; #add the allele
	}
    }
    
    close IN || die "Could not close file $filename: $!\n";
}


#for a given position region_name:position, returns in ensembl format region_id,start,end
sub get_position{
    my $position = shift;
    $region_names = shift;
    my ($region_name,$region_start) = split /:/,$position;
    if (!defined $region_names->{$region_name}){
	warn "Unabla to find slice for $region_name\n";
	die;
    }
    return ($region_names->{$region_name},$region_start,$region_start);
}

#inserts a new source in the database
sub create_source{
    my $dbVar = shift;
    my $name = shift;

    my $sth = $dbVar->dbc->prepare(qq{INSERT INTO source (name) values (?) });
    $sth->bind_param(1,$name,SQL_VARCHAR);
    $sth->execute();
    return $dbVar->dbc()->db_handle->{'mysql_insertid'}; #get the id inserted    
}

#inserts a new strain in the database: tables Sample and Population
sub insert_population{
    my $dbVar = shift;
    my $strain_name = shift;

    my $sample_id;
    my $sth = $dbVar->dbc->prepare(qq{INSERT INTO sample (name,description) values (?,?) });
    $sth->bind_param(1,$strain_name,SQL_VARCHAR);
    $sth->bind_param(2,"Population $strain_name",SQL_VARCHAR);
    $sth->execute();
    $sample_id = $dbVar->dbc()->db_handle->{'mysql_insertid'}; #get the id inserted      
    #and insert in the population table
    $sth = $dbVar->dbc->prepare(qq{INSERT INTO population (sample_id, is_strain) VALUES (?,1)});
    $sth->bind_param(1,$sample_id,SQL_INTEGER);
    $sth->execute();
    $sth->finish();
    
    return $sample_id;
}

#inserts a new strain in the database: tables Sample and Individual
sub insert_individual{
    my $dbVar = shift;
    my $strain_name = shift;

    my $sample_id;
    my $sth = $dbVar->dbc->prepare(qq{INSERT INTO sample (name,description) values (?,?) });
    $sth->bind_param(1,$strain_name,SQL_VARCHAR);
    $sth->bind_param(2,"Individual $strain_name",SQL_VARCHAR);
    $sth->execute();
    $sample_id = $dbVar->dbc()->db_handle->{'mysql_insertid'}; #get the id inserted      
    #and insert in the individual table
    $sth = $dbVar->dbc->prepare(qq{INSERT INTO individual (sample_id) VALUES (?)});
    $sth->bind_param(1,$sample_id,SQL_INTEGER);
    $sth->execute();
    $sth->finish();
      
    return $sample_id;
}

sub insert_individual_population{
    my $dbVar = shift;
    my $individual_sample_id = shift;
    my $population_sample_id = shift;
    
    my $sth = $dbVar->dbc->prepare(qq{INSERT INTO individual_population (individual_sample_id, population_sample_id) values (?,?)});
    $sth->bind_param(1,$individual_sample_id,SQL_INTEGER);
    $sth->bind_param(2,$population_sample_id,SQL_INTEGER);
    $sth->execute();
    $sth->finish();
}

sub write_file{
    my $filename = shift;
    my @values = @_;

    open FH, ">>/tmp/$filename.txt" || die "Could not open file with $filename information: $!\n";
    print FH join("\t", @values), "\n";
    close FH || die "Could not close file with $filename information: $!\n";
    
}

#for a given list of files, load them in the database
sub import_files{
    my $dbVar = shift;
    my @files_to_import = @_;
    
    foreach my $file (@files_to_import){
	debug("Importing $file");
	system("lsrcp /tmp/$file.txt ecs4a:$TMP_DIR/$TMP_FILE");
	system("rm /tmp/$file.txt");
	load($dbVar->dbc->db_handle,$file);
	
    }
}

sub get_region_names{
    my $dbCore = shift;
    
    my %region_names;
    my $slice_adaptor = $dbCore->get_SliceAdaptor();
    map {$region_names{$_->seq_region_name} = $_->get_seq_region_id} @{$slice_adaptor->fetch_all('toplevel')};
    return \%region_names;    
}

