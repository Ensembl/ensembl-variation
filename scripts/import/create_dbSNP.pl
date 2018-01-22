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
 <http://www.ensembl.org/Help/Contact>.



=head1 NAME 

  create_dbSNP.pl

=head1 DESCRIPTION

Download ddl and file listings  for new dbSNP release
and reformat for loading into mysql

=cut

use strict;
use warnings;
use Getopt::Long;
use Net::FTP;
use DBI;


my ($db, $build, $pass, $host, $db_name, $db_load_only );

GetOptions ("db=s"          => \$db,
            "build=s"       => \$build,        
            "pass=s"        => \$pass, 
            "host=s"        => \$host,
            "db_name=s"     => \$db_name,  
	    "db_load_only"  => \$db_load_only 
    );

usage() unless defined  $db && defined $build;

my @required_files ;
my $ddl ;


if($db_name =~ /shared/){
    @required_files = (qw[ Allele GtyAllele ObsVariation PopClass PopClassCode UniGty UniVariation]);
    ## get files from dbSNP ftp site
    $ddl = download_shared_files(\@required_files, $db_load_only) 
}
else{
   @required_files = (qw[ AlleleFreqBySsPop  Batch GtyFreqBySsPop Individual PedigreeIndividual PopLine Population RsMergeArch SNP SNPSubSNPLink SubInd  SubmittedIndividual SubSNP SubSNPPubmed ]);
    
    push @required_files, "b$build\_ContigInfo";
    push @required_files, "b$build\_SNPContigLoc";

    ## get files from dbSNP ftp site
   $ddl = download_organism_files(\@required_files, $db_load_only) 
}

die "Error: no ddl\n" unless defined $ddl;


## convert types in ddl files
my $new_ddl = reformat_ddl_files($ddl);

unless (defined $db_load_only){

    ## check md5sum
    check_files(\@required_files);


    ## convert nulls in data files
    reformat_data_files(\@required_files);
}

## load into mysql
load_data( $db_name, 
           $pass, 
           $host,
           \@required_files, 
           $new_ddl ) if defined $db_name;




sub download_organism_files{
    
    my $required_files = shift;
    my $db_load_only   = shift;

    print localtime(). " downloading files\n";

    my $ftp = Net::FTP->new("ftp.ncbi.nih.gov", Debug => 0 )|| die "Cannot connect to ftp.ncbi.nih.gov: $@";    
    $ftp->login("anonymous",'-anonymous@') || die "Cannot login ", $ftp->message; 
    $ftp->cwd("snp/organisms/$db/database/organism\_data")  || die "Cannot change working directory ", $ftp->message;
    $ftp->binary(); 
    
    unless (defined $db_load_only ){
	## pick up the required data files
	foreach my $file ( @$required_files ){
	    warn "looking for $file\n";
	    $ftp->get( "$file\.bcp.gz" ) || die " Failed to download file :$file\.bcp.gz\n" ;
	    #print $ftp->message() . "\n";
	    $ftp->get( "$file\.bcp.gz.md5" );
	    #print $ftp->message() . "\n";
	}
    }

    ## pick up ddl also
    $ftp->cwd("../organism_schema")  || die "Cannot change working directory ", $ftp->message;
    my $sql_file_list = $ftp->ls();
    die "No ddl\n" unless defined $sql_file_list->[0];
    
    foreach my $file (@{$sql_file_list}){    
        $ftp->get( $file ); 
        warn "picking up $file\n";   
    }
    return $sql_file_list;
}
 

sub download_shared_files{

   my $required_files = shift;
    print localtime(). " downloading files\n";

    my $ftp = Net::FTP->new("ftp.ncbi.nih.gov", Debug => 0 )|| die "Cannot connect to ftp.ncbi.nih.gov: $@";    
    $ftp->login("anonymous",'-anonymous@') || die "Cannot login ", $ftp->message; 
    $ftp->cwd("snp/database/shared_data")  || die "Cannot change working directory ", $ftp->message;
    $ftp->binary(); 
    
    ## pick up the required data files
    foreach my $file ( @$required_files ){
        warn "looking for $file\n";
        $ftp->get( "$file\.bcp.gz" ) || die " Failed to download file :$file\.bcp.gz\n" ;
        #print $ftp->message() . "\n";
        $ftp->get( "$file\.bcp.gz.md5" );
        #print $ftp->message() . "\n";
    }

    ## pick up ddl also
    $ftp->cwd("../shared_schema")  || die "Cannot change working directory ", $ftp->message;
    my $sql_file_list = $ftp->ls();
    die "No ddl\n" unless defined $sql_file_list->[0];
    
    foreach my $file (@{$sql_file_list}){    
        $ftp->get( $file ); 
        warn "picking up $file\n";   
    }
    return $sql_file_list;
}
 

sub check_files{  

    my $required_files = shift;
    print localtime() ." checking downloads\n";

    my $problem = 0;
    foreach my $file ( @$required_files ){
        
        ## extract supplied md5 - weird format
        my $md5_dbSNP; 
        open my $md5_exp, "$file\.bcp.gz.md5" || die "Failed to open $file\.bcp.gz.md5\n";
        while(<$md5_exp>){
            if(/$file/){
                chomp;
                $md5_dbSNP = (split)[1];
            }
        }
        my $md5_output = `md5sum $file\.bcp.gz`;
        my $md5_found = (split/\s+/,$md5_output)[0];
        
        if ($md5_dbSNP eq $md5_found){
            print "md5sum OK for $file\n";
        }
        else{
            $problem = 1;
            warn "Error - $file\tsupplied : $md5_dbSNP, found: $md5_output\n";
        }
    }
    die "Exiting due to errors\n" unless $problem ==0;
}


sub reformat_data_files{

    my $tables = shift;
    print localtime(). " Reformatting data files\n";
    
    foreach my $table (@$tables){
        
        print  "Reformating $table\n";
        open my $in, "gunzip -c $table\.bcp.gz | "||die "Failed to open $table\.bcp.gz for reformatting : $!\n";
        open my $out, ">$table\_new.bcp"||die "Failed to open $table\_new.bcp to write : $!\n";;

        while(<$in>){
            
            my @a = split/\t/;
            
            foreach my $a(@a){
                if($a =~/\w+|\-|\+|\?|\./ && $a =~/\n/){ print $out $a;}
                elsif( $a =~/\w+|\-|\+|\?|\./)         { print $out "$a\t";}
                elsif( $a =~/\n/)   { print $out "\\N\n";}         
                else{                 print $out "\\N\t";}
            }
        }
        close $out;
    }
}

sub reformat_ddl_files{

    my $ddl = shift;
    print localtime() ." Reformatting ddl files\n";

    my @new_ddl;
    foreach my $file(@$ddl){

        my $outname = "mysql_" . $file;
        $outname =~ s/\.gz//;
        push @new_ddl, $outname;
        print  "doing $file => $outname\n";
        open my $out, ">$outname" || die "Failed to open $outname to write: $!\n";
        
        open my $in, "gunzip -c $file |" || die "Failed to open $file to read: $!\n";

        print $out "SET storage_engine=InnoDB;\n\n" if $file=~/table/i;

        
        while(<$in>){
            
            ## remove weird indexes
            next if /getdate|DF__SNP_tax_i__statu__31583BA0/i;  
            next if /i_rsCtgMrna ON b138_SNPContigLocusId/; ## index too long but table not used
            next if /DF__/;
            next if /i_last_updated/;
	    next if /u_Allele_allele/; ## allele strings unique in mssql not unique in mysql
            
            next if /FOREIGN/;
        
            if(/CREATE\s+TABLE/){
                s/CREATE TABLE/CREATE TABLE IF NOT EXISTS/i;
                s/\[|\]//g;
            }
            else{
                
                s/\[|\]//g;
                
                ## end statement with ';' not 'GO'
                s/GO/;/;
                
                ## don't think this is needed
                s/IDENTITY\(1\,1\)//;
                
                ## fix types for table definitions
                s/ bit / TINYINT /;        
                s/uniqueidentifier/CHAR\(36\)/;
                s/image/BLOB/;
                s/smalldatetime/DATETIME/;
                
                s/float/DOUBLE/;
                s/real/FLOAT/;
                s/money|smallmoney/DECIMAL/;
                
                ## set default for empty columns
                s/NULL/DEFAULT NULL/ unless/NOT NULL/;
            
                ## for indexes
                s/NONCLUSTERED|CLUSTERED|ASC//g;
        
            }
            print $out $_;
        }
        print $out "\nSET storage_engine=MYISAM;\n" if $file=~/table/i;
    }
    return \@new_ddl;
}

sub load_data{
    
    my ( $db_name, 
         $pass, 
         $host,
         $data_files, 
         $ddl )   = @_;

    my @keys;
    foreach my $ddl_file ( @{$ddl}){
        if( $ddl_file =~ /table/){
            print localtime() .  " loading table definitions\n";
            `mysql -h$host -uensadmin -p$pass -D $db_name < $ddl_file`
        }
        else{
            push @keys,  $ddl_file;
        }
    }

    my $dir = `pwd`;
    chomp $dir;

    
    ## load data 
    foreach my $table (@{$data_files}){
        warn localtime() . " loading table $table\n";
        `mysql -h$host -uensadmin -p$pass -D $db_name -e "load data local infile '$dir/$table\_new.bcp' into table $table"` ;
    }
    
    ## add indexes
    foreach my $ddl_file ( @{$ddl}){
        unless($ddl_file =~ /table|view/){
            print localtime() . " loading ddl: $ddl_file \n";
            `mysql -h$host -uensadmin -p$pass -D $db_name < $ddl_file`;
        }
    }


    ## print counts for tables as report
    print localtime() . " Getting row counts for new tables \n";

    open my $report, ">$db_name\_report.txt"||die "Failed to open file to create report :$!\n";
    print $report "Entries\tHas PK\tTable name\n\n";

    my $dbh = DBI->connect("dbi:mysql:$db_name:$host\:3306", 'ensadmin', $pass, undef);
       
    foreach my $table (@{$data_files}){
        my $row_number   = count_rows($dbh, $table);
	## some tables never have a PK (but have lots of indexes), some occaisionally have no PK in export
	my $pk_available = "-";
	$pk_available = check_primary_key($dbh, $table) unless $table =~/SubInd|GtyFreqBySsPop|AlleleFreqBySsPop/;
        print $report "$row_number\t$pk_available\t$table\n";       
    }
 
}

sub count_rows{

   my $dbh    = shift;
   my $table_name  = shift;

   my $row_count_ext_sth  = $dbh->prepare(qq[ select count(*) from $table_name]);

   $row_count_ext_sth->execute();
   my $total_rows = $row_count_ext_sth->fetchall_arrayref();

   return $total_rows->[0]->[0]; 

}


sub check_primary_key{

   my $dbh         = shift;
   my $table_name  = shift;

   my $key_ext_sth  = $dbh->prepare(qq[ show indexes from $table_name where Key_name ='PRIMARY']);

   $key_ext_sth->execute();
   my $key = $key_ext_sth->fetchall_arrayref();

   my $found = 0;
   $found = 1 if defined $key->[0]->[0];

   return $found;
}

sub usage{

    die  "\n\tUsage: create_dbSNP.pl -db [zebrafish_7955] -build [138]

\tFor local database creation -pass [admin password]  -host [db server]  -db_name [name for empty mirror database]

\tTo skip data file download/reformatting and only create the local database use -db_load_only\n\n";

}
