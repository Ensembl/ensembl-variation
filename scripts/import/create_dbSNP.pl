#!/usr/bin/env perl

=head1 LICENSE

 Copyright (c) 1999-2012 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.



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

my ($db, $build, $pass, $host, $db_name );

GetOptions ("db=s"          => \$db,
            "build=s"       => \$build,        
	    "pass=s"        => \$pass, 
	    "host=s"        => \$host,
	    "db_name=s"       => \$db_name,  
    );

usage() unless defined  $db && defined $build;



my @required_files = (qw[ AlleleFreqBySsPop  Batch GtyFreqBySsPop Individual PedigreeIndividual PopLine Population RsMergeArch SNP SNPSubSNPLink SubInd  SubmittedIndividual SubSNP SubSNPPubmed ]);

push @required_files, "b$build\_ContigInfo";
push @required_files, "b$build\_SNPContigLoc";


## get file from dbSNP ftp site
my $ddl = download_files(\@required_files);

## check md5sum
check_files(\@required_files);

## convert types in ddl files
my $new_ddl = reformat_ddl_files($ddl);

## convert nulls in data files
reformat_data_files(\@required_files);


## load into mysql
load_data( $db_name, 
	   $pass, 
	   $host,
	   \@required_files, 
	   $new_ddl ) if defined $db_name;


sub download_files{
    
    my $required_files = shift;
    print localtime(). " downloading files\n";

    my $ftp = Net::FTP->new("ftp.ncbi.nih.gov", Debug => 0 )|| die "Cannot connect to ftp.ncbi.nih.gov: $@";    
    $ftp->login("anonymous",'-anonymous@') || die "Cannot login ", $ftp->message; 
    $ftp->cwd("snp/organisms/$db/database/organism\_data")  || die "Cannot change working directory ", $ftp->message;
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
    $ftp->cwd("../organism_schema")  || die "Cannot change working directory ", $ftp->message;
    my $sql_file_list = $ftp->ls();
    
    foreach my $file (@{$sql_file_list}){    
	$ftp->get( $file );    
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
	open my $in, "gunzip -c $table\.bcp.gz | "||die;
	open my $out, ">$table\_new.bcp"||die;

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
	
	my $outname = "mysql" . $file;
	$outname =~ s/\.gz//;
	push @new_ddl, $outname;

	open my $out, ">$outname" || die "Failed to open $outname to write: $!\n";
	
	open my $in, "gunzip -c $file |" || die "Failed to open $file to read: $!\n";
	
	while(<$in>){
	    
	    ## remove weird indexes
	    next if /getdate|DF__SNP_tax_i__statu__31583BA0/i;  
	    next if /i_rsCtgMrna ON b138_SNPContigLocusId/; ## index too long but table not used
	    next if /DF__SubSNPLin__link___56EDABB5/;
	    next if /DF__Submitted__ploid__0DBDF25B/;
	    next if /DF__SubSNPLin__link___359F647E/;
	    next if /DF__Submitted__ploid__7C06F46F/;
	    next if /i_last_updated/;
	    
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
	    `mysql -h$host -uensadmin -p$pass -D $db_name < $ddl_file`;
	}
	else{
	    push @keys,  $ddl_file;
	}
    }

    my $dir = `pwd`;
    chomp $dir;

    
    ## load data 
    foreach my $table (@{$data_files}){
	print localtime() . " loading table $table\n";
	`mysql -h$host -uensadmin -p$pass -D $db_name -e "load data local infile '$dir/$table\_new.bcp' into table $table"` ;
    }
    
    ## add indexes
    foreach my $ddl_file ( @{$ddl}){
	unless($ddl_file =~ /table|view/){
	    print localtime() . " loading ddl: $ddl_file \n";
	    `mysql -h$host -uensadmin -p$pass -D $db_name < $ddl_file`;
	}
    }
}



sub usage{

    die  "\n\tUsage: create_dbSNP.pl -db [zebrafish_7955] -build [138]

\tFor local database creation -pass [admin password]  -host [db server]  -db_name [name for empty mirror database]\n\n";
}
