#!/usr/bin/env perl

## temp script to create a table of evidence values for mart building



use strict;
use warnings;
use DBI;
use Getopt::Long;

my ($db, $host, $user, $pass);

GetOptions ("db=s"    => \$db,
            "host=s"  => \$host,
            "user=s"  => \$user,
            "pass=s"  => \$pass,
    );

die usage() unless defined $host && defined $user && defined $pass;

my $databases;
if( defined $db){
    push @{$databases}, $db;
}
else{
    ## find all variation databases on the host
    $databases = get_dbs_by_host($host, $user, $pass);
} 

foreach my $db_name (@{$databases}){
    
    my $dbh = DBI->connect( "dbi:mysql:$db_name\:$host\:3306", $user, $pass, undef);

    $dbh->do(qq[update variation set evidence_attribs = NULL where evidence_attribs = '';]);
    $dbh->do(qq[update variation_feature set evidence_attribs = NULL where evidence_attribs = '';]);
    
    $dbh->do(qq[create table MTMP_evidence (
              variation_id int(10) , 
              evidence SET('Multiple_observations','Frequency','HapMap','1000Genomes','Cited','ESP','Phenotype_or_Disease'), 
              primary key( variation_id ) ) ])||die;

    
    my $evidence_id = get_evidence_id($dbh);

    
    my $filename = "$db\.dat";
    open my $out, ">$filename"||die "Failed to open $filename to write evidence statuses : $!\n"; 
    
    my $ev_ext_sth  = $dbh->prepare(qq[ select variation_id, evidence_attribs from variation ]);
    
    my $ev_ins_sth  = $dbh->prepare(qq[ insert into  MTMP_evidence variation_id ,evidence
                                    values (?,?)
                                  ]);
    
    $ev_ext_sth->{mysql_use_result} = 1;
    $ev_ext_sth->execute()||die ;
    

    while ( my $aref = $ev_ext_sth->fetchrow_arrayref() ) {
        
        unless(defined $aref->[1]){
            print $out "$aref->[0]\t\\N\n";
            next;
        }
        
        my @ev_old = split/\,/, $aref->[1];
        
        my @ev_new;
    
        foreach my  $old(@ev_old){ 
            
            die "No id for $old\n" unless defined  $evidence_id->{$old}; 
            push @ev_new, $evidence_id->{$old};             
        }
        
        my $new = join(",", @ev_new);
        
        print $out "$aref->[0]\t$new\n";
        
    }

    close $out;

    $dbh->do( qq[ LOAD DATA LOCAL INFILE "$filename" INTO TABLE  MTMP_evidence]) || die "Error loading $filename data \n";
}


sub get_evidence_id{

    my $dbh = shift;

    my %ids;
    
    my $att_ext_sth  = $dbh->prepare(qq[ select attrib.attrib_id, attrib.value
                                     from attrib, attrib_type
                                     where  attrib_type.code ='evidence'
                                     and attrib.attrib_type_id =attrib_type.attrib_type_id
                                    ]);

    $att_ext_sth->execute()||die;
    my $attdata = $att_ext_sth->fetchall_arrayref();
    foreach my $l(@{$attdata}){

        $ids{$l->[0]} = $l->[1];
    }

    return \%ids;
}


sub usage{

    die "\n\tUsage: create_evidence_MTMP.pl -db [database name] -host [host] -user [write-user name] -pass [write-user password]\n

\t\tOptions: -db [database name]    default: all variation databases on the host\n\n";

}
    
sub get_dbs_by_host{

    my ($host, $user, $pass) = @_;

    my @databases;

    my $dbh = DBI->connect("dbi:mysql:information_schema:$host:3306", $user, $pass, undef);

    my $db_ext_sth = $dbh->prepare(qq[ show databases like '%variation%']);

    $db_ext_sth->execute()||die;
    my $db_list = $db_ext_sth->fetchall_arrayref();
    foreach my $l(@{$db_list}){
        next if $l->[0] =~/master/;
        print "Doing $l->[0] on $host\n";
        push @databases, $l->[0] ;
    }

    return \@databases;
}
