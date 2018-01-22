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

## script to create a tables of evidence values, population genotypes and 
## sets for mart building

## evidence tables are always re-created even if present (may change)
## population genotype tables are only created if missing (don't change)
##    - the population genotype table is not created for human


use strict;
use warnings;
use DBI;
use Getopt::Long;

my ($db, $host, $port, $user, $pass, $mode, $tmpdir, $filename);

GetOptions ("db=s"    => \$db,
            "host=s"  => \$host,
            "port=s"  => \$port,
            "user=s"  => \$user,
            "pass=s"  => \$pass,
            "mode:s"  => \$mode,
            "tmpdir:s" => \$tmpdir,
            "tmpfile:s" => \$filename,
    );

die usage() unless defined $host && defined $user && defined $pass && defined $mode;


$tmpdir ||= `pwd`;
chomp $tmpdir;

$port ||= 3306;

my $databases;
if( defined $db){
    push @{$databases}, $db;
}
else{
    ## find all variation databases on the host
    $databases = get_dbs_by_host($host, $port, $user, $pass);
} 

if($mode =~/evidence/){
    create_mtmp_evidence($databases);
}
elsif($mode =~/population_genotype/){
    create_mtmp_population_genotype($databases);
}
elsif($mode eq 'variation_set_variation' || $mode eq 'variation_set_structural_variation'){
    create_mtmp_variation_set($databases, $mode);
}
else{
    warn "\nERROR: Supported mode needed\n\n";
    die usage();
}

sub create_mtmp_evidence{

  my $databases =shift;

  foreach my $db_name (@{$databases}){
    
    my $dbh = DBI->connect( "dbi:mysql:$db_name\:$host\:$port", $user, $pass, undef) || die "Failed to connect to $db_name\n";

    $dbh->do(qq[drop table if exists MTMP_evidence]);    

    ## fetch current list of evidence types & database ids
    my ($type_string, $evidence_id) = get_evidence($dbh);

    $dbh->do( qq[ create table MTMP_evidence (
                  variation_id int(10), 
                  evidence SET$type_string, 
                  primary key( variation_id ) ) 
                ] )||die;

    $filename = "$db_name\.dat" unless defined $filename;
    open my $out, ">$tmpdir/$filename"||die "Failed to open $filename to write evidence statuses : $!\n"; 
    
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

    $dbh->do( qq[ LOAD DATA LOCAL INFILE "$tmpdir/$filename" INTO TABLE  MTMP_evidence]) || die "Error loading $filename data \n";
    unlink "$tmpdir/$filename" || warn "Failed to remove temp file: $tmpdir/$filename :$!\n";
  }
}

sub get_evidence{

  my $dbh = shift;

  my %ids;
  my $type_string = "(";
    
  my $att_ext_sth  = $dbh->prepare(qq[ select attrib.attrib_id, attrib.value
                                     from attrib, attrib_type
                                     where  attrib_type.code ='evidence'
                                     and attrib.attrib_type_id =attrib_type.attrib_type_id
                                    ]);

  $att_ext_sth->execute()||die "Failed to get evidence attribs\n";
  my $attdata = $att_ext_sth->fetchall_arrayref();
  foreach my $l(@{$attdata}){

    ## save id => value mapping to convert between schemaa
    $ids{$l->[0]} = $l->[1];

    ## build set
    $type_string .= "'". $l->[1] . "'," 
  }
    $type_string =~ s/\,$/\)/;

  return ($type_string, \%ids);
}

sub create_mtmp_population_genotype{
  my $databases = shift;

  foreach my $db_name (@{$databases}){

    my $dbh = DBI->connect( "dbi:mysql:$db_name\:$host\:$port", $user, $pass, undef) ||die "Failed to connect to $db_name\n";

    ## no need to re-create if the table is already present for a new import
    my $check_present_sth = $dbh->prepare(qq[show tables like 'MTMP_population_genotype']);
    $check_present_sth->execute()||die "Failed to check for existing tables \n";
    my $dat = $check_present_sth->fetchall_arrayref();
    next if $dat->[0]->[0] eq 'MTMP_population_genotype';

    print "Creating MTMP_population_genotype for $db_name\n";
    $dbh->do(qq[CREATE TABLE `MTMP_population_genotype` (
            `population_genotype_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
            `variation_id` int(10) unsigned NOT NULL,
            `subsnp_id` int(15) unsigned DEFAULT NULL,
            `allele_1` varchar(25000) DEFAULT NULL,
            `allele_2` varchar(25000) DEFAULT NULL,
            `frequency` float DEFAULT NULL,
            `population_id` int(10) unsigned DEFAULT NULL,
            `count` int(10) unsigned DEFAULT NULL,
            PRIMARY KEY (`population_genotype_id`),
            UNIQUE KEY `pop_genotype_idx` (`variation_id`,`subsnp_id`,`frequency`,`population_id`,`allele_1`(5),`allele_2`(5)),
            KEY `variation_idx` (`variation_id`),
            KEY `subsnp_idx` (`subsnp_id`),
            KEY `population_idx` (`population_id`)
            ) ENGINE=MyISAM DEFAULT CHARSET=latin1 ]);

    ## this table is not created for human databases as it is too large to use
    ## and most data is calculated on the fly from VCF
    if ($db_name =~/homo_sapiens/){
      print "Not populating MTMP_population_genotype table for human as not required\n";
      next;
    }

    $dbh->do(qq[INSERT IGNORE INTO MTMP_population_genotype
            SELECT p.population_genotype_id, p.variation_id, p.subsnp_id, 
                   ac1.allele, ac2.allele, p.frequency, p.population_id, p.count
            FROM population_genotype p, genotype_code gc1, genotype_code gc2, allele_code ac1, allele_code ac2
            WHERE p.genotype_code_id = gc1.genotype_code_id  
            AND gc1.haplotype_id = 1 
            AND gc1.allele_code_id = ac1.allele_code_id
            AND p.genotype_code_id = gc2.genotype_code_id 
            AND gc2.haplotype_id = 2 
            AND gc2.allele_code_id = ac2.allele_code_id])|| die "Failed to create MTMP_population_genotype\n";
    }
}
## variation/structural_variation sets can have parents sets
## create MTMP table linking each variation/structural_variation to each set

sub create_mtmp_variation_set{

  my $databases = shift;
  my $table     = shift;

  foreach my $db_name (@{$databases}){

    my $dbh = DBI->connect( "dbi:mysql:$db_name\:$host\:$port", $user, $pass, undef) || die "Failed to connect to $db_name\n";

    my $mtmp_table_name = 'MTMP_'. $table ;
 
    my $object_id = $table;
    $object_id    =~ s/variation_set_//;
    $object_id .= "_id";

    ## production require table to be created newly each time
    $dbh->do(qq[ drop table if exists $mtmp_table_name ]);   
    $dbh->do(qq[ create table $mtmp_table_name like $table ])|| die "Failed to create MTMP_$table" ;

    ## copy direct variant <-> set relationships
    $dbh->do(qq[ insert into $mtmp_table_name select * from $table ])
          || die "Failed to extract sets\n";;
    
    ## add variant <-> parent set relationships
    $dbh->do(qq[ insert ignore into $mtmp_table_name ( $object_id, variation_set_id )
                 select vsv.$object_id, vss.variation_set_super 
                 from $table vsv, 
                      variation_set_structure vss 
                 where vss.variation_set_sub = vsv.variation_set_id 
               ]) || die "Failed to extract first level sub sets\n";

    ## add variant <-> parent of parent set relationships
    $dbh->do(qq[ insert ignore into $mtmp_table_name ($object_id, variation_set_id )
                 select vsv.$object_id, vss2.variation_set_super 
                 from $table vsv, 
                      variation_set_structure vss,
                      variation_set_structure vss2
                 where vss.variation_set_sub = vsv.variation_set_id
                 and vss2.variation_set_sub = vss.variation_set_super
               ])|| die "Failed to extract second level sub sets\n";

  }
}


sub usage{

    die "\n\tUsage: create_MTMP_tables.pl -host [host] 
                                     -user [write-user name] 
                                     -pass [write-user password] 
                                     -mode [evidence|population_genotype|variation_set_variation|variation_set_structural_variation]\n

\t\tOptions: -db [database name]    default: all* variation databases on the host
\t\t         -tmpdir [directory for temp files]
\t\t         -tmpfile [ name for temp files]

\t\t* Note: MTMP_population_genotype is not required for human databases\n\n";

}
    
sub get_dbs_by_host{

    my ($host, $port, $user, $pass) = @_;

    my @databases;

    my $dbh = DBI->connect("dbi:mysql:information_schema:$host:$port", $user, $pass, undef) || die "Failed to look up available databases\n";

    my $db_ext_sth = $dbh->prepare(qq[ show databases like '%variation%']);

    $db_ext_sth->execute()||die "Failed to extract database list\n";
    my $db_list = $db_ext_sth->fetchall_arrayref();
    foreach my $l(@{$db_list}){
        next if $l->[0] =~/master/;
        print "Doing $l->[0] on $host\n";
        push @databases, $l->[0] ;
    }

    return \@databases;
}
