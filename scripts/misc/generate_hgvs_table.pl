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


## script to create a table of HGVS stings and variation_ids for search index building


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

use strict;
use warnings;

use DBI;
use Getopt::Long;

my ($host, $user, $pass, $db);

GetOptions ("host=s"        => \$host,	    
            "pass=s"        => \$pass,
            "db=s"          => \$db
	    );

usage() unless defined $pass;


if(defined $host && defined $db){

    create_table( $db, $host, $pass);
}
elsif(defined $host ){

    run_all_by_host($host, $pass);
}
else{

    run_all_by_host('ens-staging', $pass);
    run_all_by_host('ens-staging2', $pass);
}


sub run_all_by_host{

    my $host = shift;
    my $pass = shift;

    my $dbh = DBI->connect("dbi:mysql:information_schema:$host:3306", 'ensro', undef, undef);

    my $db_ext_sth = $dbh->prepare(qq[ show databases like '%var%']);

    $db_ext_sth->execute()||die;
    my $db_list = $db_ext_sth->fetchall_arrayref();
    foreach my $l(@{$db_list}){
	next if /master/;
	print "Doing $l->[0] on $host\n";
	create_table($l->[0], $host, $pass);
    }

}
    


sub create_table{

    my $db   = shift;
    my $host = shift;
    my $pass = shift;

    my $dbh = DBI->connect("dbi:mysql:$db:$host:3306", 'ensadmin', $pass, undef);


    my $len = get_size($dbh);


    $dbh->do(qq[ CREATE TABLE variation_hgvs (
                  variation_id int(11) unsigned NOT NULL,
	          hgvs_name    varchar(255) NOT NULL,
	          primary key(variation_id, hgvs_name)
                 ) ENGINE=MyISAM DEFAULT CHARSET=latin1
		 ]);

    my $trans_ins_sth = $dbh->prepare(qq[ insert ignore into variation_hgvs(
					  select vf.variation_id,
					  SUBSTR(tv.hgvs_transcript,$len,20)
					  from  variation_feature vf,
					  transcript_variation tv
					  where  vf.variation_feature_id = tv.variation_feature_id
					  and tv.hgvs_transcript like '%:c.%')
					  ]);
    
    my $prot_ins_sth  = $dbh->prepare(qq[ insert ignore into variation_hgvs(
					  select vf.variation_id,
					  SUBSTR(tv.hgvs_protein,$len,20)
					  from  variation_feature vf,
					  transcript_variation tv
					  where  vf.variation_feature_id = tv.variation_feature_id
					  and tv.hgvs_protein like '%:p.%')
					  ]);

    $trans_ins_sth->execute()||die;
    $prot_ins_sth->execute()||die;

}

## extract a transcript stable id from the target database 
## to work out the length of prefix to be removed

sub get_size{

    my $dbh = shift;
    my $tran_ext_sth = $dbh->prepare(qq[select feature_stable_id from transcript_variation limit 1]);
    $tran_ext_sth->execute()||die;
    my $trans = $tran_ext_sth->fetchall_arrayref();

    ## ENST00000509460 + ".2:"
    my $len = length($trans) + 3;
    return $len;
}



sub usage {

die "\n\tUsage generate_hgvs_table.pl -pass [admin password]

\toptions:   -host [host name]
\t           -db   [ database name]

\tBy default creates HGVS tables in all %var% databases on staging/staging2\n\n";
}     
