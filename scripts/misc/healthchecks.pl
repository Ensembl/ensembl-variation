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

=head healthchecks.pl

Joining human tables to check the display status is the same across variation, variation_feature
and transcript_variation is inefficient for human databases.

This checks in a less databases intensive way.

=cut

use strict;
use warnings;
use DBI;
use Getopt::Long;

## less common value to check on
our %value = ( display  => 0,
               somatic  => 1 
             );

my ($db, $host, $port, $user, $mode);

GetOptions ("db=s"    => \$db,
            "host=s"  => \$host,
            "port=s"  => \$port,
            "user=s"  => \$user,
            "mode:s"  => \$mode,
    );

die usage() unless defined $host && defined $user ;


$port ||= 3306;


my $databases;
if( defined $db){
    push @{$databases}, $db;
}
else{
    ## find all variation databases on the host
    $databases = get_dbs_by_host($host, $port, $user );
} 

if($mode =~/display/){
    check_denorm($databases, 'display');
}
elsif($mode =~/somatic/){
    check_denorm($databases, 'somatic');
}
elsif($mode =~/both/){
    check_denorm($databases, 'display');
    check_denorm($databases, 'somatic');
}
else{
    warn "\nERROR: Mode needed\n\n";
    die usage();
}


sub check_denorm{

  my $databases = shift;
  my $column    = shift;

  foreach my $db_name (@{$databases}){
    my $dbh = DBI->connect( "dbi:mysql:$db_name\:$host\:$port", $user, undef, undef) ||die "Failed to connect to $db_name\n";

    ## get the set with the least common value for each table
    my $var_ext_sth  = $dbh->prepare(qq[ select variation_id 
                                        from variation where $column = $value{$column} ],
                                     {mysql_use_result => 1});

    my $varf_ext_sth = $dbh->prepare(qq[ select variation_id, variation_feature_id 
                                         from variation_feature where $column =  $value{$column}],
                                     {mysql_use_result => 1});

    my $trv_ext_sth  = $dbh->prepare(qq[ select variation_feature_id 
                                         from transcript_variation where $column =  $value{$column}],
                                     {mysql_use_result => 1});


    ## check if variation_feature or transcript_variation records exist 
    ##if the value is not found to be the same as in variation
    my $varf_check_sth = $dbh->prepare(qq[ select variation_feature_id, $column from variation_feature where variation_id = ?]);
    my $trv_check_sth  = $dbh->prepare(qq[ select $column from transcript_variation where variation_feature_id = ?]);


    ## get all variation's with least common value
    my %variation;
    $var_ext_sth->execute()||die;
    my $var = $var_ext_sth->fetchall_arrayref();
    foreach my $l(@{$var}){
      $variation{$l->[0]} = 1;
    }


    ## get all variation_feature's with least common value to compare
    my %variation_feature;
    my %checked_var;
    $varf_ext_sth->execute()||die;
    my $varf = $varf_ext_sth->fetchall_arrayref();

    foreach my $l(@{$varf}){

      ## report if different
      print "$db_name variation_feature.$column different to variation.$column : $l->[0]\n" 
        unless defined $variation{$l->[0]};

      ## save for reverse check
      $checked_var{$l->[0]} = 1;

      ## save VF id to check TV
      $variation_feature{$l->[1]} = 1;
    }


    ## if both variation and variation_feature do not have the least common value
    ## either there is no variation_feature or there is and error
    foreach my $var(keys %variation){

      ## already checked - skip
      next if $checked_var{$var} == 1;

      $varf_check_sth->execute($var);
      my $unchecked_var = $varf_check_sth->fetchall_arrayref();
      print "$db_name variation_feature.$column different to variation.$column : $var\n" 
        if defined $unchecked_var->[0]->[0];
    }


    ## delete to save memory
    undef %variation;
    undef %checked_var;


    ## Check transcript variation
    ## get all transcript_variation's with least common value to compare
    my %checked_varf;
    $trv_ext_sth->execute()||die;
    my $trv = $trv_ext_sth->fetchall_arrayref();
    foreach my $l(@{$trv}){

      ## should not be least common value if variation_feature not least common value
      print "$db_name transcript_variation.$column different to variation_feature.$column : $l->[0]\n" 
        unless defined $variation_feature{$l->[0]};

      ## save for reverse check
      $checked_varf{$l->[0]} = 1;
    }


    ## if both transcript_variation and variation_feature do not have the least common value
    ## either there is no transcript_variation or an error
    foreach my $varf(keys %variation_feature){

      ## already checked - skip
      next if $checked_varf{$varf} == 1;

      $trv_check_sth->execute($varf);
      my $unchecked_varf = $trv_check_sth->fetchall_arrayref();
      print "$db_name transcript_variation.$column different to variation_feature.$column : $varf\n" 
        if defined $unchecked_varf->[0]->[0];
    }
  }
}

sub get_dbs_by_host{

    my ($host, $port, $user ) = @_;

    my @databases;

    my $dbh = DBI->connect("dbi:mysql:information_schema:$host:$port", $user, undef, undef) || die "Failed to look up available databases\n";

    my $db_ext_sth = $dbh->prepare(qq[ show databases like '%variation%']);

    $db_ext_sth->execute()||die;
    my $db_list = $db_ext_sth->fetchall_arrayref();
    foreach my $l(@{$db_list}){
        next if $l->[0] =~/master/;
        push @databases, $l->[0] ;
    }

    return \@databases;
}


sub usage{

  die "\n\thealthchecks.pl -host [host] 
                        -user [read-user name] 
                        -mode [display|somatic|both]\n

\t\tOptions: -db [database name]    default: all variation databases on the host
\n\n";
}

