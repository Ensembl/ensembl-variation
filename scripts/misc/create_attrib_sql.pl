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

# Read a variation config module and generate some SQL to populate the attrib, 
# attrib_type and attrib_set tables in the variation database

use strict;
use warnings;

use DBI;
use Getopt::Long;

my $config;
my $no_model;
my $host;
my $port;
my $user;
my $pass;
my $db;
my $help;

GetOptions(
    "config=s"  => \$config,
    "no_model"  => \$no_model,
    "host=s"    => \$host,
    "port=i"    => \$port,
    "user=s"    => \$user,
    "pass=s"    => \$pass,
    "db=s"      => \$db,
    "help|h"    => \$help,
);

unless ($no_model || ($host && $user && $db)) {
    print "Missing required parameter...\n" unless $help;
    $help = 1;
}

if ($help) {
    print "Usage: $0 --config <module> --host <host> --port <port> --user <user> --pass <pass> --db <database> --no_model --help > attrib_entries.sql\n";
    exit(0);
}

# pull in our config module

$config ||= 'Bio::EnsEMBL::Variation::Utils::Config';

eval qq{require $config};

die "Failed to require config module '$config':\n$@" if $@;

# and import the variables we need

our $MAX_ATTRIB_CODE_LENGTH;
our @ATTRIB_TYPES;
our @ATTRIB_SETS;
our %ATTRIBS;

eval {
    $config->import(qw(
        $MAX_ATTRIB_CODE_LENGTH
        @ATTRIB_TYPES
        @ATTRIB_SETS
        %ATTRIBS
    ));
};

die "Failed to import required data structures from config module '$config':\n$@" if $@;

# format strings for inserting into our 3 tables

my $attrib_type_fmt = 
    q{INSERT IGNORE INTO attrib_type (attrib_type_id, code, name, description) VALUES (%d, %s, %s, %s);};
my $attrib_fmt = 
    q{INSERT IGNORE INTO attrib (attrib_id, attrib_type_id, value) VALUES (%d, %d, '%s');};
my $set_fmt = 
    q{INSERT IGNORE INTO attrib_set (attrib_set_id, attrib_id) VALUES (%d, %d);};

# these hashes store mappings to our attrib and attrib_type IDs

my %attrib_ids;
my $attrib_type_ids;

# these variables store the current highest used ID for each table 

my $last_attrib_type_id = 0;
my $last_attrib_id      = 0;
my $last_attrib_set_id  = 0;

# these hashes store the existing IDs from the database, if --no_model
# is used then these will just be empty and new IDs will be generated

my $existing_attrib_type;
my $existing_attrib;
my $existing_set;

unless ($no_model) {

    # prefetch existing IDs from the database

    my $dbh = DBI->connect(
        "DBI:mysql:database=$db;host=$host;port=$port",
        $user,
        $pass,
    );

    my $get_types_sth = $dbh->prepare(qq{
        SELECT code, attrib_type_id FROM attrib_type
    });

    $get_types_sth->execute;

    while (my ($code, $id) = $get_types_sth->fetchrow_array) {
        $existing_attrib_type->{$code} = $id;
        $last_attrib_type_id = $id if $id > $last_attrib_type_id;
    }

    my $get_attribs_sth = $dbh->prepare(qq{
        SELECT attrib_type_id, value, attrib_id FROM attrib
    });

    $get_attribs_sth->execute;

    while (my ($type_id, $value, $id) = $get_attribs_sth->fetchrow_array) {
        $existing_attrib->{$type_id}->{$value} = $id;
        $last_attrib_id = $id if $id > $last_attrib_id;
    }

    my $get_sets_sth = $dbh->prepare(qq{
        SELECT attrib_set_id, attrib_id FROM attrib_set
    });

    $get_sets_sth->execute;

    while (my ($set_id, $attrib_id) = $get_sets_sth->fetchrow_array) {
        $existing_set->{$set_id}->{$attrib_id} = 1;

        $last_attrib_set_id = $set_id if $set_id > $last_attrib_set_id;
    }
}

# the following subroutines are used to get the corresponding IDs for
# our tables, they will return existing IDs where possible and generate
# new ones when required

sub get_attrib_type_id {
    my ($code) = @_;

    warn "$code is > $MAX_ATTRIB_CODE_LENGTH characters, have you changed the schema to match?"
        if length($code) > $MAX_ATTRIB_CODE_LENGTH;

    my $id = $existing_attrib_type->{$code};

    unless (defined $id) {
        $id = ++$last_attrib_type_id;
        $existing_attrib_type->{$code} = $id
    }

    return $id;
}

sub get_attrib_id {
    my ($type_id, $value) = @_;

    my $id = $existing_attrib->{$type_id}->{$value};

    unless (defined $id) {
        $id = ++$last_attrib_id;
        $existing_attrib->{$type_id}->{$value} = $id
    }

    return $id;
}

sub get_attrib_set_id {

    my $new_set = { map {$_ => 1} @_ };

    # we need to check if the new set is a sub or super set 
    # of an existing set, i.e. it is just adding or removing
    # members of the set, in which case we can reuse the set id
    # otherwise we assign a new id

    my $is_subset = sub {
        my ($s1, $s2) = @_;
        for my $e (keys %$s2) {
            return 0 unless $s1->{$e};
        }
        return 1;
    };

    for my $set_id (keys %$existing_set) {
        my $set = $existing_set->{$set_id};
        if ($is_subset->($set, $new_set) || $is_subset->($new_set, $set)) {
            return $set_id;
        }
    }

    # assigne a nee set id
  
    $last_attrib_set_id++;

    map { $existing_set->{$last_attrib_set_id}->{$_} = 1 } keys %$new_set;

    return $last_attrib_set_id;
}

# the SQL string we are building

my $SQL;

# first define the attrib type entries

for my $attrib_type (@ATTRIB_TYPES) {
    
    my $code        = delete $attrib_type->{code} or die "code required for attrib_type";
    my $name        = delete $attrib_type->{name};
    my $description = delete $attrib_type->{description};
    
    my $attrib_type_id = get_attrib_type_id($code);

    die "Unexpected entries in attrib_type definition: ".(join ',', keys %$attrib_type)
        if keys %$attrib_type;

    $SQL .= sprintf($attrib_type_fmt, 
        $attrib_type_id, 
        "'$code'",
        ($name ? "'$name'" : "''"),
        ($description ? "'$description'" : 'NULL'),
    )."\n";

    $attrib_type_ids->{$code} = $attrib_type_id;
}

# second, take the entries from the ATTRIBS and add them as single-element hashes to the @ATTRIB_SETS array
while (my ($type,$values) = each(%ATTRIBS)) {
    
    map {push(@ATTRIB_SETS,{$type => $_})} @{$values};
} 

# third, loop over the ATTRIB_SETS array and add attribs and assign them to attrib_sets as necessary
for my $set (@ATTRIB_SETS) {
    
    # Keep the attrib_ids
    my @attr_ids;
    
    # Iterate over the type => value entries in the set
    while (my ($type,$value) = each(%{$set})) {
        
        # Lookup the attrib_type
        my $attrib_type_id = $attrib_type_ids->{$type} or next;
        
        # insert a new attrib if we haven't seen it before
        my $attrib_id = $attrib_ids{$type . "_" . $value};
        
        unless (defined($attrib_id)) {
            $attrib_id = get_attrib_id($attrib_type_id, $value);
            $SQL .= sprintf($attrib_fmt, $attrib_id, $attrib_type_id, $value)."\n"; 
            $attrib_ids{$type . "_" . $value} = $attrib_id;
        }
        
        push(@attr_ids,$attrib_id);   
    }
    
    # If the set had more than one attribute, group them into a set
    if (scalar(@attr_ids) > 1) {
        
        my $attrib_set_id = get_attrib_set_id(@attr_ids);
        map {$SQL .= sprintf($set_fmt, $attrib_set_id, $_)."\n"} @attr_ids;
        
    }
}

# print out our SQL

print $SQL . "\n" if $SQL;


