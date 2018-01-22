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

## given a list of short name attribs for sets
##  - find variants which are newly added to the set
##  - find variants which have been removed from the set
##  - re-calculate the summary set list in variation_feature for all changed variants.

## quicker than recalculating everything for human

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Getopt::Long;

my ($registry_file, $updated_sets, $species);
GetOptions ("registry=s" => \$registry_file,
            "species=s"  => \$species,
	    "sets=s"     => \$updated_sets
            );

usage() unless defined $registry_file && defined $species && defined $updated_sets;

my @changed_sets = split/\,/, $updated_sets;


my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all($registry_file);

my $db_adaptor = $registry->get_DBAdaptor($species, 'variation') 
    || die "Could not get variation DBAdaptor";

my $dbh = $db_adaptor->dbc;
 

## find all variants in the set with their summary lines to check for consistency
my $get_new_sth = $dbh->prepare(qq[ select vf.variation_id, vf.variation_set_id
                                    from  variation_set_variation vsv, variation_feature vf
                                    where vsv.variation_set_id =?
                                    and vsv.variation_id = vf.variation_id
                                   ]);

## find variants where the denormalised data says is in the set
## NB. human sets with id < 10 don't change
my $get_old_sth = $dbh->prepare(qq[ select vf.variation_id, vf.variation_set_id
                                    from variation_feature vf
                                    where variation_set_id like ?
                                   ]);

my %recalculate_var;
foreach my $set (@changed_sets){

  my $set_id = get_set($dbh, $set);

  my %in_set;
  ## find everything in the set and check the set is in the summary 
  $get_new_sth->execute($set_id) ||die;;
  my $current = $get_new_sth->fetchall_arrayref();
  foreach my $l(@{$current}){
    ## save to compare to vf summary line
    $in_set{$l->[0]} = 1;
    my @summ = split/\,/,$l->[1];

    $recalculate_var{$l->[0]} = 1 unless grep $set_id,@summ;
  }
  my $total_wrong = scalar( keys %recalculate_var); 
  print "Found current: $total_wrong look wrong for set $set\n";


  ## get those with summary info and check they are still in the set
  my $pattern = '%' . $set_id . '%';
  $get_old_sth->execute($pattern) ||die;
  my $previous = $get_old_sth->fetchall_arrayref();
  foreach my $l(@{$previous}){
    $recalculate_var{$l->[0]} = 1 unless $in_set{$l->[0]} == 1;
  }
  $total_wrong = scalar( keys %recalculate_var);
  print "Found previous: $total_wrong now look wrong for set $set\n";

}

recalculate($dbh, \%recalculate_var);


## find parent sets
sub get_structure{
  my $dbh = shift;

  my $get_struc_sth = $dbh->prepare(qq[ select variation_set_sub, variation_set_super from variation_set_structure]);

  my %parent;
  $get_struc_sth->execute() ||die;
  my $dat = $get_struc_sth->fetchall_arrayref();
  foreach my $l(@{$dat}){
    $parent{$l->[0]} = $l->[1] ;
  }
  return \%parent;
}

## find set id from short name attrib  
sub get_set{
  my $dbh        = shift;
  my $short_name = shift;
 
  my $get_set_sth = $dbh->prepare(qq[ select variation_set_id
                                      from variation_set, attrib, attrib_type
                                      where attrib_type.code ='short_name'
                                      and attrib_type.attrib_type_id = attrib.attrib_type_id
                                      and attrib.value = ?
                                      and variation_set.short_name_attrib_id = attrib.attrib_id
                                    ]);

  $get_set_sth->execute($short_name) || die;
  my $set_id = $get_set_sth->fetchall_arrayref();
  die "No set found for $short_name\n" unless defined $set_id->[0]->[0];
  return $set_id->[0]->[0];

}

sub recalculate{
  my $dbh  = shift;
  my $vars = shift;

  ## get structure
  my $parent = get_structure($dbh);

  ## find all sets for a variant
  my $get_vs_sth = $dbh->prepare(qq[ select variation_set_id
                                     from variation_set_variation vsv
                                     where vsv.variation_id =?
                                   ]);

  my $vf_upd_sth = $dbh->prepare(qq[ update variation_feature set variation_set_id = ? where variation_id = ?]);


  foreach my $var (keys %{$vars}){
    my @sets;
    $get_vs_sth->execute($var)||die;
    my $dat = $get_vs_sth->fetchall_arrayref();
    foreach my $set(@{$dat}){
      push @sets, $set->[0];
      if (exists $parent->{$set->[0]}){
        push @sets, $parent->{$set->[0]};
        push @sets, $parent->{$parent->{$set->[0]}} if exists $parent->{$parent->{$set->[0]}};
      }
    }
    my $set_summary = join(",", @sets);
    $vf_upd_sth->execute( $set_summary,  $var )|| die;
  }

}

sub usage{

  die "\n\n\tUsage: update_changed_sets.pl -registry [registry file] -species [species name] -sets [short_name_attrib1,short_name_attrib2]\n\n";
}

