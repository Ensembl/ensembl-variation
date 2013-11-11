#!/usr/bin/env perl
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


### post-processing of clinvar import

## copies clinical_significance statuses to variation table
## creates a set of associated variants (Removes any previous variation_set_variation entries)
## removes empty frequency data from allele table

use strict;
use warnings;

use Getopt::Long;
use Bio::EnsEMBL::Registry;

our $SOURCENAME = 'dbSNP_ClinVar';

my $registry_file;

GetOptions ("registry_file=s"       => \$registry_file);
die "\nERROR: registry file for human db needed\n\n" unless defined $registry_file;

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all($registry_file);

my $db_adaptor = $registry->get_DBAdaptor('homo_sapiens','variation') 
    or die ("Could not get variation DBAdaptor");

my $dbh = $db_adaptor->dbc;


my $variation_ids = update_variation($dbh);

add_to_set($dbh, $variation_ids );

fix_allele($dbh);

check_phenotype_names($dbh);

delete_pheno_less($dbh);

check_counts($dbh);


## copy clinical significance to variation table from phenotype_feature 
sub update_variation{

    my $dbh = shift;
    warn "starting     update_variation\n";
    my $clinvar_ext_sth = $dbh->prepare(qq[ select variation.variation_id, phenotype_feature_attrib.value
                                            from   variation, phenotype_feature, phenotype_feature_attrib, source
                                            where source.name = ?
                                            and   phenotype_feature.source_id = source.source_id
                                            and   phenotype_feature.phenotype_feature_id = phenotype_feature_attrib.phenotype_feature_id
                                            and   phenotype_feature_attrib.attrib_type_id = 10
                                            and   variation.name = phenotype_feature.object_id
                                            ]);



    my $var_upd_sth = $dbh->prepare(qq[ update variation
                                        set clinical_significance = ?
                                        where variation_id = ?
                                       ] );



    $clinvar_ext_sth->execute( $SOURCENAME );
    my $var_list = $clinvar_ext_sth->fetchall_arrayref();

    my %assoc;  ## variants to put in set
    my %class;  ## all clinsig statuses for variation table


    foreach my $l (@{$var_list}){

	## only variants with associated statuses to go into set
	$assoc{$l->[0]} = 1  if $l->[1] =~ /pathogenic|drug-response|histocompatibility/ && $l->[1] !~ /non/;

	push @{$class{$l->[0]}}, $l->[1];
    }

    foreach my $var (keys %class){
	my @statuses = unique( @{$class{$var}});
	my $statuses = join",", @statuses;
	$var_upd_sth->execute($statuses, $var);
    }

    return \%assoc;

}
sub add_to_set{

    my ($dbh, $variation_ids ) = @_;

    my $set_id  = get_set($dbh);

    my $vsv_ins_sth = $dbh->prepare(qq[ insert into variation_set_variation
                                       (variation_id, variation_set_id)
                                        values (?,?)] );


    foreach my $var ( keys %{$variation_ids} ){
   
	$vsv_ins_sth->execute( $var, $set_id );
    }
}

## find old set id ( and remove linked variants)
##  or enter new one
sub get_set{

    my $dbh = shift;
    
    my $set_ext_sth = $dbh->prepare(qq[ select variation_set_id from variation_set where name ='clinically associated']);

    my $set_ins_sth = $dbh->prepare(qq[insert into variation_set (  name, description, short_name_attrib_id) 
                                        values ( 'clinically associated', 
                                       'Variants described by ClinVar as being probable-pathogenic, pathogenic, drug-response or histocompatibility',
                                        345) ]);

    my $set_del_sth = $dbh->prepare(qq[ delete from variation_set_variation where variation_set_id = ?]);


    ### look for old set record
    $set_ext_sth->execute();
    my $id = $set_ext_sth->fetchall_arrayref();

    if (defined $id->[0]->[0] ){
	## remove old set content
	$set_del_sth->execute($id->[0]->[0]);
	return $id->[0]->[0] ;
    }

    ### enter new set record
    $set_ins_sth->execute(); 
    $set_ext_sth->execute();
    $id = $set_ext_sth->fetchall_arrayref();

    return $id->[0]->[0] if defined $id->[0]->[0] ;


    ### give up
    die "ERROR: variation set could not be entered\n"; 

}


## remove empty frequency and population info 
sub fix_allele{

    my $dbh = shift;
    
    my $allele_upd_sth = $dbh->prepare(qq[ update allele 
                                           set frequency = \\N, population_id = \\N 
                                           where population_id = (select population_id from population where name ='clinvar')
                                           and frequency = 0
                                        ]);


    $allele_upd_sth->execute($SOURCENAME );
}

## check for unsupported characters in phenotype names 
sub check_phenotype_names{

    my $dbh = shift;

    my $pheno_ext_sth = $dbh->prepare(qq[ select description from phenotype ]);
    $pheno_ext_sth->execute()||die;

    my $ph =  $pheno_ext_sth->fetchall_arrayref();
    foreach my $l (@{$ph}){

        next unless defined $l->[0];
	my $full = $l->[0];
	$l->[0] =~ s/\w+|\-|\,|\(|\)|\s+|\/|\.|\;|\+|\'|\:|\@|\*|\%//g;

	warn "Phenotype : $full looks suspect\n" if(length($l->[0]) >0);
    }
}

sub check_counts{
 
   my $dbh = shift;

   my $var_count_ext_sth   = $dbh->prepare(qq[ select count(*) from variation where clinical_significance is not null ]);
   
   my $set_count_ext_sth   = $dbh->prepare(qq[ select count(*) from variation_set, variation_set_variation
                                               where variation_set.variation_set_id = variation_set_variation.variation_set_id
                                               and   variation_set.name ='clinically associated' ]);

   my $pheno_count_ext_sth = $dbh->prepare(qq[ select count(*) from phenotype_feature, source 
                                               where source.name = ?
                                               and source.source_id = phenotype_feature.source_id ]);

   my $featless_count_ext_sth = $dbh->prepare(qq[ select count(*) from phenotype 
                                                  where phenotype_id not in (select phenotype_id from phenotype_feature) ]);

   my $class_ext_sth       = $dbh->prepare(qq[ select value from attrib where attrib_type_id = 10 ]);

   my $class_count_ext_sth = $dbh->prepare(qq[ select clinical_significance, count(*) 
                                               from variation where clinical_significance is not null 
                                               group by clinical_significance]);


   $pheno_count_ext_sth->execute( $SOURCENAME )||die;
   my $ph =  $pheno_count_ext_sth->fetchall_arrayref();
   warn "$ph->[0]->[0] phenotype_feature entries from ClinVar\n";

   $var_count_ext_sth->execute()||die;
   my $var =  $var_count_ext_sth->fetchall_arrayref();
   warn "$var->[0]->[0] variation entries with ClinVar statuses\n";

   $set_count_ext_sth->execute()||die;
   my $set =  $set_count_ext_sth->fetchall_arrayref();
   warn "$set->[0]->[0] variation entries in ClinVar associated set\n";


   $featless_count_ext_sth->execute()||die;
   my $featless =  $featless_count_ext_sth->fetchall_arrayref();
   warn "$featless->[0]->[0] phenotype entries have no phenotype features\n";

   my %class_count;
   warn "\nGetting counts by class\n";
   $class_count_ext_sth->execute()||die;
   my $class_num = $class_count_ext_sth->fetchall_arrayref();
   foreach my $l(@{$class_num}){
       warn "$l->[1]\t$l->[0]\n";
       $class_count{$l->[0]} = $l->[1];
   }


   $class_ext_sth->execute()||die;
   my $classes = $class_ext_sth->fetchall_arrayref();
   foreach my $l(@{$classes}){
       warn "\nNo variants with class: $l->[0]\n"  unless defined  $class_count{$l->[0]} ;
   }
}

sub delete_pheno_less{

    my $dbh = shift;

    my $pheno_ext_sth   = $dbh->prepare(qq[ select phenotype_id from phenotype where description = "." ]);
    
    my $pheno_attrib_del_sth   = $dbh->prepare(qq[ delete from phenotype_feature_attrib 
                                                   where  phenotype_feature_id in
                                                    (select phenotype_feature_id from phenotype_feature 
                                                      where source_id in (select source_id  from source where name = ? )
                                                    and phenotype_id = ?)
                                               ]);
   
    my $pheno_feature_del_sth   = $dbh->prepare(qq[ delete from phenotype_feature  
                                                    where source_id in (select source_id  from source where name = ? )
                                                    and phenotype_id = ? ]);
 
    #my $pheno_del_sth   = $dbh->prepare(qq[  delete from phenotype where phenotype_id = ?  ]);

    $pheno_ext_sth->execute()||die;
    my $id =  $pheno_ext_sth->fetchall_arrayref();

   
    die "Error - 2 phenotypes called . - not cleaning up\n"  if defined $id->[1]->[0] ;
 
    die "Error - no phenotypes called . - not cleaning up\n"  unless defined $id->[0]->[0] ;

    $pheno_attrib_del_sth->execute( $SOURCENAME,  $id->[0]->[0] )||die ;
    $pheno_feature_del_sth->execute( $SOURCENAME,  $id->[0]->[0] )||die ;

}


sub unique {
  my %a;
  map { $a{$_} = 1; } @_;
  return  keys %a;
}
