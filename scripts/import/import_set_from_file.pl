#!/usr/bin/perl


use strict;
use warnings;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use File::Spec;

=head1 import_set_from_file.pl 

populate variation_set_variation from a list of variants & pre-loaded variation_set  

import_set_from_file.pl -load_file [rs list] -variation_set [short name attrib] -species [species] -registry [ensembl registry]  -done_file [previously loaded rs id]

sql to pre-load set & attrib if not present:

insert into attrib (attrib_id, attrib_type_id, value ) values (332,477,'Affy_500K');

insert into variation_set (name, description, short_name_attrib_id) values ( 'Affy GeneChip 500K Array','Variants from the Affymetrix GeneChip Human Mapping 500K Array Set',332 );

=cut

my ($load_file, $registry_file, $species, $set_attrib, $done_file);

GetOptions ("load_file=s"     => \$load_file,
            "registry=s"      => \$registry_file,
	    "species=s"       => \$species,
 	    "variation_set=s" => \$set_attrib,
	    "done_file=s"     => \$done_file
    );

unless(defined $load_file && defined $set_attrib){
die "\n Usage: import_set_from_file.pl -load_file [rs list] -variation_set [short name attrib] -species [species] -registry [ensembl registry]

\tOptions: -done_file [previously loaded rs id]\n\n";
}

## read variants already entered if recovering from a crash
my $done = load_done($done_file) if $done_file;

my $registry = 'Bio::EnsEMBL::Registry';
$registry_file = File::Spec->rel2abs($registry_file);
$registry->load_all($registry_file);

my $dba  = $registry->get_DBAdaptor($species, 'variation');
my $dbh  = $dba->dbc->db_handle();

## quicker than API look up
my $var_ext_sth         =  $dbh->prepare(qq [select variation_id from variation where name =?]);

my $old_var_ext_sth     =  $dbh->prepare(qq [select variation_id
                                             from variation_synonym
                                             where name = ?
                                            ]);

## Using ignore not ideal but  variants may be in the list twice, under primary name and synonym
my $insert_varset_var   =  $dbh->prepare(qq [insert ignore into variation_set_variation
		  	                     (variation_set_id, variation_id)
				             values (?, ?)
                                            ]); 

## find pre-existing set
my $vs_adaptor = $registry->get_adaptor($species,'variation','variationset');
my $varset = $vs_adaptor->fetch_by_short_name($set_attrib);
die "No set found for $set_attrib\n\n" unless defined $varset; 

my $load_var;
if($load_file =~ /\.gz$/){
  open $load_var, "gunzip -c $load_file |" ||die "Failed to open $load_file :$!\n";
}
else{
  open $load_var, $load_file ||die "Failed to open $load_file :$!\n";
}

while(<$load_var>){
  chomp;
  my $rs = (split)[0];

  ##look up on current id
  $var_ext_sth->execute(  $rs ) ||die "Error geting var id $rs\n";
  my $id = $var_ext_sth->fetchall_arrayref();

  unless(defined $id->[0]->[0]){ 
    ##look up on previous id
    $old_var_ext_sth->execute($rs) ||die "Error geting synon var id $rs\n";
    $id = $old_var_ext_sth->fetchall_arrayref();
    
    unless(defined $id->[0]->[0]){ 
      warn "Could not find id for $rs\n";
      next;
    }
  }  
  $insert_varset_var->execute( $varset->dbID() , $id->[0]->[0] ) ||die "Error inserting var_set var $rs\n"; 
}

## read variants already entered if recovering from a crash
sub load_done{

  my $done = shift;
  my %done;

  open my $done_list, $done ||die "Failed to open $done :$!\n";
  while(<$done_list>){
    my $rs = (split)[0];
    $done{$rs} = 1;
  }
  return \%done;
}


