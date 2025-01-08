#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2025] EMBL-European Bioinformatics Institute
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

use Try::Tiny;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Utils::QCUtils qw(get_reference_base check_variant_size check_for_ambiguous_alleles check_illegal_characters);

our $DEBUG = 0;

our $SOURCENAME = 'ClinVar';

my $registry_file;

GetOptions ("registry_file=s"       => \$registry_file,
            "debug"                 => \$DEBUG);
die "\nERROR: registry file for human database needed\n\n" unless defined $registry_file;

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all($registry_file);

my $db_adaptor = $registry->get_DBAdaptor('homo_sapiens','variation') 
    or die ("Could not get variation DBAdaptor");
my $slice_adaptor = $registry->get_adaptor('homo_sapiens','core', 'slice')
        or die ("Could not get slice DBAdaptor");

my $dbh = $db_adaptor->dbc;


#do Variation QC for the variants added from ClinVar
run_variation_checks($dbh);

my $variation_ids = update_variation($dbh);

add_to_set($dbh, $variation_ids );

check_phenotype_names($dbh);

#delete_pheno_less($dbh);

check_counts($dbh);


sub run_variation_checks {
    my $dbh = shift;

    ## count pass & fail for basic reporting
    my %status;

    #get ClinVar - source_id
    my $source_id_sth = $dbh->prepare(qq[ SELECT source_id
                                    FROM source
                                    WHERE name = ?]);
    #get ClinVar inserted variants
    my $variant_ext_sth = $dbh->prepare(qq[SELECT vf.variation_id, vf.variation_name,
                                                  vf.seq_region_id, vf.seq_region_start, vf.seq_region_end, vf.seq_region_strand,
                                                  vf.source_id,
                                                  sr.name,
                                                  vf.flags, vf.allele_string
                                           FROM variation_feature vf,
                                                seq_region sr,
                                                variation v
                                          WHERE  v.source_id = ?
                                                AND vf.seq_region_id = sr.seq_region_id
                                                AND vf.variation_id = v.variation_id
                                          ]);
    my $failed_set_ext_sth = $dbh->prepare(qq[ SELECT variation_set_id
                                               FROM variation_set, attrib
                                               WHERE variation_set.short_name_attrib_id = attrib.attrib_id
                                               AND attrib.value = ? ]);
    my $insert_failed_var_sth = $dbh->prepare(qq[INSERT IGNORE INTO failed_variation(variation_id, failed_description_id)
                                                VALUES (?, ?) ]);
    my $insert_failed_var_set_var_sth = $dbh->prepare(qq[INSERT IGNORE INTO variation_set_variation (variation_id, variation_set_id)
                                              VALUES (?, ?) ]);
    my $update_failed_var_feat_sth = $dbh->prepare(qq[UPDATE variation_feature
                                                      SET variation_set_id = concat(variation_set_id,',', ?)
                                                      WHERE variation_id = ? ]);
    ## get ClinVar source id
    $source_id_sth->execute($SOURCENAME);
    my $source_id =  $source_id_sth->fetchall_arrayref();
    ## get ClinVar inserted variants
    $variant_ext_sth-> execute($source_id->[0]->[0]);
    my $var_feat_list =  $variant_ext_sth->fetchall_arrayref();

    ## look for fail_all set record
    $failed_set_ext_sth->execute( "fail_all");
    my $set_id = $failed_set_ext_sth->fetchall_arrayref();

    my %clinvar_vars = ();
    foreach my $l (@{$var_feat_list}){
        my %var = ( "start"       =>  $l->[3],
                    "end"         =>  $l->[4],
                    "strand"      =>  $l->[5],
                    "seqreg_name" =>  $l->[7],
                    "allele"      =>  $l->[9],
                    "source"      =>  $l->[6],
                    "name"        =>  $l->[1],
                    "id"          =>  $l->[0],
            );

        $clinvar_vars{$var{name}} =1;
        #check for rsID sanity:
        #the expectation is that any new rsID that has to be inserted is a new one (not in dbSNP import)
        #given that approx. 650mil variation is already known, the new rsID should have a bigger number than that
        my $number = $var{name}  =~ s/rs//r;
        warn "WARNING: clinvar rsID less than 650mil, likely typo! $var{name} \n" if ($number < 650000000 && $DEBUG == 1);

        $var{fail_reasons} = run_checks(\%var);

        $status{all}++;
        # update failed variation if needed
        if (defined $var{fail_reasons} && scalar(@{$var{fail_reasons}}) > 0 ){
            warn $var{name},"\tfailed: ", join(",",@{$var{fail_reasons}}), "\n" if $DEBUG == 1;
            ## keep a count for fail rate
            $status{fail}++ ;

            # add failed reasons
            while (my $fail = shift(@{$var{fail_reasons}}) ) {
              $insert_failed_var_sth->execute($var{id}, $fail);
            }
            # add to failed variation_set
            $insert_failed_var_set_var_sth->execute($var{id}, $set_id->[0]->[0]);
            # update variation_feature variation_set_id
            $update_failed_var_feat_sth->execute($set_id->[0]->[0], $var{id});
        } elsif($DEBUG == 1) {
            warn $var{name}, "\n";
        }
    }
    print "ClinVar imported variants: ", defined $status{all} ? $status{all} : 0, ", failed: ", defined $status{fail} ? $status{fail} : 0, "\n";
    print "ClinVar imported variant names: \n", join "\n", keys %clinvar_vars, "\n";
}

## call standard QC checks & return string of failure reasons
sub run_checks{

    my $var = shift;

    my @fail;

    # Extract reference sequence to run ref checks [ compliments for reverse strand multi-mappers]
    my $ref_seq = get_reference_base($var, $slice_adaptor, "coredb") ;

    ## Type 15: Mapped position is not compatible with reported alleles
    unless(defined $ref_seq){
        ## don't check further if obvious coordinate error
        push @fail, 15;
        return \@fail;
    }
    ## Type 2: None of the variant alleles match the reference allele - is ref base in agreement?
    my $exp_ref = (split/\//, $var->{allele} )[0] ;
    push @fail, 2 unless (defined $exp_ref && "\U$exp_ref" eq "\U$ref_seq"); ## using soft masked seq

    ## is either allele of compatible length with given coordinates?
    my $ref = (split/\//, $var->{allele} )[0];
    my $match_coord_length = check_variant_size( $var->{start}, $var->{end}, $ref);
    push @fail, 15 unless  ($match_coord_length == 1);

    ## Type 14 resolve ambiguities before reference check - flag variants & alleles as fails
    my $is_ambiguous = check_for_ambiguous_alleles( $var->{allele} );
    push @fail, 14  if(defined $is_ambiguous) ;

    ## Type 13 Alleles contain non-nucleotide characters
    my $is_illegal = check_illegal_characters( $var->{allele} );
    push @fail, 13  if(defined $is_illegal->[0]) ;

    return \@fail;
}

## copy clinical significance to variation table from phenotype_feature 
sub update_variation{
  my $dbh = shift;

  print "Starting update_variation...\n" if $DEBUG == 1;

  my $clinvar_ext_sth = $dbh->prepare(qq[ select variation.variation_id, phenotype_feature_attrib.value
                                            from   variation, phenotype_feature, phenotype_feature_attrib, source, attrib_type
                                            where source.name = ?
                                            and phenotype_feature.source_id = source.source_id
                                            and phenotype_feature.phenotype_feature_id = phenotype_feature_attrib.phenotype_feature_id
                                            and phenotype_feature_attrib.attrib_type_id = attrib_type.attrib_type_id 
                                            and attrib_type.code = 'clinvar_clin_sig'
                                            and variation.name = phenotype_feature.object_id
                                            ]);

  my $var_upd_sth = $dbh->prepare(qq[ update variation
                                        set clinical_significance = ?
                                        where variation_id = ?
                                     ] );

  my $varf_upd_sth = $dbh->prepare(qq[ update variation_feature
                                         set clinical_significance = ?
                                         where variation_id = ?
                                     ] );

  $clinvar_ext_sth->execute( $SOURCENAME );
  my $var_list = $clinvar_ext_sth->fetchall_arrayref();

  my %assoc;  ## variants to put in set
  my %class;  ## all clinsig statuses for variation table

  foreach my $l (@{$var_list}){
    ## save ids to make 2 variation_sets - 'All' and 'Clinically significant'
    $assoc{$l->[0]}{A} = 1;
    ## only variants with associated statuses to go into set
    $assoc{$l->[0]}{C} = 1  if $l->[1] =~ /pathogenic|drug-response|histocompatibility/i && $l->[1] !~ /non/;

    $l->[1]  =~ s/\//\,/g ; # convert 'pathogenic/likely pathogenic' to 'pathogenic,likely pathogenic'
    $l->[1]  =~ s/\;\s+/\,/g ; # convert 'uncertain significance; risk factor' to 'uncertain significance,risk factor'
    $l->[1]  =~ s/\,\s+/\,/g ; # convert 'likely benign, other' to 'likely benign,other' or 'benign, association, risk factor' to 'benign,association,risk factor'
    #replace 'conflicting interpretations of pathogenicity' with 'uncertain significance'
    #for the purpose of variation, variation_feature clin_sig entry
    $l->[1]  =~ s/conflicting interpretations of pathogenicity/uncertain significance/g;
    #similar replace of 'association not found' to 'other'
    $l->[1]  =~ s/association not found/other/g;

    # remove commas from some values
    $l->[1]  =~ s/pathogenic,low penetrance/pathogenic low penetrance/g;

    push @{$class{$l->[0]}}, $l->[1];
  }

  foreach my $var (keys %class){
    my @statuses = unique( @{$class{$var}});
    my $statuses = join",", @statuses;
      try {
        $var_upd_sth->execute($statuses, $var);
        $varf_upd_sth->execute($statuses, $var);
      } catch {
        warn "variation_id:$var<>statuses:$statuses<\n";
        warn "issue with:@_\n";
        next;
      }
  }

  print "Starting update_variation... done\n" if $DEBUG == 1;

  return \%assoc;
}

sub add_to_set{

    my ($dbh, $variation_ids ) = @_;

    my $clin_set_id  = get_set($dbh, 'clin_assoc');
    my $all_set_id   = get_set($dbh, 'ClinVar');

    my $vsv_ins_sth = $dbh->prepare(qq[ insert into variation_set_variation
                                       (variation_id, variation_set_id)
                                        values (?,?)] );


    foreach my $var ( keys %{$variation_ids} ){
    ## update 'all clinvar' set
    $vsv_ins_sth->execute( $var, $all_set_id );
    ## update 'clinically associated' set
    $vsv_ins_sth->execute( $var, $clin_set_id ) 
      if defined $variation_ids->{$var}{C} && $variation_ids->{$var}{C} == 1;
    }
}

## find old set id ( and remove linked variants)
##  or enter new one
sub get_set{

    my $dbh = shift;
    my $set = shift; ## short attrib name for set

    ## info on sets
    my $data =      {

   'clin_assoc' => {
      'desc' => 'Variants described by ClinVar as being probable-pathogenic, pathogenic, drug-response or histocompatibility',
      'name' => 'Clinically associated'},

      'ClinVar'           => {
      'desc' => 'Variants with ClinVar annotation',
      'name' => 'All ClinVar'},
    };

    my $attrib_ext_sth = $dbh->prepare(qq[ select attrib_id from attrib 
                                           where value = ? and attrib_type_id = 477]);
    
    my $set_ext_sth = $dbh->prepare(qq[ select variation_set_id 
                                        from variation_set, attrib
                                        where variation_set.short_name_attrib_id = attrib.attrib_id
                                        and attrib.value =? ]);

    my $set_ins_sth = $dbh->prepare(qq[ insert into variation_set 
                                       (name, description, short_name_attrib_id) 
                                        values ( ?,?,?  ) ]);

    my $set_del_sth = $dbh->prepare(qq[ delete from variation_set_variation where variation_set_id = ?]);


    ### look for old set record
    $set_ext_sth->execute( $set);
    my $set_id = $set_ext_sth->fetchall_arrayref();

    if (defined $set_id->[0]->[0] ){
      ## remove old set content
      $set_del_sth->execute($set_id->[0]->[0]);
      return $set_id->[0]->[0] ;
    }

    ### enter new set record
    $attrib_ext_sth->execute( $set);
    my $attrib_id = $attrib_ext_sth->fetchall_arrayref();
    die "Exiting: attrib not available for $set\n" unless defined $attrib_id->[0]->[0];

    $set_ins_sth->execute( $data->{$set}->{name}, $data->{$set}->{desc}, $attrib_id->[0]->[0] ); 
    $set_id = $dbh->db_handle->last_insert_id(undef, undef, qw(variation_set set_id)) or die "no insert id for set $set\n";

    return $set_id;
}



## check for unsupported characters in phenotype names 
sub check_phenotype_names{
  my $dbh = shift;

  my $pheno_ext_sth = $dbh->prepare(qq[ select phenotype_id, description from phenotype ]);
  $pheno_ext_sth->execute() or die;

  my $ph =  $pheno_ext_sth->fetchall_arrayref();
  foreach my $l (@{$ph}){
    if (!defined $l->[1]) {
      warn "Phenotype id:$l->[0] has no description\n";
    }

    my $full = $l->[1];
    $l->[1] =~ s/\w+|\-|\,|\(|\)|\s+|\/|\.|\;|\+|\'|\:|\@|\*|\%//g;

    if(length($l->[1]) > 0 && $DEBUG == 1) {
      warn "Phenotype : $full (id:$l->[0]) looks suspect\n";
    }
  }
}

## basic report on imported data
sub check_counts{
 
   my $dbh = shift;

   print "*** Imported data report ***\n";

   my $var_count_ext_sth   = $dbh->prepare(qq[ select count(*) from variation 
                                               where clinical_significance is not null 
                                               ]);

   my $varf_count_ext_sth  = $dbh->prepare(qq[ select count(*) from variation_feature 
                                               where clinical_significance is not null 
                                               ]);
   
   my $set_count_ext_sth   = $dbh->prepare(qq[ select count(*) from variation_set, variation_set_variation
                                               where variation_set.variation_set_id = variation_set_variation.variation_set_id
                                               and   variation_set.name = ? ]);

   my $pheno_count_ext_sth = $dbh->prepare(qq[ select count(*) from phenotype_feature, source 
                                               where source.name = ?
                                               and source.source_id = phenotype_feature.source_id ]);

   my $featless_count_ext_sth = $dbh->prepare(qq[ select count(*) from phenotype 
                                                  where phenotype_id not in (select phenotype_id from phenotype_feature) ]);

   my $class_ext_sth       = $dbh->prepare(qq[ select attrib.value 
                                               from attrib, attrib_type
                                               where attrib.attrib_type_id = attrib_type.attrib_type_id 
                                               and attrib_type.code = 'clinvar_clin_sig' ]);

   my $class_count_ext_sth = $dbh->prepare(qq[ select clinical_significance, count(*) 
                                               from variation where clinical_significance is not null 
                                               group by clinical_significance]);

  my $pheno_attrib_ext_sth = $dbh->prepare(qq[ SELECT a.code, count(*)
                                               FROM attrib_type a, phenotype_feature_attrib pfa, phenotype_feature pf, source s
                                               WHERE s.name='ClinVar'
                                               AND pf.source_id = s.source_id AND pf.phenotype_feature_id = pfa.phenotype_feature_id AND pfa.attrib_type_id = a.attrib_type_id
                                               group by a.code]);

  my $omim_set_ext_sth = $dbh->prepare(qq[ SELECT count(*)
                                           FROM variation_set_variation vsv, variation_set vs, attrib a
                                           WHERE vs.short_name_attrib_id = a.attrib_id
                                           AND vsv.variation_set_id = vs.variation_set_id
                                           AND a.value = ? ]);

  $pheno_count_ext_sth->execute( $SOURCENAME ) or die;
  my $ph =  $pheno_count_ext_sth->fetchall_arrayref();
  print "$ph->[0]->[0] phenotype_feature entries from ClinVar\n";

  $var_count_ext_sth->execute() or die;
  my $var =  $var_count_ext_sth->fetchall_arrayref();
  print "$var->[0]->[0] variation entries with ClinVar statuses\n";

  $varf_count_ext_sth->execute() or die;
  my $varf =  $varf_count_ext_sth->fetchall_arrayref();
  print "$varf->[0]->[0] variation_feature entries with ClinVar statuses\n";

  foreach my $set ('Clinically associated variants', 'All ClinVar'){
    $set_count_ext_sth->execute($set) or die;
    my $set_count =  $set_count_ext_sth->fetchall_arrayref();
    print "$set_count->[0]->[0] variation entries in $set set\n";
  }

  $featless_count_ext_sth->execute() or die;
  my $featless =  $featless_count_ext_sth->fetchall_arrayref();
  print "$featless->[0]->[0] phenotype entries have no phenotype features\n";

  my %class_count;
  print "\nGetting counts by class\n";
  $class_count_ext_sth->execute() or die;
  my $class_num = $class_count_ext_sth->fetchall_arrayref();
  foreach my $l(@{$class_num}){
    print "$l->[1]\t$l->[0]\n";
    $class_count{$l->[0]} = $l->[1];
  }

  $class_ext_sth->execute() or die;
  my $classes = $class_ext_sth->fetchall_arrayref();
  foreach my $l(@{$classes}){
    print "\nNo variants with class: $l->[0]\n"  unless defined $class_count{$l->[0]} ;
  }

  #get OMIM variant set count:
  $omim_set_ext_sth->execute( "ph_omim"); #get OMIM set id
  my $omim_set_count = $omim_set_ext_sth->fetchrow_array();
  print "OMIM set variants: ", $omim_set_count, "\n";

  print "Getting ClinVar phenotype_attrib counts\n";
  $pheno_attrib_ext_sth->execute() or die;
  my $rows = $pheno_attrib_ext_sth->fetchall_arrayref();
  foreach my $row(@{$rows}) {
    print join(': ', @{$row}), "\n";
  }

  print "******\n";
}

## check for and remove 'not provided' phenotypes
sub delete_pheno_less{

    my $dbh = shift;

    my $pheno_ext_sth   = $dbh->prepare(qq[ select phenotype_id from phenotype where description = "ClinVar: phenotype not specified" ]);
    
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

    $pheno_ext_sth->execute() or die;
    my $id =  $pheno_ext_sth->fetchall_arrayref();

   
    die "Error - 2 phenotypes called . - not cleaning up\n"  if defined $id->[1]->[0] ;
 
    die "Error - no phenotypes called . - not cleaning up\n"  unless defined $id->[0]->[0] ;

    $pheno_attrib_del_sth->execute( $SOURCENAME,  $id->[0]->[0] ) or die ;
    $pheno_feature_del_sth->execute( $SOURCENAME,  $id->[0]->[0] ) or die ;

}


sub unique {
  my %a;
  map { $a{$_} = 1; } @_;
  return  keys %a;
}
