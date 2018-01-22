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

=head1 import_Phencode.pl

   - extracts data from local phencode database
        - variants with multiple or inconsistent positions are skipped 
   - runs QC
        - bad alleles flagged
        - 4-allele substitutions flagged
        - reference miss-matches flagged
        - inconsistent coordinates flagged
   - imports data to ensembl schema
        - variants =< 50 bases are variations
        - variants > 50 bases are structural variations
        - variants already imported from dbSNP with Phencode ids are not re-entered
   - creates variation set


requirements in ensembl db:

    - loaded attribs for fail statuses
    - seq_region table
    - dbSNP load with Phencode synonyms for redundancy checking

=cut


use strict;
use warnings;
use Getopt::Long;

use Bio::DB::Fasta;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp );
use Bio::EnsEMBL::Variation::Utils::Sequence qw( get_hgvs_alleles);
use Bio::EnsEMBL::Variation::Utils::QCUtils qw( check_four_bases get_reference_base check_illegal_characters check_for_ambiguous_alleles remove_ambiguous_alleles find_ambiguous_alleles check_variant_size);

my ($registry_file,  $fasta_file,  $phenu, $phenp);

GetOptions ( "fasta=s"      => \$fasta_file,
             "registry=s"   => \$registry_file,
             "phenu=s"      => \$phenu,
             "phenp=s"      => \$phenp,

    );

usage() unless defined $registry_file && defined $fasta_file && defined $phenu && defined $phenp;

### write to and read from tmp file
our $DATA_FILE = "phencode_data_QC.txt";

my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_all($registry_file);

my $varfeat_adaptor   = $reg->get_adaptor('homo_sapiens', 'variation', 'variationfeature');
my $variation_adaptor = $reg->get_adaptor('homo_sapiens', 'variation', 'variation');
my $allele_adaptor    = $reg->get_adaptor('homo_sapiens', 'variation', 'allele');


## create index on fasta for reference check or reference allele determination
## unless one is supplied
my $db = Bio::DB::Fasta->new( $fasta_file );
  
## writes all data to a file to allow checking & returns the number of times each name is seen
my $var_counts = extract_data($phenu, $phenp);  

## find or enter source info
my %source_data  = ( "name"    => "PhenCode",
                     "version" => "14-Nov-2012",
                     "desc"    => "PhenCode is a collaborative project to better understand the relationship between genotype and phenotype in humans",
                     "url"     => "http://phencode.bx.psu.edu/",
                     "somatic" => 0
    );

my $source_ob = get_source($reg,  \%source_data );

## get seq_region_ids for quick variation_feature creation
my $seq_ids   = get_seq_ids($varfeat_adaptor->dbc );


## get names of PhenCode variants already imported from dbSNP
my $already_loaded = get_dbSNP_data($varfeat_adaptor->dbc );


## save variant ids to add to set
my %new_var;
my %new_svar;


## read the file just written  & enter data 
open my $infile, $DATA_FILE ||die "Failed to open data file to read : $!\n";
while(<$infile>){

    chomp;

    my @a = split/\t/;

    ## ~40 with same name, different locations (in repeat) - skip these
    next if  $var_counts->{$a[0]} > 1;

    ## skip if already imported via dbSNP
    next if $already_loaded->{$a[0]};

    my $ref_len = $a[4] - $a[3];
    my $al_len  = length($a[1]);

    if( $ref_len> 50 || $al_len  > 52){
        ## long deletions to go in structural variation table for display purposes
        my $dbID = insert_svar($varfeat_adaptor->dbc, \@a, $source_ob ,$seq_ids );
        ## save structural var ids to add to variation_set later
        $new_svar{$dbID} = 1;
    }
    else{
        my $dbID = insert_var(\@a);
        ## save var ids to add to variation_set later
        $new_var{$dbID} = 1;

    }
}


## add all variants & structural variants to phencode set
add_to_set( $varfeat_adaptor->dbc, \%new_var, \%new_svar);


sub insert_var{

    my $line = shift;

    if ($line->[5] eq "-1" && $line->[1] !~/Phencode|del|ins/i){

        ## compliment individual alleles keeping ref/alt order
        my @al = split(/\//,$line->[1]);
        $line->[1] = '';
        foreach my $al(@al){
            reverse_comp(\$al) ;
            $line->[1] .= $al . "/";
        }
        $line->[1] =~ s/\/$//;
        $line->[5] = 1;
    }

    my $var = enter_var($line, 
                        $variation_adaptor,
                        $source_ob 
        );
   
    enter_varfeat($line, 
                  $varfeat_adaptor,
                  $source_ob,
                  $var,
                  $seq_ids
        );
    
    enter_alleles($line,
                  $allele_adaptor,
                  $var
        );

    return $var->dbID();
  
}



## export data from local Phencode database and run QC checks
## write tmp file for checking
sub extract_data{

    my $phenu = shift;
    my $phenp = shift;

    ## save the number of times a variant name is see & don't import those seen twice ( position uncertain)
    my %count;

    ## count pass & fail for basic reporting
    my %status; 

    my $dbh = DBI->connect('dbi:mysql:phencode:ens-variation2:3306', $phenu, $phenp, undef);
    
    my $variant_ext_sth = $dbh->prepare(qq[ SELECT gv.id,
                                              gv.name,
                                              map.label,
                                              map.chrom,
                                              map.chromStart,
                                              map.chromEnd,
                                              map.strand,
                                              gv.srcId                                                 
                                        FROM gv, gvPosHg19 map 
                                        WHERE gv.id = map.name
                                        ]);       
    

    open my $out, ">$DATA_FILE"|| die "Failed to open data file to write: $!\n";
   
    $variant_ext_sth->execute() ||die;
    my $data = $variant_ext_sth->fetchall_arrayref();
    
    
    foreach my $l (@{$data}){
        
        ## skip unless hgvs in DNA terms  
        next unless $l->[2] =~ /\:g|\:c/;
        
        $l->[3] =~ s/chr//;
        $l->[4]++;
        
        ## hgvs but on non-reference seq => switch seq name & coords
        my ($reported_seq,$change) = split(/\:g|\:c\./,$l->[2], 2);
        
        ## sort out alleles 
        my ($ref_allele, $alt_allele);
        
        eval{
	    ## this expects dups to have alleles declared
            ($ref_allele, $alt_allele) = get_hgvs_alleles( $l->[2] );
        };
        if($@){
            warn "skipping $l->[2]  - $@\n";
            next;
        }
        ## inserted bases may not be supplied
        $alt_allele = "PhenCode_variation" if $change =~/ins/ && ! $alt_allele ;
        
        ## take deletion from reference
        if ($change =~/del/  && ! $ref_allele ){
            $ref_allele =  $db->seq($l->[3], $l->[4], $l->[5]) ;
            reverse_comp(\$ref_allele )if $l->[6] eq "-";
        }
	if ($change =~/dup/  && ! $alt_allele ){
            my $dup_allele =  $db->seq($l->[3], $l->[4], $l->[5]) ;
            reverse_comp(\$ref_allele )if $l->[6] eq "-";
	    $alt_allele = "dup"  . $dup_allele
        }
        unless (defined $ref_allele && defined $alt_allele){
            warn "Skipping $l->[2] from $l->[7] as could not determine alleles\n";
            next;
        }
        my $len = length($ref_allele);
        $ref_allele = "$len\_base_deletion" if ($len > 4000);
        my $allele_string = "$ref_allele/$alt_allele";
        
        
        ## resolve start/end coordinates
        my $start = $l->[4];
        my $end   = $l->[5];
        if($change =~/ins/ && $change !~/del/ ){
            ($start, $end) = ($end, $start);
        }
        $start = $end + 1 if $change =~/dup/ ;

        ## an insertion should be between only 2 bases
        next if $end +1 < $start;
        

        ## resolve strand
        my $strand ;
        if ($l->[6] eq "+"){
            $strand = 1;
        }
        elsif ($l->[6] eq "-"){
            $strand = -1;
        }
        else{
            die "Unknown strand: $l->[6] \n";
        }

        my %var = ( "start"       =>  $start, 
                    "end"         =>  $end,
                    "strand"      =>  $strand, 
                    "seqreg_name" =>  $l->[3], 
                    "allele"      =>  $allele_string,
                    "label"       =>  $l->[2], 
                    "name"        =>  $l->[0],
                    "source"      =>  $l->[7],

        );
        ## check for duplicates
        $count{$l->[0]}++;
        $var{fail_reasons} = " ";
        unless  ($ref_allele =~/_base_deletion/){  ## can't do much with these
            $var{fail_reasons} = run_checks(\%var);
        }
        
        ## keep a count for fail rate
        $status{fail}++ if defined $var{fail_reasons} &&  $var{fail_reasons}=~/\d+/;        
        $status{all}++;
        
        ## print to tmp file 
        print $out "$var{name}\t$var{allele}\t$var{seqreg_name}\t$var{start}\t$var{end}\t$var{strand}\t$var{label}\t$var{source}\t$var{fail_reasons}\n";
    }

    close $out;

    ## print out fail rate
    my $per_fail = 100 * $status{fail} / $status{all};
    print "\nWARNING: $status{fail} variants($per_fail %) fail from $status{all})\n";

    return \%count;

}

## call standard QC checks & return string of failure reasons
sub run_checks{
    
    my $var = shift;
    
    my @fail;
    
    ## Type 3  flag variation as fail if it has [A/T/G/C] allele string 
        
    my $all_possible_check = check_four_bases($var->{allele});
    push @fail, 3 if ($all_possible_check ==1);
    
    ## Type 14 resolve ambiguities before reference check - flag variants & alleles as fails
    
    my $is_ambiguous = check_for_ambiguous_alleles( $var->{allele} );
    push @fail, 14  if(defined $is_ambiguous) ;
    
    
    # Extract reference sequence to run ref checks [ compliments for reverse strand multi-mappers]    
    
    my $ref_seq = get_reference_base($var, $db, "fasta_seq") ;
    
    unless(defined $ref_seq){ 
        ## don't check further if obvious coordinate error
        push @fail, 15;
        
        return ( join",", @fail );    
    }

    
    ## is ref base in agreement?
    my $exp_ref = (split/\//, $var->{allele} )[0] ; 
    push @fail, 2 unless (defined $exp_ref && "\U$exp_ref" eq "\U$ref_seq"); ## using soft masked seq
       
    
    ## is either allele of compatible length with given coordinates?
    my $ref = (split/\//, $var->{allele} )[0];
    my $match_coord_length = check_variant_size( $var->{start}, $var->{end}, $ref);
    push @fail, 15 unless  ($match_coord_length == 1);
    
    
    return (join",", @fail);
}


## look up or insert source
sub get_source {

    my $reg         = shift;
    my $source_data = shift;

    my $source_adaptor = $reg->get_adaptor('homo_sapiens', 'variation', 'source');

    my $source = $source_adaptor->fetch_by_name($source_data->{name});

    ### source already loaded
    return $source if defined $source;

    ## add new source
    $source = Bio::EnsEMBL::Variation::Source->new
       (-name        => $source_data->{name},
        -version     => $source_data->{version},
        -description => $source_data->{desc},
        -url         => $source_data->{url},
     );

    $source_adaptor->store($source);

    return $source if defined $source;

    die "Failed to find or enter source for $source_data->{name} \n";

}
 
## look up seq region ids for quick variation_feature creation
sub get_seq_ids{

    my $dbh = shift;

    my %seq_ids;

    my $seq_ext_sth = $dbh->prepare(qq[ select seq_region_id, name from seq_region]);
    $seq_ext_sth->execute();
    my $dat = $seq_ext_sth->fetchall_arrayref();

    foreach my $l(@{$dat}){
        $seq_ids{$l->[1]} = $l->[0];
    }

    return \%seq_ids;
}

sub enter_var{

    my $line    = shift;
    my $adaptor = shift;
    my $source  = shift;

    my $var = Bio::EnsEMBL::Variation::Variation->new_fast({
        name             => $line->[0],
        _source_id       => $source_ob->dbID,
        is_somatic       => 0
                                                           });
    $adaptor->store($var);

    if(defined $line->[8] && $line->[8] =~ /\d+/){
        ## add fail status
        my @fails;
        if( $line->[8] =~/\,/){  @fails = split /\,/, $line->[8];}
        else{ push @fails , $line->[8] ;} 
        my $var_id = $var->dbID();

        foreach my $type (@fails){
               $adaptor->dbc->do(qq[ insert into failed_variation (variation_id, failed_description_id) values ($var_id, $type ) ]);
         }
    }

    return $var;
}

sub enter_varfeat{      

    my $line    = shift;
    my $adaptor = shift;
    my $source  = shift;
    my $var     = shift;
    my $seq_ids = shift;

    my $varfeat = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
        variation_name   => $line->[0],
        _source_id       => $source_ob->dbID,
        allele_string    => $line->[1], 
        _variation_id    => $var->dbID(),                                          
        seq_region_id    => $seq_ids->{$line->[2]},
        start            => $line->[3],
        end              => $line->[4],
        strand           => $line->[5],
        is_somatic       => 0
                                                                      });

    $adaptor->store($varfeat);
}

sub enter_alleles{      

    my $line    = shift;
    my $adaptor = shift;
    my $var     = shift;

    my @alleles = split/\//, $line->[1];

    foreach my $allele (@alleles){

        my $allele_to_insert;
        my $len = length($allele);
        if($len > 100) { 
           $allele_to_insert = "$len\_base_deletion";
        }
        elsif(  $allele =~/deletion|insertion/){
           $allele_to_insert = $allele;
        }
        else{
           $allele_to_insert = "\U$allele";
        }
        my $al = Bio::EnsEMBL::Variation::Allele->new_fast({
            allele         =>$allele_to_insert ,
            variation      => $var
                                                           });
        $adaptor->store($al);

    }
}

## attrib ids
#| deletion       |        12 |
#| duplication    |       253 |
#| insertion      |        10 |

sub insert_svar{


    my $dbh     = shift;
    my $line    = shift;
    my $source  = shift;
    my $seq_ids = shift;

    my $svar_ins_sth = $dbh->prepare(qq[ insert into structural_variation 
                                         (variation_name, source_id, class_attrib_id)
                                         values ( ?,?,? )
                                       ]);  

    my $svar_ext_sth = $dbh->prepare(qq[ select structural_variation_id from structural_variation 
                                         where variation_name =? and source_id =?
                                       ]);

   my $svarf_ins_sth = $dbh->prepare(qq[ insert into structural_variation_feature 
                                         (seq_region_id, seq_region_start, seq_region_end, seq_region_strand,
                                         structural_variation_id, variation_name, source_id,  
                                         class_attrib_id, allele_string)
                                         values ( ?,?,?,?,?,?,?,?,?)
                                       ]);  


    $svar_ins_sth->execute($line->[0], $source, 12)||die;
    $svar_ext_sth->execute($line->[0], $source)||die;

    my $svar_id = $svar_ext_sth->fetchall_arrayref();
    die "Problem importing $line->[0] as sv\n" unless defined $svar_id->[0]->[0];

    $line->[1] = '\\N' if $line->[1] =~/deletion/;  ## no point entering descriptions like 9353_base_deletion as allele strings

    $svarf_ins_sth->execute($seq_ids->{$line->[2]}, $line->[3], $line->[4], $line->[5],$svar_id->[0]->[0], $line->[0], $source, 12, $line->[1])||die;

    return $svar_id->[0]->[0];

}

## add variation ids to set
sub add_to_set{

    my ($dbh, $variation_ids, $struct_variation_ids ) = @_;

    ## get PhenCode set
    my $set_id  = get_set($dbh);

    ## add variants
    my $vsv_ins_sth = $dbh->prepare(qq[ insert ignore into  variation_set_variation
                                       (variation_id, variation_set_id)
                                        values (?,?)] );

    foreach my $var ( keys %{$variation_ids} ){  
        $vsv_ins_sth->execute( $var, $set_id );
    }

    
    ## add structural variants
    my $vssv_ins_sth = $dbh->prepare(qq[ insert ignore into  variation_set_structural_variation
                                         (structural_variation_id, variation_set_id)
                                         values (?,?)] );

    foreach my $svar ( keys %{$struct_variation_ids} ){  
        $vssv_ins_sth->execute( $svar, $set_id );
    }


}

## look up or enter variation set
sub get_set{

    my $dbh = shift;
    
    my $set_ext_sth = $dbh->prepare(qq[ select variation_set_id from variation_set where name ='PhenCode']);

    my $set_ins_sth = $dbh->prepare(qq[insert into variation_set (  name, description, short_name_attrib_id) 
                                        values ( 'PhenCode', 
                                       'Variants from the PhenCode Project',
                                        355) ]);

    ### look for old set record
    $set_ext_sth->execute();
    my $id = $set_ext_sth->fetchall_arrayref();

    return $id->[0]->[0] if defined $id->[0]->[0] ;


    ### enter new set record
    $set_ins_sth->execute(); 
    $set_ext_sth->execute();
    $id = $set_ext_sth->fetchall_arrayref();

    return $id->[0]->[0] if defined $id->[0]->[0] ;


    ### give up
    die "ERROR: variation set could not be entered\n"; 

}

## extract list of variants already imported from dbSNP 
##    - no need to re-import
sub get_dbSNP_data{

    my $dbh = shift;

    my %already_loaded;

    my $phencode_ext_sth = $dbh->prepare(qq[ select vs.name
                                             from  variation_synonym vs, source
                                             where source.name = 'PhenCode'                                         
                                             and vs.source_id = source.source_id
                                           ]);

    $phencode_ext_sth->execute()||die;
    my $dat = $phencode_ext_sth->fetchall_arrayref();

    foreach my $l (@{$dat}){
        $already_loaded{$l->[0]} = 1;
    }

    return \%already_loaded;
}

sub usage{

    die "Usage:\n\timport_Phencode.pl  -fasta [genomic sequence file for QC] -registry [registry file]\n\n";
}
