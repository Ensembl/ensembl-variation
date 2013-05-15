#!/usr/bin/env perl

## script to read data export from Europe PMC and import as variation_citations

BEGIN{
    $ENV{http_proxy} = 'http://wwwcache.sanger.ac.uk:3128'; 

}
use strict;
use warnings;

use HTTP::Tiny;
use XML::Simple;
use Getopt::Long;
use Data::Dumper;


use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Publication;
use Bio::EnsEMBL::Variation::DBSQL::PublicationAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor;
use Bio::EnsEMBL::Variation::Utils::QCUtils qw(count_rows);

our $DEBUG = 0;

my ($registry_file, $data_file, $species, $check_dbSNP);

GetOptions ("data=s"        => \$data_file,
            "species=s"     => \$species,
            "check_dbSNP:s" => \$check_dbSNP,
            "registry=s"    => \$registry_file,
    );

usage() unless defined $registry_file && defined $data_file && defined  $species;


my $xs   = XML::Simple->new();
my $http = HTTP::Tiny->new();


## only import papers to non-human databases if the species is mentioned
my %check_species = ("bos_taurus"      =>   [ "bos taurus",   "cow",   "bovine", "cattle" ],
                     "mus_musculus"    =>   [ "mus musculus", "mouse", "murine", "mice"],
                     "homo_sapiens"    =>   [ "homo sapiens", "human", "child", "patient" ],
                     "gallus_gallus"   =>   [ "gallus gallus", "chicken", ],
                     "canis_familiaris"=>   [ "canis familiaris", "dog", "canine"],
                     "rattus_norvegicus"=>  [ "rattus norvegicus", "rat"],
    );

our $species_string = join "|", @{$check_species{$species}} if defined $check_species{$species};
print "Checking species $species_string\n";

my $reg = 'Bio::EnsEMBL::Registry';
$reg->no_version_check(1); 
$reg->load_all($registry_file);
my $dba = $reg->get_DBAdaptor($species, 'variation') || die "Error getting db adaptor\n";



## if new dbSNP release has been imported, pull out full info for citations
my $dbSNP_data = check_dbSNP($dba) unless defined $check_dbSNP && $check_dbSNP ==0;

## check production db for suspected errors
my $avoid_list = get_avoid_list($reg);

## read ePMC file, ommitting suspected errors
my $file_data = parse_file($data_file, $avoid_list);

##import ePMC info
import_citations($reg, $file_data, $species_string);


## update variations & variation_features to have Cited status
update_evidence($dba);

## create report on curent status
report_summary($dba, $species);

sub import_citations{

    my $reg            = shift;
    my $data           = shift;
    my $species_string = shift;

    open my $error_log, ">>$data_file\_$species\.log"|| die "Failed to open log file: $!\n";


    my $var_ad = $reg->get_adaptor($species, 'variation', 'variation');
    my $pub_ad = $reg->get_adaptor($species, 'variation', 'publication');

    
    foreach my $pub(keys %$data){
        
        ### get a set of variation objects
        my @var_obs;
        ### remove duplicate ids
        my @var_id = unique(@{$data->{$pub}->{rsid}});
        
        foreach my $rsid ( @var_id ){  
            my $v = $var_ad->fetch_by_name($rsid);
            if (defined $v){
                push @var_obs, $v;
            }
            else{
                ### write file of variants not found in this species to use as input file for next
                print $error_log "$rsid,$data->{$pub}->{pmcid},$data->{$pub}->{pmid},  No variant record\n";
            }
        }


        ## don't enter publication if there are no variants found for this species
        unless (defined $var_obs[0] ){
            next;
        }
        
        ### get data on publication from ePMC
        my ($ref, $mined);
        if( defined $data->{$pub}->{pmid} ){
            $ref   = get_epmc_data( "search/query=ext_id:$data->{$pub}->{pmid}%20src:med" );
            $mined = get_epmc_data( "MED/$data->{$pub}->{pmid}/textMinedTerms/ORGANISM" );
        }
        else{
            $ref   = get_epmc_data( "search/query=$data->{$pub}->{pmcid}" );
            $mined = get_epmc_data( "PMC/$data->{$pub}->{pmcid}/textMinedTerms/ORGANISM" );
        }
        

        ### check title available       
        unless (defined $ref->{resultList}->{result}->{title}){
            print $error_log "ALL\t$data->{$pub}->{pmid}\t$data->{$pub}->{pmcid} - as no title\n";
            next ;
        }


        ### check is species mentioned?
        my $looks_ok = check_species($mined ,$data) ;
        
        if ($looks_ok == 0 && $ref->{resultList}->{result}->{title} !~ /$species_string/){
            print $error_log "ALL\t$data->{$pub}->{pmid}\t$ref->{resultList}->{result}->{title} - species not mentioned\n";
            ## can't really check human
            next unless $species_string =~/human/;
        }

      
        ### Everything passes - import publication (if new) & links to variants

        my $publication;
        if(defined $data->{$pub}->{pmid} ){
            $publication = $pub_ad->fetch_by_pmid( $data->{$pub}->{pmid} );
        }
        else{
            $publication = $pub_ad->fetch_by_pmcid( $data->{$pub}->{pmcid} );
        }
        



        if(defined $publication){
#       warn "Linkings vars to existing pub\n";
            $pub_ad->update_variant_citation( $publication,\@var_obs, );
        }
        else{
            ### create new object
            my $publication = Bio::EnsEMBL::Variation::Publication->new( -title    => $ref->{resultList}->{result}->{title},
                                                                         -authors  => $ref->{resultList}->{result}->{authorString},
                                                                         -pmid     => $ref->{resultList}->{result}->{pmid},
                                                                         -pmcid    => $ref->{resultList}->{result}->{pmcid},
                                                                         -variants => \@var_obs,
                                                                         -adaptor  => $pub_ad
                );
        
            $pub_ad->store( $publication);
        }
    }
}


sub get_epmc_data{

    my $id = shift; ## specific part of URL including pmid or pmcid 

    my $server = 'http://www.ebi.ac.uk/europepmc/webservices/rest/';
    my $request  = $server . $id;

    my %data;

    print "Looking for $request\n\n"  if $DEBUG == 1;
    my $response = $http->get($request, {
        headers => { 'Content-type' => 'application/xml' }
                              });
    warn "Failed! :$!\n" unless $response->{success};
   
    return $xs->XMLin($response->{content} );

}

## check if species is available in mined information
## essential for non-humans

sub check_species{

    my ($mined ,$data)=@_ ;
    my $looks_ok = 0;

    if($mined->{semanticTypeList}->{semanticType}->{total} ==1){
        print "found spec ". $mined->{semanticTypeList}->{semanticType}->{tmSummary}->{term} . "\n"  if $DEBUG == 1;
        if ($mined->{semanticTypeList}->{semanticType}->{tmSummary}->{term} =~ /$species_string/i){
            $looks_ok = 1;
        }
    }
    else{
        foreach my $found (@{$mined->{semanticTypeList}->{semanticType}->{tmSummary}}  ){
            print "found spec ". $found->{term} . "\n"  if $DEBUG == 1;
            if ($found->{term} =~ /$species_string/i){
                $looks_ok = 1;
            }
        }
    }     
    return $looks_ok;
}

## read input
##
## expected format:  rsID,PMCID,PMID
##                   rs10016388,PMC1513515,16770606
##                   RS10028945,PMC1513515,16770606

sub parse_file{

    my $pubmed_file = shift;
    my $avoid_list  = shift;

    my %data;
    open my $file, $pubmed_file || die "Failed to open file of PubMed ids $pubmed_file : $!\n";
    while (<$file>){
        ## basic check that format is as expected
        next unless /rs\d+\,/i;
        chomp;
        s/RS/rs/;

        my ($rs, $pmcid, $pmid) = split/\,/;

        ## remove known errors
         if ($avoid_list->{$rs}->{$pmcid} || $avoid_list->{$rs}->{$pmid}){
            warn "Removing suspected error : $rs in $pmcid\n";  
            next;
        }
        ## Group by publication; may have pmid, pmcid or both
        my $tag = $pmcid . "_" .$pmid;


        $data{$tag}{pmid} = $pmid;
        $data{$tag}{pmid} = undef unless $pmid =~ /\d+/;

        $data{$tag}{pmcid} = $pmcid;
        push @{$data{$tag}{rsid}}, $rs;
    }

    return \%data;
}

## look up manually curated list of suspect citations
sub get_avoid_list{

    my $reg = shift;

    my %avoid;
    my $intdba = $reg->get_DBAdaptor('multi', 'intvar') || warn "Error getting db adaptor\n";

    unless (defined $intdba){
        warn "No avoid list will be used as db adaptor for internal db not found\n";
        return \%avoid;
    }
    my $dat_ext_sth = $intdba->dbc->prepare(qq[ select pmid, pmc_id, variation_name from failed_citations]);
    $dat_ext_sth->execute();
    my $dat = $dat_ext_sth->fetchall_arrayref();

    foreach my $l (@{$dat}){
        $avoid{$l->[2]}{$l->[0]} = 1;
        $avoid{$l->[2]}{$l->[1]} = 1;
    }
    return \%avoid;

}

### check for citations from raw dbSNP import & add detail
sub check_dbSNP{

   
    my $dba       = shift;


    open my $error_log, ">>$data_file\_$species\_file.log"|| die "Failed to open log file:$!\n";

    my $pub_ext_sth = $dba->dbc()->prepare(qq[ select publication.publication_id, publication.pmid
                                               from publication
                                               where publication.pmcid not like 'PMC%'                      
                                              ]);

    my $pub_upd_sth = $dba->dbc()->prepare(qq[ update publication
                                               set publication.title = ?,
                                               publication.pmcid=?,
                                               publication.authors = ?
                                               where publication.publication_id = ?
                                             ]);


    $pub_ext_sth->execute();
    my $dat = $pub_ext_sth->fetchall_arrayref();
    return unless defined $dat->[0]->[0];

    foreach my $l (@{$dat}){ 

        my $mined = get_epmc_data( "MED/$l->[1]/textMinedTerms/ORGANISM" );
        my $ref   = get_epmc_data("search/query=ext_id:$l->[1]%20src:med");

        ## dbSNP may list publications which have not been mined
        my $looks_ok  = check_species($mined ,$ref) ;
        if ($looks_ok == 0 && $ref->{resultList}->{result}->{title} !~/$species_string/i){
            print $error_log "WARN dbSNP data:\t$l->[1]\t($ref->{resultList}->{result}->{title}) - species not mentioned\n"
        }

        $pub_upd_sth->execute( $ref->{resultList}->{result}->{title},
                               $ref->{resultList}->{result}->{pmcid},
                               $ref->{resultList}->{result}->{authorString},
                               $l->[0]
            ) if defined $ref->{resultList}->{result}->{title};
    }  
    close $error_log;  
}

sub update_evidence{

    my $dba = shift;

    my $ev_ext_sth = $dba->dbc()->prepare(qq[ select variation.variation_id, variation.evidence 
                                              from variation, variation_citation
                                              where variation.variation_id = variation_citation.variation_id
                                              and variation.evidence not like '%Cited%'
                                             ]);

    my $var_upd_sth     = $dba->dbc()->prepare(qq[ update variation set evidence = ? where variation_id = ?]);
    my $varfeat_upd_sth = $dba->dbc()->prepare(qq[ update variation_feature set evidence = ? where variation_id = ?]);

    $ev_ext_sth->execute()||die;
    my $dat =  $ev_ext_sth->fetchall_arrayref();

    my $n = scalar @{$dat};
    warn "$n variant evidence statuses to update\n";

    foreach my $l (@{$dat}){

        my $evidence;
        if($l->[1] =~ /\w+/){
            $evidence .= "$l->[1],Cited";
        }
        else{
            $evidence = "Cited";
        }
        $var_upd_sth->execute($evidence, $l->[0]);
        $varfeat_upd_sth->execute($evidence, $l->[0]);
    }
}

sub report_summary{

    my $dba     = shift;
    my $species = shift;

    open my $report, ">import_EPMC_$species\.txt" ||die "Failed to open report file to write: $!\n";

    my $publication_count = count_rows($dba, 'publication');
    my $citation_count    = count_rows($dba, 'variation_citation');

    print $report "Total publications:\t$publication_count\n";
    print $report "Total citations:\t$citation_count\n";

    my $dup1_ext_sth = $dba->prepare(qq[ select p1.publication_id, p2.publication_id, p2.pmid
                                         from publication p1, publication p2
                                         where p1.pmid = p2.pmid
                                         and p1.publication_id < p2.publication_id
                                         and p1.pmid is not null
                                       ]);

    my $dup2_ext_sth = $dba->prepare(qq[ select p1.publication_id, p2.publication_id, p2.pmcid
                                         from publication p1, publication p2
                                         where p1.pmcid = p2.pmcid
                                         and p1.publication_id < p2.publication_id
                                         and p1.pmid is null
                                       ]);

    $dup1_ext_sth->execute()||die;
    my $dup1 =  $dup1_ext_sth->fetchall_arrayref();

    $dup2_ext_sth->execute()||die;
    my $dup2 =  $dup1_ext_sth->fetchall_arrayref();


    if (defined $dup1->[0]->[0] || defined $dup1->[0]->[0]){
        print $report "Duplicated publications:\n";

        foreach my $l (@{$dup1}){
            print $report "$l->[0]\t$l->[1]\t$l->[2]\n";
        }
        foreach my $k (@{$dup2}){
            print $report "$k->[0]\t$k->[1]\t$k->[2]\n";
        }   
    }
}


sub usage {

    die "\n\tUsage: import_EPMC.pl -data [EPMC file] -species [name] -registry [registry file]

Options:  -check_dbSNP [0/1]   - add detail to citations from dbSNP import (default:1)\n\n";

}

sub unique {
  my %a;
  map { $a{$_} = 1; } @_;
  return sort keys %a;
}

