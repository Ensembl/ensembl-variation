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


## script to read data export from Europe PMC & UCSC and import as variation_citations

use strict;
use warnings;

use HTTP::Tiny;
use XML::Simple;
use Getopt::Long;
use Data::Dumper;


use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Publication;
use Bio::EnsEMBL::Variation::DBSQL::PublicationAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::PhenotypeAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor;
use Bio::EnsEMBL::Variation::Utils::QCUtils qw(count_rows);

our $DEBUG = 0;

my ($registry_file, $data_file, $species, $check_dbSNP, $type, $do_disease, $no_evidence);

GetOptions ("data=s"        => \$data_file,
            "type=s"        => \$type,
            "species=s"     => \$species,
            "check_dbSNP:s" => \$check_dbSNP,
            "registry=s"    => \$registry_file,
            "no_evidence"   => \$no_evidence,
    );

## an export file is needed for EPMC data
usage() unless defined $registry_file && defined $type && defined  $species && (defined $data_file || $type eq "UCSC" || $type eq "disease") ;


my $http = HTTP::Tiny->new();


## only import papers to non-human databases if the species is mentioned
my %check_species = ("bos_taurus"      =>   [ "bos taurus",   "cow",   "bovine", "cattle" ],
                     "mus_musculus"    =>   [ "mus musculus", "mouse", "murine", "mice"],
                     "homo_sapiens"    =>   [ "homo sapiens", "human", "child", "patient" ],
                     "gallus_gallus"   =>   [ "gallus gallus", "chicken", "hen" ],
                     "canis_familiaris"=>   [ "canis familiaris", "dog", "canine"],
                     "rattus_norvegicus"=>  [ "rattus norvegicus", "rat"],
		     "sus_scrofa"      =>   [ "sus scrofa", "pig", "porcine" ],
		     "ovis_aries"      =>   [ "ovis aries", "sheep", "ovine"],
                     "equus_caballus"  =>   [ "equus caballus", "horse", "equine"],
    );

our $species_string = join "|", @{$check_species{$species}} if defined $check_species{$species};
print "Checking species $species_string\n";

my $reg = 'Bio::EnsEMBL::Registry';
$reg->no_version_check(1); 
$reg->load_all($registry_file);

my $dba = $reg->get_DBAdaptor($species, 'variation') || die "Error getting db adaptor\n";
## extract all variants - cited variants failing QC are still displayed
$dba->include_failed_variations(1);
our $pheno_adaptor = $reg->get_adaptor($species, 'variation', 'phenotype');


## if new dbSNP release has been imported, pull out full info for citations
my $dbSNP_data = check_dbSNP($dba) unless defined $check_dbSNP && $check_dbSNP ==0;

## check production db for suspected errors
my $avoid_list = get_avoid_list($reg);


## parse input data
my $file_data;

if( $type eq "EPMC"){

    ## read ePMC file, ommitting suspected errors
    $file_data = parse_EPMC_file($data_file, $avoid_list);

    ##import ny new publications & citations 
    import_citations($reg, $file_data);

    ## update variations & variation_features to have Cited status
    update_evidence($dba) unless defined $no_evidence;
}
elsif($type eq "UCSC"){

    ## download current data from UCSC public server unless already extracted
    $data_file = get_current_UCSC_data() unless defined $data_file;
    
    ## read UCSC file, ommitting suspected errors
    $file_data = parse_UCSC_file($data_file, $avoid_list);

    ##import any new publications & citations 
    import_citations($reg, $file_data);

    ## update variations & variation_features to have Cited status
    update_evidence($dba) unless defined $no_evidence;

}

elsif($type eq "disease"){
    ## look for phenotype citations in EPMC for all current publications 
    do_disease($reg);
}

else{
    die "Type $type is not recognised - must be EPMC or UCSC\n";
}


## create report on curent status
report_summary($dba, $species);

sub import_citations{

    my $reg            = shift;
    my $data           = shift;

    open my $error_log, ">>$species\.log"|| die "Failed to open log file: $!\n";


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
                no warnings ;
                ### write file of variants not found in this species to use as input file for next
                print $error_log "$rsid,$data->{$pub}->{pmcid},$data->{$pub}->{pmid},  No variant record\n";
                use warnings ;
            }
        }


        ## don't enter publication if there are no variants found for this species
        next unless defined $var_obs[0] ;
        
        
      
        ### Check if publication already known & enter if not  
        my $publication;

        ## looking up missing data from EPMC before redundancy check
        my $ref = get_publication_info_from_epmc($data, $pub, $error_log);
	
        my $title = (defined $ref->{resultList}->{result}->{title} ? $ref->{resultList}->{result}->{title}  : $data->{$pub}->{title});
        next unless defined $title;       
        next if $title =~/Erratum/i; 
 
        ## save ids 
        my $pmid   = $ref->{resultList}->{result}->{pmid}   || $data->{$pub}->{pmid};
        my $pmcid  = $ref->{resultList}->{result}->{pmcid}  || undef;
        my $doi    = $ref->{resultList}->{result}->{DOI}    || $data->{$pub}->{doi};

        ## try looking up on doi first
        $publication = $pub_ad->fetch_by_doi( $doi )     if defined $doi;
        ## then PMID
        $publication = $pub_ad->fetch_by_pmid( $pmid )   if defined $pmid && ! defined $publication;
        ## then PMCID   
        $publication = $pub_ad->fetch_by_pmcid( $pmcid ) if defined $pmcid && ! defined $publication;


        if(defined $publication){
            ##  warn "Linkings vars to existing pub ". $publication->dbID() ."\n";
            $pub_ad->update_variant_citation( $publication,\@var_obs );
            $pub_ad->update_ucsc_id( $publication,  $data->{$pub}->{ucsc} ) if defined $data->{$pub}->{ucsc};
        }
        else{
            ## add new publication

            ### create new object
            my $publication = Bio::EnsEMBL::Variation::Publication->new( 
                -title    => $title,
                -authors  => $ref->{resultList}->{result}->{authorString}   || $data->{$pub}->{authors},
                -pmid     => $ref->{resultList}->{result}->{pmid}           || $data->{$pub}->{pmid},
                -pmcid    => $ref->{resultList}->{result}->{pmcid}          || undef,
                -year     => $ref->{resultList}->{result}->{pubYear}        || $data->{$pub}->{year},
                -doi      => $ref->{resultList}->{result}->{DOI}            || $data->{$pub}->{doi},
                -ucsc_id  => $data->{$pub}->{ucsc}                          || undef,
                -variants => \@var_obs,
                -adaptor  => $pub_ad
                );
        
            $pub_ad->store( $publication);
        }
    }
}

sub get_publication_info_from_epmc{

    my $data      = shift;
    my $pub       = shift;
    my $error_log = shift;

    my $ref;

    ### check is species mentioned if not human?
    if ($species_string !~/human|homo/ && defined $data->{$pub}->{pmcid} ){         
        my $mined = get_epmc_data( "PMC/$data->{$pub}->{pmcid}/textMinedTerms/ORGANISM" );        
        my $looks_ok = check_species($mined ,$data) ;
        
        if ($looks_ok == 0 && $ref->{resultList}->{result}->{title} !~ /$species_string/){
            print $error_log "ALL\t$data->{$pub}->{pmid}\t$ref->{resultList}->{result}->{title} - species not mentioned\n";
            return undef;
        }
    }


    ### get data on publication from ePMC

    if( defined $data->{$pub}->{pmid} ){
        $ref   = get_epmc_data( "search/query=ext_id:$data->{$pub}->{pmid}%20src:med" );
    }
    elsif( defined $data->{$pub}->{doi} ){
        $ref   = get_epmc_data( "search/query=$data->{$pub}->{doi}" );
        ## check results of full text query
         return undef unless defined  $data->{$pub}->{doi} &&
	     $ref->{resultList}->{result}->{doi} eq $data->{$pub}->{doi}; 
    }
    elsif(defined $data->{$pub}->{pmcid}){
        $ref   = get_epmc_data( "search/query=$data->{$pub}->{pmcid}" );
        ## check results of full text query
         return undef unless defined $data->{$pub}->{pmcid} &&
	     $ref->{resultList}->{result}->{pmcid} eq $data->{$pub}->{pmcid};  
        
    }
    else{
        print $error_log "ALL\t$pub - nothing to search on\n";
         return undef;
    }
    
    ### check title available       
    unless (defined $ref->{resultList}->{result}->{title}){
        my $mess = "ALL\t";
        $mess .= "pmid:$data->{$pub}->{pmid}\t"   if defined $data->{$pub}->{pmid};
        $mess .= "pmcid:$data->{$pub}->{pcmid}\t" if defined $data->{$pub}->{pmcid};
        print $error_log $mess ." as no title\n";

         return undef;
    }        
   
    return $ref;
}
sub get_epmc_data{

    my $id = shift; ## specific part of URL including pmid or pmcid 

    return undef unless defined $id && $id =~/\d+/;

    my $xs   = XML::Simple->new();
    my $server = 'http://www.ebi.ac.uk/europepmc/webservices/rest/';
    my $request  = $server . $id;

    my %data;

    print "Looking for $request\n\n"  if $DEBUG == 1;
    my $response = $http->get($request, {
        headers => { 'Content-type' => 'application/xml' }
                              });
    unless ($response->{success}){
        warn "Failed request: $request :$!\n" ;
        return;
    }
    return $xs->XMLin($response->{content} );

}

## check if species is available in mined information
## essential for non-humans

sub check_species{

    my ($mined ,$data)=@_ ;
    my $looks_ok = 0;

    if(defined $mined->{semanticTypeList}->{semanticType}->{total}  &&
       $mined->{semanticTypeList}->{semanticType}->{total} ==1){

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

sub parse_EPMC_file{

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

    open my $error_log, ">>$species\_file.log"|| die "Failed to open log file:$!\n";

    my $pub_ext_sth = $dba->dbc()->prepare(qq[ select publication.publication_id, publication.pmid
                                               from publication
                                               where publication.title is null
                                               and publication.pmid is not null                     
                                              ]);

    my $pub_upd_sth = $dba->dbc()->prepare(qq[ update publication
                                               set publication.title = ?,
                                               publication.pmcid=?,
                                               publication.authors = ?,
                                               publication.year = ?,
                                               publication.doi = ?
                                               where publication.publication_id = ?
                                             ]);


    $pub_ext_sth->execute();
    my $dat = $pub_ext_sth->fetchall_arrayref();
    return unless defined $dat->[0]->[0];

    foreach my $l (@{$dat}){ 

        my $mined = get_epmc_data( "MED/$l->[1]/textMinedTerms/ORGANISM" );
        my $ref   = get_epmc_data("search/query=ext_id:$l->[1]%20src:med");
        unless(defined $mined && defined $ref){
            print $error_log "NO EPMC data for PMID:$l->[1] - skipping\n";
            next;
        }
        ## dbSNP may list publications which have not been mined
        my $looks_ok  = check_species($mined ,$ref) ;
        if ($looks_ok == 0 && defined $ref->{resultList}->{result}->{title} &&
	    $ref->{resultList}->{result}->{title} !~/$species_string/i){
            print $error_log "WARN dbSNP data:\t$l->[1]\t($ref->{resultList}->{result}->{title}) - species not mentioned\n"
        }

        $pub_upd_sth->execute( $ref->{resultList}->{result}->{title},
                               $ref->{resultList}->{result}->{pmcid},
                               $ref->{resultList}->{result}->{authorString},
                               $ref->{resultList}->{result}->{pubYear},
                               $ref->{resultList}->{result}->{doi},
                               $l->[0]
            ) if defined $ref->{resultList}->{result}->{title};
    }  
    close $error_log;  
}

sub update_evidence{

    my $dba = shift;

    ## find cited attrib
    my $attrib_ext_sth = $dba->dbc()->prepare(qq[ select attrib_id from attrib where value ='Cited']);
    $attrib_ext_sth->execute()||die;
    my $attrib =  $attrib_ext_sth->fetchall_arrayref();
    die "Not updating evidence as no attrib found\n" unless defined $attrib->[0]->[0];

    my $ev_ext_sth = $dba->dbc()->prepare(qq[ select variation.variation_id, variation.evidence_attribs 
                                              from variation, variation_citation
                                              where variation.variation_id = variation_citation.variation_id
                                             ]);

    my $var_upd_sth     = $dba->dbc()->prepare(qq[ update variation set evidence_attribs = ? where variation_id = ?]);
    my $varfeat_upd_sth = $dba->dbc()->prepare(qq[ update variation_feature set evidence_attribs = ? where variation_id = ?]);

    $ev_ext_sth->execute()||die;
    my $dat =  $ev_ext_sth->fetchall_arrayref();

    my $n = scalar @{$dat};
    warn "$n variant evidence statuses to update\n";

    foreach my $l (@{$dat}){

        my $evidence;
        if(defined $l->[1] && $l->[1] =~ /\w+/){
            $evidence .= "$l->[1],$attrib->[0]->[0]";
        }
        else{
            $evidence = "$attrib->[0]->[0]";
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

    my $dup1_ext_sth = $dba->dbc->prepare(qq[ select p1.publication_id, p2.publication_id, p2.pmid
                                         from publication p1, publication p2
                                         where p1.pmid = p2.pmid
                                         and p1.publication_id < p2.publication_id
                                         and p1.pmid is not null
                                       ]);

    my $dup2_ext_sth = $dba->dbc->prepare(qq[ select p1.publication_id, p2.publication_id, p2.pmcid
                                         from publication p1, publication p2
                                         where p1.pmcid = p2.pmcid
                                         and p1.publication_id < p2.publication_id
                                       ]);

   my $dup3_ext_sth = $dba->dbc->prepare(qq[ select p1.publication_id, p2.publication_id, p2.doi
                                         from publication p1, publication p2
                                         where p1.doi = p2.doi
                                         and p1.publication_id < p2.publication_id
                                       ]);




    my $fail_ext_sth = $dba->dbc->prepare(qq[ select count(*) from publication
                                              where title is null
                                             ]);

    $dup1_ext_sth->execute()||die;
    my $dup1 =  $dup1_ext_sth->fetchall_arrayref();

    $dup2_ext_sth->execute()||die;
    my $dup2 =  $dup2_ext_sth->fetchall_arrayref();

    $dup3_ext_sth->execute()||die;
    my $dup3 =  $dup3_ext_sth->fetchall_arrayref();

    $fail_ext_sth->execute() ||die;
    my $fail =  $fail_ext_sth->fetchall_arrayref();

    if (defined $dup1->[0]->[0] || defined $dup1->[0]->[0]){
        print $report "Duplicated publications:\n";

        foreach my $l (@{$dup1}){
            print $report "$l->[0]\t$l->[1]\t$l->[2]\n";
        }
        foreach my $k (@{$dup2}){
            print $report "$k->[0]\t$k->[1]\t$k->[2]\n";
        }   
        foreach my $m (@{$dup3}){
            print $report "$m->[0]\t$m->[1]\t$m->[2]\n";
        }

    }

    print $report "$fail->[0]->[0] publications without a title - to be deleted\n";

}


sub get_current_UCSC_data{

    my $filename = "UCSC_export";
    open my $out, ">$filename" || die "Failed to open file to write data : $!\n";

    my $dbh = DBI->connect('dbi:mysql:hgFixed:genome-mysql.cse.ucsc.edu:3306:max_allowed_packet=1MB', 'genome', '', undef);

    my $cit_ext_sth = $dbh->prepare(qq[ SELECT pma.markerId, pa.pmid, pma.section, pa.doi, pa.title, pa.authors, pa.year, pa.extId
                                        FROM pubsMarkerAnnot pma JOIN pubsArticle pa USING
                                        (articleId) WHERE pma.markerType="snp" 
                                      ]);


    $cit_ext_sth->execute()||die;

    while( my $line = $cit_ext_sth->fetchrow_arrayref()){

        next if $line->[6] < 1999;             ## pre-dbSNP - must be random match
        next if $line->[7] eq "PMC$line->[1]"; ## incorrect/ missing PMID 
        next if $line->[4] eq "NotFound" || $line->[4] eq "TOC" || $line->[4] eq "Highlights"
            || $line->[4]  eq "Contents" || $line->[4] eq "Table of Contents" || $line->[4]  eq "Volume Contents"
            || $line->[4]  eq "Cannabis" || $line->[4] eq "Index"  || $line->[4] eq "Author Index"
            || $line->[4]  =~/This Issue/i      || $line->[4]  =~/This Month in/
            || $line->[4] eq "Ensembl 2011"; 

        print $out join("\t", @{$line}) . "\n";
    }
    close $out;
    return $filename;
}

sub parse_UCSC_file{

   my $ucsc_file  = shift;
   my $avoid_list = shift;

    my %data;
    open my $file, $ucsc_file || die "Failed to open file of UCSC citations $ucsc_file : $!\n";
    while (<$file>){

        next unless /rs\d+/i;
        chomp;
        s/RS/rs/;
        s/^\|\s+//;

        #my ($rs, $pmid, $section, $doi, $title, $authors, $year ) = split/\s+\|\s+|\t/;
        my ($rs, $pmid, $section, $doi, $title, $authors, $year, $extId ) = split/\t/;
        next if $section =~/refs|ack/; ## not in this publication

        $pmid = "" if $pmid eq "0";  ## nulls are set to 0 in UCSC database

        ## remove known errors
         if ( $avoid_list->{$rs}->{$pmid}){
            warn "Removing suspected error : $rs in $pmid\n";  
            next;
        }

        ## Group by publication; may have pmid, doi or both
        my $tag = $pmid . "_" .$doi;


        $data{$tag}{pmid} = $pmid;
        $data{$tag}{pmid} = undef unless $pmid =~ /\d+/;

        $data{$tag}{doi}     = $doi;
        $data{$tag}{title}   = $title;
        $data{$tag}{authors} = $authors;
        $data{$tag}{year}    = $year;
        $data{$tag}{ucsc}    = $extId;

        push @{$data{$tag}{rsid}}, $rs;
    }

    return \%data;
}
 

sub do_disease{

    my $reg = shift;

    my $dba = $reg->get_DBAdaptor($species, 'variation') || die "Error getting db adaptor\n";
   # my $pheno_adaptor = $reg->get_adaptor($species, 'variation', 'phenotype');


    my $phencit_ins_sth = $dba->dbc()->prepare(qq[ insert into phenotype_citation ( publication_id, phenotype_id ) values ( ?, ?) ]);
    my $phencit_ext_sth = $dba->dbc()->prepare(qq[ select publication_id from phenotype_citation where publication_id = ? and  phenotype_id = ? ]);


    my $pheno_ext_sth = $dba->dbc()->prepare(qq[ select phenotype_id, description from phenotype]);


    my $pub_ext_sth = $dba->dbc()->prepare(qq[ select pmid, publication_id
                                               from publication 
                                               where publication_id not in 
                                              (select publication_id from phenotype_citation)
                                               and pmid is not null
                                             ]);

    my %pheno_ids;
    $pheno_ext_sth->execute()||die;
    my $curr_pheno = $pheno_ext_sth->fetchall_arrayref();
    foreach my $l (@{$curr_pheno}){
        $pheno_ids{"\U$l->[1]"} = $l->[0];
    }


    $pub_ext_sth->execute()||die;
    my $dat =  $pub_ext_sth->fetchall_arrayref();

    foreach my $l (@{$dat}){
        #warn "looking for disease for $l->[0]\n";
        ## extract disease info

=head  This bit does one disease only
        my $top_disease = find_disease($l->[0]);
        ## skip paper if none found      
        next unless defined $top_disease

=cut
        my $multi_disease = find_disease($l->[0]);
        ## skip paper if none found     
        next unless defined $multi_disease->[0];
        foreach my $top_disease(@{$multi_disease}){
            $top_disease =~ s/\'//g;
            $top_disease = "\U$top_disease";
            #warn "looking for/adding $top_disease\n";
            
            if(defined  $pheno_ids{$top_disease}){
                
                $phencit_ext_sth->execute( $l->[1], $pheno_ids{$top_disease} )||die;
                my $already_done = $phencit_ext_sth->fetchall_arrayref();
                
                unless(defined $already_done ->[0]->[0]){
                    ## create phenotype citation
                    $phencit_ins_sth->execute($l->[1], $pheno_ids{$top_disease} )||die;
                }
            }
            else{
            ## enter if new phenotype - warn for quick rubbish scanning
                #warn "Entering new phenotype  $top_disease\n";
                my $pheno = Bio::EnsEMBL::Variation::Phenotype->new(-DESCRIPTION => $top_disease);
                $pheno_adaptor->store($pheno );
                $pheno_ids{$top_disease} =  $pheno->dbID() ;
                ## create phenotype citation
                $phencit_ins_sth->execute($l->[1], $pheno_ids{$top_disease})||die;
            }               
        }
    }        
}

sub find_disease{

    my $pmid = shift;

    my %count;
    my @diesase;

    my $mined = get_epmc_data( "MED/$pmid/textMinedTerms/DISEASE" );
    if( ref($mined->{semanticTypeList}->{semanticType}->{tmSummary}) eq 'ARRAY' ){
        foreach my $found (@{$mined->{semanticTypeList}->{semanticType}->{tmSummary}}  ){
            next if $found->{count} < 2;
            push @diesase, $found->{term};

        }
        
        
        
    }
    else{
        ## if there is only one disease found, return it unless it is only mentioned once
        return undef if $mined->{semanticTypeList}->{semanticType}->{tmSummary}->{count} < 2;
        push @diesase, $mined->{semanticTypeList}->{semanticType}->{tmSummary}->{term};

    }
    return \@diesase;
    
}


sub find_disease_max{

    my $pmid = shift;

    my %count;
    my @diesase;

    my $mined = get_epmc_data( "MED/$pmid/textMinedTerms/DISEASE" );
    if( ref($mined->{semanticTypeList}->{semanticType}->{tmSummary}) eq 'ARRAY' ){
        foreach my $found (@{$mined->{semanticTypeList}->{semanticType}->{tmSummary}}  ){
            next if $found->{count} < 2;
            $count{ $found->{term}} = $found->{count} ;
        }
        
        my $max = 0;
        ## find highest number of disease mentions in publication
        foreach my $disease (keys %count){
            $max = $count{$disease} if $max < $count{$disease};     
        }

        ## store diseases reaching this threshhold
        foreach my $disease (keys %count){
            push @diesase, $disease if $max == $count{$disease};
            
        }
        return undef if scalar(@diesase) >1 ; ## give up - not specific enough
        return $diesase[0];
        
    }
    else{
        ## if there is only one disease found, return it unless it is only mentioned once
        return undef if $mined->{semanticTypeList}->{semanticType}->{tmSummary}->{count} < 2;
        return $mined->{semanticTypeList}->{semanticType}->{tmSummary}->{term};

    }

    
}

sub usage {
    
    die "\n\tUsage: import_EPMC.pl -type [ EPMC or UCSC or disease] -species [name] -registry [registry file]

Options:  
          -data [file of citations]   - *required* for EPMC import
          -check_dbSNP [0/1]          - add detail to citations from dbSNP import (default:1)
          -no_evidence                - don't update variation & variation_feature evidence statuses\n\n";

}

sub unique {
  my %a;
  map { $a{$_} = 1; } @_;
  return sort keys %a;
}

