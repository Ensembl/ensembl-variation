#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2023] EMBL-European Bioinformatics Institute
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
## the script also imports citations from phenotype tables (phenotype_feature and phenotype_feature_attrib) only for human

use strict;
use warnings;

use HTTP::Tiny;
use XML::Simple;
use Getopt::Long;
use utf8;
use Text::Unidecode;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Publication;
use Bio::EnsEMBL::Variation::DBSQL::PublicationAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::PhenotypeAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor;
use Bio::EnsEMBL::Variation::Utils::QCUtils qw(count_rows);
use Bio::EnsEMBL::Variation::Utils::Date qw( run_date log_time);
use Bio::EnsEMBL::Variation::Utils::SpecialChar qw(replace_hex);

our $DEBUG = 0;

my ($registry_file, $data_file, $species, $check_dbSNP, $clean, $type, $no_evidence);

GetOptions ("data=s"        => \$data_file,
            "type=s"        => \$type,
            "species=s"     => \$species,
            "check_dbSNP:s" => \$check_dbSNP,
            "clean:s"       => \$clean,
            "registry=s"    => \$registry_file,
            "no_evidence"   => \$no_evidence,
    );

## an export file is needed for EPMC data
usage() unless defined $registry_file && defined $type && defined  $species && (defined $data_file || $type eq "UCSC" || $type eq "phenotype");

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
                     "capra_hircus"    =>   [ "capra hircus", "goat", "caprine"],
    );

our $species_string = join "|", @{$check_species{$species}} if defined $check_species{$species};
print "Checking species $species_string\n";

# phenotype option is only available for human
# print error message
if($type eq 'phenotype' && $species_string !~/human|homo/) {
  die "Type $type is only available for human\n";
}

my $reg = 'Bio::EnsEMBL::Registry';
$reg->no_version_check(1); 
$reg->load_all($registry_file);

my $dba = $reg->get_DBAdaptor($species, 'variation') or die "Error getting db adaptor\n";
## extract all variants - cited variants failing QC are still displayed
$dba->include_failed_variations(1);

## if new dbSNP release has been imported, pull out full info for citations
my $dbSNP_data = check_dbSNP($dba) unless defined $check_dbSNP && $check_dbSNP ==0;

## check production db for suspected errors
my $avoid_list = get_avoid_list($reg);

## parse input data
my $file_data;

if( $type eq "EPMC"){

    ## read ePMC file, ommitting suspected errors
    $file_data = parse_EPMC_file($data_file, $avoid_list);

    ##import any new publications & citations
    import_citations($reg, $file_data, $type);

}
elsif($type eq "UCSC"){

    ## download current data from UCSC public server unless already extracted
    $data_file = get_current_UCSC_data() unless defined $data_file;
    
    ## read UCSC file, ommitting suspected errors
    $file_data = parse_UCSC_file($data_file, $avoid_list);

    ##import any new publications & citations
    import_citations($reg, $file_data, $type);

}
elsif($type eq "phenotype"){

  print "Publications from Phenotype feature table:\n";

  my $var_ad = $reg->get_adaptor($species, 'variation', 'variation');
  my $pub_ad = $reg->get_adaptor($species, 'variation', 'publication');
  my $source_ad = $reg->get_adaptor($species, 'variation', 'source');

  # Fetch all citation source attribs
  my $citation_attribs = get_citation_attribs($dba);

  my $citations_pheno_feature = process_phenotype_feature($dba, $source_ad, $citation_attribs);
  import_citations($reg, $citations_pheno_feature, $type);

  my $citations_pheno_feature_attrib = process_phenotype_feature_attrib($dba, $source_ad, $citation_attribs);
  import_citations($reg, $citations_pheno_feature_attrib, $type);

  print "Removing outdated citations...\n";
  remove_outdated_citations($dba, $var_ad, $pub_ad, $citation_attribs, $citations_pheno_feature, $citations_pheno_feature_attrib);
  print "Removing outdated citations... done!\n";
}
else{
    die "Type $type is not recognised - must be EPMC, UCSC or phenotype\n";
}

# create report on curent status
my $title_null = report_summary($dba, $species);

# add run date to meta table
update_meta($dba, $type) unless $type eq "phenotype";

# clean publications after import
if(defined $clean && $clean == 1){
  clean_publications($dba, $title_null);
}

# update variations & variation_features to have Cited status
update_evidence($dba) unless defined $no_evidence;

sub import_citations{

    my $reg            = shift;
    my $data           = shift;
    my $type           = shift;

    open (my $error_log, ">>$species\_$type\_errors\_" . log_time() . "\.log") or die "Failed to open log file: $!\n";
    ## store variants not from this species
    open (my $not_found, ">>rsID\_$type\_notfound\_$species\_" . log_time() . ".csv") or die "Failed to open file: $!\n";

    if($type eq "EPMC"){ 
      print $not_found "refsnp,PMCID,EXTID,SOURCE\n";
    }
    elsif($type eq "UCSC"){
      print $not_found "rsID\tPMID\tDOI\tTitle\tAuthors\tYear\tUCSC\n";
    }
    elsif($type eq "phenotype"){
      print $not_found "rsID,PMID\n";
    }

    my $var_ad = $reg->get_adaptor($species, 'variation', 'variation');
    my $pub_ad = $reg->get_adaptor($species, 'variation', 'publication');

    # get list of citations already in the db
    my $dba = $reg->get_DBAdaptor($species, 'variation') or die "Error getting db adaptor\n";
    my $done_list = get_current_citations($dba);

    # Get attrib id for source EPMC and UCSC
    my $source_attrib_id = get_source_attrib_id($reg, $type) unless ($type eq 'phenotype');

    foreach my $pub(keys %$data){

        ### get a set of variation objects
        my @var_obs;
        ### remove duplicate ids
        my @var_id = unique(@{$data->{$pub}->{rsid}});

        foreach my $rsid ( @var_id ){
            my $pub_pmid = $data->{$pub}->{pmid};
            $source_attrib_id = $data->{$pub}->{source} if($type eq 'phenotype');

            my $v = $var_ad->fetch_by_name($rsid);
            if (defined $v){
                push @var_obs, $v;

                # Update column data_source_attrib in table variation_citation
                if(defined $pub_pmid && $done_list->{$rsid}{$pub_pmid}){
                  my $source_data = $done_list->{$rsid}{$pub_pmid};
                  my $var_id = $v->dbID();
                  $pub_ad->update_citation_data_source($source_attrib_id, $var_id, $pub_pmid) unless $source_data =~ $source_attrib_id;
                }
            }
            else{
                no warnings ;
                ### write file of variants not found in this species to use as input file for next
                if($type eq "EPMC"){
                  print $not_found "\"$rsid\",$data->{$pub}->{pmcid},$pub_pmid,MED\n";
                }
                elsif($type eq "UCSC") {
                  print $not_found $rsid ."\t". $pub_pmid."\t-\t". $data->{$pub}->{doi} ."\t". $data->{$pub}->{title}."\t". $data->{$pub}->{authors} ."\t". $data->{$pub}->{year} ."\t".  $data->{$pub}->{ucsc} . "\n";
                }
                elsif($type eq "phenotype") {
                  print $not_found "$rsid,$pub_pmid\n";
                }
               use warnings ;
            }
        }

        # don't enter publication if there are no variants found for this species
        next unless defined $var_obs[0];
        
        ### Check if publication already known & enter if not
        my $publication;

        ## looking up missing data from EPMC before redundancy check
        my $ref = get_publication_info_from_epmc($data, $pub, $error_log);
        
        my $title = (defined $ref->{resultList}->{result}->{title} ? $ref->{resultList}->{result}->{title}  : $data->{$pub}->{title});
        next unless defined $title;

        # Publications with invalid titles are skipped
        # UCSC imports publications where title is not valid, example: 'P1-200'
        next if ($title =~/Erratum/i || $title =~/Errata/i || $title =~/Not Available/i || $title =~/^P[0-9]+\-[0-9]+$/i);

        # Delete brackets from title, e.g. [title].
        if($title =~ /^\[ ?[A-Za-z]{2}/){
          $title =~ s/^\[//;
          $title =~ s/\]\.//;
        }  

        # Some publication have newline in the title
        # Clean the title before it's inserted
        # example: PMID = 33498513
        if($title =~ /\n/) {
          $title =~ s/\n//;
          $title =~ s/\s\s+/ /;
        }
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
            $pub_ad->update_variant_citation( $publication,$source_attrib_id,\@var_obs );
            $pub_ad->update_ucsc_id( $publication,  $data->{$pub}->{ucsc} ) if defined $data->{$pub}->{ucsc};
        }
        else{
            ## add new publication
            ## check title size
          my $new_title = $title;
          if(length($title) >= 300){
            # truncate title
            my $aux = substr($title, 0, 296);
            my @list_title = split(' ', $aux);
            pop @list_title;
            $new_title = join(' ', @list_title);
            $new_title .= '...';
          }
          
          
          $new_title =~ s|<.+?>||g;
            ### create new object
          my $publication = Bio::EnsEMBL::Variation::Publication->new(
            -title    => $new_title,
            -authors  => unidecode($ref->{resultList}->{result}->{authorString})   || unidecode($data->{$pub}->{authors}),
            -pmid     => $ref->{resultList}->{result}->{pmid}           || $data->{$pub}->{pmid},
            -pmcid    => $ref->{resultList}->{result}->{pmcid}          || $data->{$pub}->{pmcid},
            -year     => $ref->{resultList}->{result}->{pubYear}        || $data->{$pub}->{year},
            -doi      => $ref->{resultList}->{result}->{DOI}            || $data->{$pub}->{doi},
            -ucsc_id  => $data->{$pub}->{ucsc}                          || undef,
            -variants => \@var_obs,
            -adaptor  => $pub_ad
          );
        
          $pub_ad->store( $publication,$source_attrib_id );
        }
    }
    close $not_found;
}

sub get_source_attrib_id{
  my $reg = shift;
  my $source = shift;

  my $attrib_adaptor = $reg->get_adaptor($species, 'variation', 'attribute');

  my $attrib_id = $attrib_adaptor->attrib_id_for_type_value('citation_source',$source);

  die "No attribute '$source' found\n" unless defined $attrib_id;

  return $attrib_id;
}

sub get_publication_info_from_epmc{

    my $data      = shift;
    my $pub       = shift;
    my $error_log = shift;

    my $ref;

    ### check is species mentioned if not human?
    if ($species_string !~/human|homo/ && defined $data->{$pub}->{pmcid} ){         
        my $pmcid = $data->{$pub}->{pmcid};
        $pmcid  =~ s/PMC//;
        my $mined = get_epmc_data( "annotations_api/annotationsByArticleIds?articleIds=PMC:$pmcid&type=Organisms&format=XML" );        
        my $looks_ok = check_species($mined ,$data) ;
        
        if ($looks_ok == 0 && $ref->{resultList}->{result}->{title} !~ /$species_string/){
            print $error_log "ALL\t$data->{$pub}->{pmid}\t$ref->{resultList}->{result}->{title} - species not mentioned\n";
            return undef;
        }
    }


    ### get data on publication from ePMC
    if( defined $data->{$pub}->{pmid} ){
      $ref = get_epmc_data( "webservices/rest/search?query=ext_id:$data->{$pub}->{pmid}%20src:med" );
    }
    elsif( defined $data->{$pub}->{doi} ){
      $ref = get_epmc_data( "webservices/rest/search?query=$data->{$pub}->{doi}" );
      ## check results of full text query
      return undef unless defined $data->{$pub}->{doi} && defined $ref->{resultList}->{result}->{doi} &&
      $ref->{resultList}->{result}->{doi} eq $data->{$pub}->{doi};
    }
    elsif(defined $data->{$pub}->{pmcid}){
      $ref = get_epmc_data( "webservices/rest/search?query=$data->{$pub}->{pmcid}" );
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
        $mess .= "pmid:$data->{$pub}->{pmid}\t" if defined $data->{$pub}->{pmid};
        $mess .= "pmcid:$data->{$pub}->{pcmid}\t" if defined $data->{$pub}->{pmcid};
        print $error_log $mess ." as no title\n";

         return undef;
    }

    ### format authors list
    my $trimmed_author_list = trim_author_list($ref->{resultList}->{result}->{authorString});
    $ref->{resultList}->{result}->{authorString} = $trimmed_author_list if $trimmed_author_list;

    return $ref;
}

# Format authors list
# only store the first 4 authors
sub trim_author_list{
  my $authors_final = shift;

  my $trimmed_authors;

  if(defined $authors_final){
    my @author_list = split /, /, $authors_final;
    if(scalar @author_list > 4) {
      $trimmed_authors = join(', ', @author_list[0..3]) . ', et al';
    }
  }
  return $trimmed_authors;
}

sub get_epmc_data{

    my $id = shift; ## specific part of URL including pmid or pmcid

    return undef unless defined $id && $id =~/\d+/;

    my $xs   = XML::Simple->new();
    my $server = 'https://www.ebi.ac.uk/europepmc/';
    my $request  = $server . $id . '&format=xml';

    my %data;

    print "Looking for $request\n\n"  if $DEBUG == 1;
    my $response = $http->get($request, {
        headers => {
         'Accept' => 'application/xml'
        }
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

    if(defined $mined->{item}->{annotations}->{annotation}  &&
      ref($mined->{item}->{annotations}->{annotation}) eq 'ARRAY'){
      foreach my $found (@{$mined->{item}->{annotations}->{annotation}}) {
        if ($found->{exact} =~ /$species_string/i){
          $looks_ok = 1;
        }
      }
    }
    elsif(defined $mined->{item}->{annotations}->{annotation} &&
      ref($mined->{item}->{annotations}->{annotation}) eq 'HASH'){
      foreach my $key (keys %{$mined->{item}->{annotations}->{annotation}}){
        if(ref($mined->{item}->{annotations}->{annotation}->{$key}) eq 'HASH' &&
          $mined->{item}->{annotations}->{annotation}->{$key}->{exact}){
          if($mined->{item}->{annotations}->{annotation}->{$key}->{exact} =~ /$species_string/i){
            $looks_ok = 1;
          }
        }
        if($key eq 'exact' && $mined->{item}->{annotations}->{annotation}->{$key} =~ /$species_string/i){
          $looks_ok = 1;
        }
      }
    }
    return $looks_ok;
}

## read input
##
## expected format:  refsnp,PMCID,EXTID,SOURCE
##                   "rs1914711",PMC6909262,31872004,MED

sub parse_EPMC_file{

    my $pubmed_file = shift;
    my $avoid_list  = shift;

    my %data;
    open (my $file, $pubmed_file) or die "Failed to open file of PubMed ids $pubmed_file : $!\n";
    while (<$file>){
      ## basic check that format is as expected
      next unless /\"rs\d+\"\,/i;
      chomp;
      s/"//g;

      my ($rs, $pmcid, $pmid, $source) = split/\,/;

      next if($source ne "MED");

      my $is_pmcid;
      if ($pmcid =~ /PMC\d+/) {
        $is_pmcid = 1;
      }

      ## remove known errors
      if ($is_pmcid && ($avoid_list->{$rs}->{$pmcid} || $avoid_list->{$rs}->{$pmid})){
        print "Removing suspected error : $rs in $pmcid\n";
        next;
      }
      ## Group by publication; may have pmid, pmcid or both
      my $tag = $pmcid . "_" .$pmid;

      $data{$tag}{pmid} = $pmid;
      $data{$tag}{pmid} = undef unless $pmid =~ /\d+/;
      $data{$tag}{pmcid} = $is_pmcid ? $pmcid : undef;

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

    my $dba = shift;

    open (my $error_log, ">>$species\_dbSNP.log") or die "Failed to open log file:$!\n";

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

        my $mined = get_epmc_data( "annotations_api/annotationsByArticleIds?articleIds=MED:$l->[1]&type=Organisms&format=XML" );
        my $ref   = get_epmc_data("webservices/rest/search?query=ext_id:$l->[1]%20src:med");
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

        # format authors list
        my $trimmed_author_list = trim_author_list($ref->{resultList}->{result}->{authorString});
        $ref->{resultList}->{result}->{authorString} = $trimmed_author_list if $trimmed_author_list;

        if (defined $ref->{resultList}->{result}->{title}) {
          my $title = $ref->{resultList}->{result}->{title};
          # Truncate title
          if(length($title) >= 300){
            my $aux = substr($title, 0, 296);
            my @list_title = split(' ', $aux);
            pop @list_title;
            $title = join(' ', @list_title);
            $title .= '...';
          }

          $title =~ s|<.+?>||g;

          $pub_upd_sth->execute( $title,
                                 $ref->{resultList}->{result}->{pmcid},
                                 unidecode($ref->{resultList}->{result}->{authorString}),
                                 $ref->{resultList}->{result}->{pubYear},
                                 $ref->{resultList}->{result}->{doi},
                                 $l->[0]
          );
        }
    }
    close $error_log;
}

## citation is one form of evidence used to support a variant
## add this evidence attrib to the variation & variation_feature records
sub update_evidence{

    my $dba = shift;

    ## find cited attrib
    my $attrib_ext_sth = $dba->dbc()->prepare(qq[ select attrib_id from attrib where value ='Cited']);
    $attrib_ext_sth->execute() or die "Failed to select attrib_id from attrib where value 'Cited'\n";
    my $attrib =  $attrib_ext_sth->fetchall_arrayref();
    die "Not updating evidence as no attrib found\n" unless defined $attrib->[0]->[0];

    my $ev_ext_sth = $dba->dbc()->prepare(qq[ select variation.variation_id, variation.evidence_attribs 
                                              from variation, variation_citation
                                              where variation.variation_id = variation_citation.variation_id
      ]);

    my $var_upd_sth     = $dba->dbc()->prepare(qq[ update variation set evidence_attribs = ? where variation_id = ?]);
    my $varfeat_upd_sth = $dba->dbc()->prepare(qq[ update variation_feature set evidence_attribs = ? where variation_id = ?]);
    
    $ev_ext_sth->execute() or die "Failed to select variation_id and evidence_attribs from table variation_citation\n";
    my $dat =  $ev_ext_sth->fetchall_arrayref();
    
    my $n = scalar @{$dat};
    print "\n$n variants with citation evidence\n";

    foreach my $l (@{$dat}){

        next if defined $l->[1] && $l->[1] =~ /$attrib->[0]->[0]/;

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

    my $title_null; 

    open (my $report, ">import_EPMC_$species\_"  . log_time() . ".txt") or die "Failed to open report file to write: $!\n";

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

    # Query is taking a few hours to run
    my $dup2_ext_sth = $dba->dbc->prepare(qq[ select p1.publication_id, p2.publication_id, p2.pmcid
                                         from publication p1, publication p2
                                         where p1.pmcid = p2.pmcid
                                         and p1.publication_id < p2.publication_id
                                         and p1.pmcid is not null
                                         ]);

   # Query is taking a few hours to run
   my $dup3_ext_sth = $dba->dbc->prepare(qq[ select p1.publication_id, p2.publication_id, p2.doi
                                         from publication p1, publication p2
                                         where p1.doi = p2.doi
                                         and p1.publication_id < p2.publication_id
                                         and p1.doi is not null
                                         ]);

    my $dup4_ext_sth = $dba->dbc->prepare(qq[ select title, authors, count(*)
                                              from publication 
                                              group by upper(title), upper(authors) having count(*) > 1
                                          ]);
 
    my $fail_ext_sth = $dba->dbc->prepare(qq[ select count(*) from publication
                                              where title is null
                                          ]);

    $dup1_ext_sth->execute() or die;
    my $dup1 = $dup1_ext_sth->fetchall_arrayref();
    my $duplicated_pub = $dup1->[0]->[0];

    $dup2_ext_sth->execute() or die;
    my $dup2 = $dup2_ext_sth->fetchall_arrayref();
    my $duplicated_pub2 = $dup2->[0]->[0];

    $dup3_ext_sth->execute() or die;
    my $dup3 = $dup3_ext_sth->fetchall_arrayref();
    my $duplicated_pub3 = $dup3->[0]->[0];

    $dup4_ext_sth->execute() or die;
    my $dup4 = $dup4_ext_sth->fetchall_arrayref();
    my $duplicated_pub4 = $dup4->[0]->[0];

    $fail_ext_sth->execute() or die;
    my $fail = $fail_ext_sth->fetchall_arrayref();
    $title_null = $fail->[0]->[0]; 

    if (defined $duplicated_pub){
      print $report "\nDuplicated publications (publication 1, publication 2, pmid):\n";
      foreach my $l (@{$dup1}){
        print $report "$l->[0]\t$l->[1]\t$l->[2]\n";
      }
    }
    if (defined $duplicated_pub2){
      print $report "\nDuplicated publications (publication 1, publication 2, pmcid):\n";
      foreach my $k (@{$dup2}){
        print $report "$k->[0]\t$k->[1]\t$k->[2]\n";
      }
    }
    if (defined $duplicated_pub3){
    print $report "\nDuplicated publications (publication 1, publication 2, doi):\n";
      foreach my $m (@{$dup3}){
        print $report "$m->[0]\t$m->[1]\t$m->[2]\n";
      }
    }
    if (defined $duplicated_pub4){
      print $report "\nDuplicated publications (title, authors, number of publications):\n";
      foreach my $i (@{$dup4}){
        print $report "$i->[0]\t$i->[1]\t$i->[2]\n"; 
      }
    }

    print $report "\n$title_null publications without a title - to be deleted\n";

    return $title_null;
}

## store rundate in meta table for each source 
sub update_meta{

    my $dba  = shift;
    my $type = shift;

    my $meta_ins_sth = $dba->dbc->prepare(qq[ insert ignore into meta
                                              ( meta_key, meta_value) values (?,?)
                                            ]);


    $meta_ins_sth->execute($type . "_citation_update", run_date() );
}

## clean publication and variation_citation tables after import
## Clean brackets from titles
## Delete publications with unspecific titles including titles 'Not Available'
## Delete publications without title
sub clean_publications{

    my $dba  = shift; 
    my $n_title_null = shift;

    open (my $report, ">import_EPMC_$species\_clean_publications_"  . log_time() . ".txt") or die "Failed to open report file to write: $!\n";

    my $publication_count = count_rows($dba, 'publication');
    my $citation_count    = count_rows($dba, 'variation_citation');

    print $report "Total publications before cleaning:\t$publication_count\n";
    print $report "Total citations before cleaning:\t$citation_count\n\n";

    my $pub_ad = $reg->get_adaptor($species, 'variation', 'publication');

    my $title_sth = $dba->dbc->prepare(qq[ select publication_id, title from publication where title like '\[%' and (title like '%\]' or title like '%\.\]' or title like '%\]\.') ]);

    my $title_cr_sth = $dba->dbc->prepare(qq[ select publication_id, title from publication where title like '%\n%' ]);
    my $authors_cr_sth = $dba->dbc->prepare(qq[ select publication_id, authors from publication where authors like '%\n%' ]);

    my $title_hex_char_sth = $dba->dbc->prepare(qq[ select publication_id, title from publication where title like '%&#x%' ]);
    my $authors_hex_char_sth = $dba->dbc->prepare(qq[ select publication_id, authors from publication where authors like '%&#x%' ]);

    my $wrong_title_sth = $dba->dbc->prepare(qq[ select publication_id, title from publication where title like '%Errata%'
                   or title like '%Erratum%' or title like '%In This Issue%' or title like '%Oral abstracts%' 
                   or title like '%Oral Presentations%' or title like '%Proffered paper%' or title like '%Subject Index%' 
                   or title like '%Summaries of Key Journal Articles%' or title like '%This Month in The Journal%' 
                   or title like 'Index%' or title like '%Table of Contents%' or title like '%Not Available%' 
                   or title like 'Beyond Our Pages%' or title like 'EP News%' or title like 'ACTS Abstracts%' 
                   or title like 'Poster %' or title like 'Abstracts.%' or title like 'Abstract.%' ]);

    my $empty_sth = $dba->dbc->prepare(qq[ select publication_id, title from publication where (authors = '' or authors is null) and pmid is null and pmcid is null ]);

    $title_sth->execute() or die;
    my $title_brackets = $title_sth->fetchall_arrayref();
    
    $title_cr_sth->execute() or die;
    my $get_title_cr = $title_cr_sth->fetchall_arrayref();
    my $title_cr = $get_title_cr->[0]->[0];

    $authors_cr_sth->execute() or die;
    my $get_authors_cr = $authors_cr_sth->fetchall_arrayref();
    my $authors_cr = $get_authors_cr->[0]->[0];

    $title_hex_char_sth->execute() or die;
    my $get_title_hex_char = $title_hex_char_sth->fetchall_arrayref();
    my $title_hex_char = $get_title_hex_char->[0]->[0];

    $authors_hex_char_sth->execute() or die;
    my $get_authors_hex_char = $authors_hex_char_sth->fetchall_arrayref();
    my $authors_hex_char = $get_authors_hex_char->[0]->[0];

    $wrong_title_sth->execute() or die;
    my $wrong_title = $wrong_title_sth->fetchall_arrayref();

    $empty_sth->execute() or die;
    my $empty_fields = $empty_sth->fetchall_arrayref();

    # Clean brackets from publication title 
    my $count_titles = 0;
    if (defined $title_brackets->[0]->[0]){
      print $report "Publications with brackets in the title - clean the brackets:\n";
      foreach my $l (@{$title_brackets}){
        $count_titles += 1;
        my $pub_id = $l->[0];
        my $title = $l->[1];
        print $report "$pub_id\t$title\n";
        $title =~ s/^\[//;
        $title =~ s/\]//;

        my $pub_clean_sth = $dba->dbc()->prepare(qq[ update publication set title = ? where publication_id = $pub_id ]);
        $pub_clean_sth->execute($title);
      }
    }

    # Clean carriage returns in title
    if (defined $title_cr){
      print $report "\nPublications with carriage return in the title - clean:\n";
      foreach my $t_cr (@{$get_title_cr}){
        my $pub_id = $t_cr->[0];
        my $title = $t_cr->[1];
        print $report "$pub_id\t$title\n";
        $title =~ s/\n/ /g;
        
        my $clean_title_sth = $dba->dbc()->prepare(qq[ update publication set title = ? where publication_id = $pub_id ]);
        $clean_title_sth->execute($title);
      }
    }
    
    # Clean carriage returns in authors
    if (defined $authors_cr){
      print $report "\nPublications with carriage return in the authors - clean:\n";
      foreach my $a_cr (@{$get_authors_cr}){
        my $pub_id = $a_cr->[0];
        my $authors = $a_cr->[1];
        print $report "$pub_id\t$authors\n";
        $authors =~ s/\n/ /g;
        
        my $clean_authors_sth = $dba->dbc()->prepare(qq[ update publication set authors = ? where publication_id = $pub_id ]);
        $clean_authors_sth->execute($authors);
      }
    }

    # Replace hexadecimal characters in title
    if (defined $title_hex_char){
      print $report "\nPublications with hexadecimal characters in the title - clean:\n";
      foreach my $t_hex_char (@{$get_title_hex_char}){
        my $pub_id = $t_hex_char->[0];
        my $title = $t_hex_char->[1];
        my $new_title = replace_hex($title);

        # Double check if title still contains hexadecimal characters
        if($new_title =~ /&#x/) {
          print $report "$pub_id\t$title\tWARNING: Title with hexadecimal characters\n";
        }
        # Continue update if title string is clean
        else {
          print $report "$pub_id\t$title\n";
          my $clean_title_sth = $dba->dbc()->prepare(qq[ update publication set title = ? where publication_id = $pub_id ]);
          $clean_title_sth->execute($new_title);
        }
      }
    }

    # Replace hexadecimal characters in authors
    if (defined $authors_hex_char){
      print $report "\nPublications with hexadecimal characters in the authors - clean:\n";
      foreach my $a_hex_char (@{$get_authors_hex_char}){
        my $pub_id = $a_hex_char->[0];
        my $authors = $a_hex_char->[1];
        my $new_authors = replace_hex($authors);

        # Double check if authors still contain hexadecimal characters
        if($new_authors =~ /&#x/) {
          print $report "$pub_id\t$authors\tWARNING: Authors with hexadecimal characters\n";
        }
        # Continue update if authors string is clean
        else {
          print $report "$pub_id\t$authors\n";
          my $clean_authors_sth = $dba->dbc()->prepare(qq[ update publication set authors = ? where publication_id = $pub_id ]);
          $clean_authors_sth->execute($new_authors);
        }
      }
    }

    # Delete publications with not acceptable titles
    if(defined $wrong_title->[0]->[0]){
      print $report "\nDeleted publications with title not acceptable (publication_id, title):\n";
      remove_publications($pub_ad, $report, $wrong_title);
    }

    # If there is publications without title, delete them 
    if($n_title_null != 0){
      my $title_null_sth = $dba->dbc->prepare(qq[ select publication_id, pmid from publication where title is null or title = '' ]);
      $title_null_sth->execute();
      my $titles_null = $title_null_sth->fetchall_arrayref();

      if(defined $titles_null->[0]->[0]){
        print $report "\nPublications without title deleted (publication_id, title):\n";
        remove_publications($pub_ad, $report, $titles_null);
      }
    }

    # Delete publications without authors, pmid and pmcid
    if(defined $empty_fields->[0]->[0]){
      print $report "\nDeleted publications with empty authors, pmid and pmcid (publication_id, title)";
      remove_publications($pub_ad, $report, $empty_fields);
    }

    $publication_count = count_rows($dba, 'publication');
    $citation_count    = count_rows($dba, 'variation_citation');

    print $report "\nTotal publications after cleaning:\t$publication_count\n";
    print $report "Total citations after cleaning:\t$citation_count\n";
}

# Delete publication from publication and variation_citation tables
sub remove_publications{
  my $pub_ad = shift;
  my $report = shift;
  my $publications = shift;

  foreach my $p (@{$publications}){
    my $pub_id = $p->[0];
    my $title  = $p->[1];
    $pub_ad->remove_publication_by_dbID($pub_id) or
      die "ERROR while removing publication $pub_id : $title\n";
    print $report "$pub_id\t$title\n";
  }
}

## export current snapshot from database
sub get_current_UCSC_data{

    my $filename = "UCSC_export";
    open (my $out, ">$filename") or die "Failed to open file to write data : $!\n";

    my $dbh = DBI->connect('dbi:mysql:hgFixed:genome-mysql.cse.ucsc.edu:3306:max_allowed_packet=1MB', 'genome', '', undef);

    my $cit_ext_sth = $dbh->prepare(qq[ SELECT pma.markerId, pa.pmid, pma.section, pa.doi, pa.title, pa.authors, pa.year, pa.extId
                                        FROM pubsMarkerAnnot pma JOIN pubsArticle pa USING
                                        (articleId) WHERE pma.markerType="snp"
                                      ]);


    $cit_ext_sth->execute() or die "Failed to select data from UCSC database\n";

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
    open (my $file, $ucsc_file) or die "Failed to open file of UCSC citations $ucsc_file : $!\n";
    while (<$file>){

        next unless /rs\d+/i;
        chomp;
        s/RS/rs/;
        s/^\|\s+//;

        #my ($rs, $pmid, $section, $doi, $title, $authors, $year ) = split/\s+\|\s+|\t/;
        my ($rs, $pmid, $section, $doi, $title, $authors, $year, $extId ) = split/\t/;
        next if $section =~/refs|ack/; ## not in this publication

        # format authors list
        my $trimmed_author_list = trim_author_list($authors);
        $authors = $trimmed_author_list if $trimmed_author_list;

        $pmid = "" if $pmid eq "0";  ## nulls are set to 0 in UCSC database

        ## remove known errors
         if ( $avoid_list->{$rs}->{$pmid}){
            print "Removing suspected error : $rs in $pmid\n";
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
        $data{$tag}{pmcid}   = undef;
        if($extId =~ /^PMC\d+/) {
          $data{$tag}{pmcid} = $extId;
        }
        $data{$tag}{ucsc}    = $extId;

        push @{$data{$tag}{rsid}}, $rs;
    }

    return \%data;
}

## find current citation in db
## and skip these lines in the file.
sub get_current_citations{

    my $dba = shift;

    my $cit_ext_sth = $dba->dbc()->prepare(qq[ select variation.name, publication.pmid, variation_citation.data_source_attrib
                                               from publication, variation, variation_citation
                                               where publication.publication_id = variation_citation.publication_id
                                               and variation_citation.variation_id = variation.variation_id
                                              ]);

    my %citations;

    $cit_ext_sth->execute() or die "Failed to select current citations from database\n";
    my $data =  $cit_ext_sth->fetchall_arrayref();
    foreach my $l(@{$data}){
        $citations{$l->[0]}{$l->[1]} = $l->[2] if defined $l->[1];
    }

    return \%citations;
}

###
# Get citations from phenotype tables

sub process_phenotype_feature {
  my $dba = shift;
  my $source_ad = shift;
  my $citation_attribs = shift;

  my %list_citations_pheno_feature;

  ## Get citations from phenotype_feature
  my $pheno_citations_sth = $dba->dbc()->prepare(qq[ select s.study_id, s.source_id, s.external_reference, s.study_type, p.object_id
                                                from study s
                                                inner join phenotype_feature p on s.study_id = p.study_id
                                                where p.type = 'variation' and p.study_id is not null and s.external_reference is not null
                                                group by s.study_id, s.source_id, s.external_reference, s.study_type, p.object_id ]);
  $pheno_citations_sth->execute() or die "Failed to select citations from table phenotype_feature\n";
  my $data = $pheno_citations_sth->fetchall_arrayref();

  my %source_id_list;

  foreach my $pheno_data (@{$data}){
    my $source_id = $pheno_data->[1];
    my $study_type = $pheno_data->[3];

    next if(defined $study_type);

    if(!defined $source_id_list{$source_id}) {
      my $source_obj = $source_ad->fetch_by_dbID($source_id);
      my $source_name = $source_obj->name();

      my $source_attrib_id = $citation_attribs->{$source_name};
      die "No attrib of type 'citation_source' was found for '$source_name'!\n" unless defined $source_attrib_id;

      $source_id_list{$source_id} = $source_attrib_id;
    }
  }

  foreach my $l (@{$data}){
    my $study_id = $l->[0];
    my $source_id = $l->[1];
    my $external_reference = $l->[2];
    $external_reference =~ s/PMID://;
    my $study_type = $l->[3];
    my $rsid = $l->[4];

    # Get attrib id for source - some are null 
    my $source_attrib_id;
    if(defined $study_type){
      $source_attrib_id = $citation_attribs->{$study_type};
      die "No attrib of type 'citation_source' was found for '$study_type'!\n" unless defined $source_attrib_id;
    }
    else{
      # Get source name from source table (dbGaP)
      $source_attrib_id = $source_id_list{$source_id};
    }

    my $tag = $source_attrib_id . "_" . $external_reference;

    $list_citations_pheno_feature{$tag}{pmid} = $external_reference;

    $list_citations_pheno_feature{$tag}{pmcid} = undef;
    push @{$list_citations_pheno_feature{$tag}{rsid}}, $rsid;

    $list_citations_pheno_feature{$tag}{source} = $source_attrib_id;
  }

  return \%list_citations_pheno_feature;
}

sub process_phenotype_feature_attrib {
  my $dba = shift;
  my $source_ad = shift;
  my $citation_attribs = shift;

  my %list_citations_pheno_feature_attrib;

  my $pheno_feature_sth = $dba->dbc()->prepare(qq[ select pfa.phenotype_feature_id, pfa.value, pf.source_id, pf.object_id
                                                   from phenotype_feature_attrib pfa
                                                   inner join phenotype_feature pf on pfa.phenotype_feature_id = pf.phenotype_feature_id
                                                   join attrib_type att on pfa.attrib_type_id = att.attrib_type_id
						   left join source s on s.source_id = pf.source_id
                                                   where att.code = 'pubmed_id' and s.name = 'ClinVar' ]);

  $pheno_feature_sth->execute() or die "Failed to select citations from table phenotype_feature_attrib\n";
  my $pheno_feature_data = $pheno_feature_sth->fetchall_arrayref();

  my %source_id_list;

  foreach my $pheno_feat_data (@{$pheno_feature_data}){
    my $source_id = $pheno_feat_data->[2];

    if(!defined $source_id_list{$source_id}) {
      my $source_obj = $source_ad->fetch_by_dbID($source_id);
      my $source_name = $source_obj->name();
      my $source_attrib_id = $citation_attribs->{$source_name};
      die "No attrib of type 'citation_source' was found for '$source_name'!\n" unless defined $source_attrib_id;

      $source_id_list{$source_id} = $source_attrib_id;
    }
  }

  # Create map to associate PMID with all phenotype features it's linked to
  foreach my $pheno_feat_data (@{$pheno_feature_data}){
    # my $pheno_feat_id = $pheno_feat_data->[0];
    my $value_pubid = $pheno_feat_data->[1];
    my $source_id = $pheno_feat_data->[2];
    my $var_name = $pheno_feat_data->[3];

    my @value_pubid_splited = split /,/, $value_pubid;

    # Get publication with same PMID from study table
    foreach my $pmid (@value_pubid_splited){
      # PMID=25806920 is the same paper as PMID=25806919
      # PMID=25806919 is already in the publication table - faster to skip this PMID and avoid duplicated data
      next if $pmid eq '25806920';

      # Get source name from source table (ClinVar)
      my $source_attrib_id = $source_id_list{$source_id};

      my $tag = $source_attrib_id . "_" . $pmid;

      $list_citations_pheno_feature_attrib{$tag}{pmid} = $pmid;

      $list_citations_pheno_feature_attrib{$tag}{pmcid} = undef;
      push @{$list_citations_pheno_feature_attrib{$tag}{rsid}}, $var_name;

     $list_citations_pheno_feature_attrib{$tag}{source} = $source_attrib_id;
    }
  }

  return \%list_citations_pheno_feature_attrib;
}

# To avoid outdated citations, we remove citations from ClinVar, GWAS or dbGaP that are not in the Phenotype tables.
sub remove_outdated_citations {
  my $dba = shift;
  my $var_ad = shift;
  my $pub_ad = shift;
  my $citation_attribs = shift;
  my $citations_pheno_feature = shift;
  my $citations_pheno_feature_attrib = shift;

  open (my $wrt, ">Outdated_Phenotype_citations_$species\_"  . log_time() . ".txt") or die "Failed to open file to write: $!\n";
  print $wrt "RSID\tPMID\tSource\n";

  # get all citations from the sources 'ClinVar', 'dbGaP', 'GWAS' and 'G2P' - imported from the phenotype tables
  my $attrib_id_clinvar = $citation_attribs->{'ClinVar'};
  my $attrib_id_gwas = $citation_attribs->{'GWAS'};
  my $attrib_id_dbgap = $citation_attribs->{'dbGaP'};
  my $attrib_id_ddg2p = $citation_attribs->{'G2P'};
  my @attrib_ids = ($attrib_id_clinvar, $attrib_id_gwas, $attrib_id_dbgap, $attrib_id_ddg2p);

  my $citations_sth = $dba->dbc()->prepare(qq[
    select variation_id, publication_id, data_source_attrib
    from variation_citation
    where data_source_attrib like '%$attrib_id_clinvar%'
       or data_source_attrib like '%$attrib_id_gwas%'
       or data_source_attrib like '%$attrib_id_dbgap%'
       or data_source_attrib like '%$attrib_id_ddg2p%' ]);

  $citations_sth->execute() or die "Failed to fetch outdated citations from database\n";
  my $citations_data = $citations_sth->fetchall_arrayref();

  my $rm_citations_sth = $dba->dbc()->prepare(qq[
        delete from variation_citation
        where data_source_attrib = ? and
              variation_id = ? and 
              publication_id = ?
  ]);
        
  my $update_citations_sth = $dba->dbc()->prepare(qq[
        update variation_citation
        set data_source_attrib = ?
        where variation_id = ? and 
              publication_id = ?
  ]);

  my $cited_attrib_id = $dba->get_AttributeAdaptor->attrib_id_for_type_value('evidence', 'Cited');
  my $rm_cited_evidence = qq[
        update %s v
        set v.evidence_attribs = NULLIF(TRIM(BOTH ',' FROM REPLACE(
            CONCAT(',', v.evidence_attribs, ','), ',$cited_attrib_id,', ',')), '')
        where v.variation_id = ?
          and not exists (select * from variation_citation vc
                          where vc.variation_id = v.variation_id)
  ];
  my $rm_variation_cited_evidence_sth = $dba->dbc()->prepare(sprintf $rm_cited_evidence, 'variation');
  my $rm_variation_feature_cited_evidence_sth = $dba->dbc()->prepare(sprintf $rm_cited_evidence, 'variation_feature');

  foreach my $c (@{$citations_data}){
    my $variation_id = $c->[0];
    my $publication_id = $c->[1];
    my $attrib_id = $c->[2];

    my $variation = $var_ad->fetch_by_dbID($variation_id);
    die "No variation found for '$variation_id'!\n" unless defined $variation;

    my $publication = $pub_ad->fetch_by_dbID($publication_id);
    die "No publication found for '$publication_id'!\n" unless defined $publication;

    my $variation_rsid = $variation->name();
    my $publication_pmid = $publication->pmid();

    my @split_attrib_id = split /,/, $attrib_id;

    my @outdated_attribs;
    my @current_attribs;
    foreach my $attrib (@split_attrib_id) {
      my $is_selected_attrib_id = grep(/^$attrib$/, @attrib_ids) ? 1 : 0;

      my $value = $attrib . '_' . $publication_pmid;
      my $is_cited_in_pheno = $citations_pheno_feature->{$value} || $citations_pheno_feature_attrib->{$value};
      if (!$is_selected_attrib_id || $is_cited_in_pheno) {
        push @current_attribs, $attrib;
      } else {
        push @outdated_attribs, $attrib;
        print $wrt "$variation_rsid\t$publication_pmid\t$attrib\n";
      }
    }

    if (@outdated_attribs eq @split_attrib_id) {
      # remove citation if outdated in all sources
      $rm_citations_sth->execute($attrib_id, $variation_id, $publication_id) or
        die "Error: cannot remove $variation_id, $publication_id, $attrib_id from variation_citation\n";
      # remove 'Cited' evidence if variant is not in variation_citation
      $rm_variation_cited_evidence_sth->execute($variation_id) or
        die "Error: cannot remove 'Cited' evidence for $variation_id in variation";
      $rm_variation_feature_cited_evidence_sth->execute($variation_id) or
        die "Error: cannot remove 'Cited' evidence for $variation_id in variation_feature";
    } elsif (@outdated_attribs) {
      # discard outdated sources
      my $new_attrib_id = join(",", @current_attribs);
      $update_citations_sth->execute($new_attrib_id, $variation_id, $publication_id) or
        die "Error: cannot update $variation_id, $publication_id, $attrib_id from variation_citation with new value $new_attrib_id\n";
    }
  }

  # Remove publications not associated with any variants
  my $orphan_publications_sth = $dba->dbc()->prepare(qq[
    select publication_id, title
    from publication
    where publication_id NOT IN
      (select distinct publication_id from variation_citation) ]);
  $orphan_publications_sth->execute() or die;
  my $orphan_publications = $orphan_publications_sth->fetchall_arrayref();

  if ($orphan_publications) {
    print $wrt "\nDeleted publications with no associated variants (publication_id, title):\n";
    remove_publications($pub_ad, $wrt, $orphan_publications);
  } else {
    print $wrt "\nNo publications with no associated variants to remove.\n";
  }

  close($wrt);
}

# Fetch all attribs that have type citation_source
# Method used by type = 'phenotype'
sub get_citation_attribs {
  my $dba = shift;

  my %attribs;

  my $attrib_value_sth = $dba->dbc()->prepare(qq[ SELECT a.attrib_id, a.value FROM attrib a
                                                  JOIN attrib_type at ON a.attrib_type_id = at.attrib_type_id
                                                  WHERE at.code = 'citation_source' ]);
  $attrib_value_sth->execute();
  my $citation_attribs = $attrib_value_sth->fetchall_arrayref();

  die "No attribute of type 'citation_source' was found in attrib table!\n" unless defined $citation_attribs->[0];

  foreach my $data (@{$citation_attribs}){
    my $id = $data->[0];
    my $name = $data->[1];
    $attribs{$name} = $id;
  }

  return \%attribs;
}

sub usage {
    
    die "\n\tUsage: import_EPMC.pl -type [ EPMC or UCSC or phenotype] -species [name] -registry [registry file]

Options:  
          -data [file of citations]   - *required* for EPMC import
          -check_dbSNP [0/1]          - add detail to citations from dbSNP import (default:1)
          -clean [0/1]                - clean publications after import (default:0)
          -no_evidence                - don't update variation & variation_feature evidence statuses\n\n";

}

sub unique {
  my %a;
  map { $a{$_} = 1; } @_;
  return sort keys %a;
}
