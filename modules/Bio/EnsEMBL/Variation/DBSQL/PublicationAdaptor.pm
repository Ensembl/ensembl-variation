=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut



=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::PublicationAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $pa = $reg->get_adaptor("human", "variation","publication");

  # fetch a publication by its PubMed ID
  $publication = $sta->fetch_by_pmid(17554300); 

  
  

=head1 DESCRIPTION

This adaptor provides database connectivity for Publication objects.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::PublicationAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Publication;

use base qw{Bio::EnsEMBL::DBSQL::BaseAdaptor};




=head2 fetch_by_pmid

  Arg [1]    : string $pmid
  Example    : $publication = $publication_adaptor->fetch_by_pmid(17554300);
  Description: Retrieves a publication object via its pubmed id
  Returntype : Bio::EnsEMBL::Variation::Publication
  Exceptions : throw if name argument is not defined
  Caller     : general
  Status     : At Risk

=cut

sub fetch_by_pmid {
  my $self = shift;
  my $pmid = shift;

  throw('pmid argument expected') if(!defined $pmid );

  my $result = $self->generic_fetch("p.pmid = $pmid");

  return ($result ? $result->[0] : undef);
}

=head2 fetch_by_pmcid

  Arg [1]    : string $pmcid
  Example    : $publication = $publication_adaptor->fetch_by_pmcid(17554300);
  Description: Retrieves a publication object via its PubMed Central id
  Returntype : Bio::EnsEMBL::Variation::Publication
  Exceptions : throw if name argument is not defined
  Caller     : general
  Status     : At Risk

=cut

sub fetch_by_pmcid {
  my $self = shift;
  my $pmcid = shift;

  throw('pmcid argument expected') if(!defined $pmcid );

  my $result = $self->generic_fetch("p.pmcid = \'$pmcid\'");

  return ($result ? $result->[0] : undef);
}

=head2 fetch_by_doi

  Arg [1]    : string $doi
  Example    : $publication = $publication_adaptor->fetch_by_doi("10.2337/db08-0569");
  Description: Retrieves a Publication object via its doi identifier.
               If no such publication exists undef is returned.
  Returntype : Bio::EnsEMBL::Variation::Publication
  Exceptions : throw if dbID arg is not defined
  Caller     : general
  Status     : At Risk

=cut

sub fetch_by_doi {
  my $self = shift;
  my $doi = shift;

  throw('doi argument expected') if(!defined($doi));

  my $result = $self->generic_fetch("p.doi=\"$doi\"");
 
  return ($result ? $result->[0] : undef);
}

=head2 fetch_by_dbID

  Arg [1]    : int $dbID
  Example    : $publication = $publication_adaptor->fetch_by_dbID(254);
  Description: Retrieves a Publication object via its internal identifier.
               If no such publication exists undef is returned.
  Returntype : Bio::EnsEMBL::Variation::Publication
  Exceptions : throw if dbID arg is not defined
  Caller     : general
  Status     : At Risk

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  throw('dbID argument expected') if(!defined($dbID));

  my $result = $self->generic_fetch("p.publication_id=$dbID");
 
  return ($result ? $result->[0] : undef);
}

        
=head2 fetch_all_by_dbID_list

  Arg [1]    : listref $list
  Example    : $publication = $publication_adaptor->fetch_all_by_dbID_list([907,1132]);
  Description: Retrieves a listref of publication objects via a list of internal
               dbID identifiers
  Returntype : listref of Bio::EnsEMBL::Variation::Publication objects
  Exceptions : throw if list argument is not defined
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_dbID_list {
    my $self = shift;
    my $list = shift;
    
    if(!defined($list) || ref($list) ne 'ARRAY') {
        throw("list reference argument is required");
    }
    
    my $id_str = (@$list > 1)  ? " IN (".join(',',@$list).")"   :   ' = \''.$list->[0].'\'';
    
    my $result = $self->generic_fetch("p.publication_id $id_str");
    
    return ($result ? $result : undef);
}

=head2 fetch_all_by_Variation

  Arg [1]    : listref $list
  Example    : $publication = $publication_adaptor->fetch_all_by_Variation( $var_object]);
  Description: Retrieves a listref of publication objects via a variation object
  Returntype : listref of Bio::EnsEMBL::Variation::Publication objects
  Exceptions : throw if variation argument is not defined
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_Variation {
    my $self    = shift;
    my $var_obj = shift;
    
    if(!defined($var_obj) || ref($var_obj ne 'Variation')) {
        throw("variation argument is required");
    }
    
    my @pub;
    my $publication_id;

    my $sth = $self->prepare(qq{SELECT  publication_id from variation_citation where variation_id = ?  });
    $sth->execute($var_obj->dbID);
    $sth->bind_columns(\$publication_id);
    while ($sth->fetch()){
        push @pub, $self->fetch_by_dbID($publication_id)
    }
    $sth->finish;
    
    return \@pub;

}



sub _columns {
  return qw(p.publication_id p.title p.authors p.pmid p.pmcid p.year p.doi p.ucsc_id );
}

sub _tables { return (['publication', 'p']); }

sub _default_where_clause {
  my $self = shift;
  return ;
}

#
# private method, creates publication objects from an executed statement handle
# ordering of columns must be consistant
#
sub _objs_from_sth {
  my $self = shift;
  my $sth  = shift;

  my @publication;

  my ($publication_id,$title,$authors,$pmid,$pmcid,$year,$doi,$ucsc_id);
  
  $sth->bind_columns(\$publication_id, \$title, \$authors, \$pmid, \$pmcid, \$year, \$doi, \$ucsc_id);
  
  while($sth->fetch()) {
                
                
      push @publication, Bio::EnsEMBL::Variation::Publication->new
          (-dbID       => $publication_id,
           -ADAPTOR    => $self,
           -TITLE      => $title,
           -AUTHORS    => $authors,
           -PMID       => $pmid,
           -PMCID      => $pmcid,
           -YEAR       => $year,
           -DOI        => $doi,
           -UCSC_ID    => $ucsc_id
          );
  }

  return \@publication;
}


sub store{

    my $self = shift;
    my $pub  = shift;
    
    my $dbh = $self->dbc->db_handle;
    
    my $pub_ins_sth = $dbh->prepare(qq[ insert into publication( title, authors, pmid, pmcid, year, doi, ucsc_id ) values ( ?,?,?,?,?,?,? ) ]);   
    
    $pub_ins_sth->execute( $pub->{title},
                           $pub->{authors},
                           $pub->{pmid}    || undef,
                           $pub->{pmcid}   || undef,
                           $pub->{year},
                           $pub->{doi}     || undef, 
                           $pub->{ucsc_id} || undef, 
                         );


    $pub->{dbID} = $dbh->last_insert_id(undef, undef, 'publication', 'publication_id');

    ## link cited variants to publication
    $self->update_variant_citation($pub) if defined $pub->{variants}->[0] ;

}

sub update_variant_citation {

    my $self = shift;       
    my $pub  = shift;
    my $var  = shift;   

    my $dbh = $self->dbc->db_handle;

    my $citation_ext_sth = $dbh->prepare(qq[ select count(*) from variation_citation where variation_id =? and publication_id =?  ]);   
    
    my $citation_ins_sth = $dbh->prepare(qq[ insert into variation_citation( variation_id, publication_id) values ( ?,?) ]);   
    
    ## ensure any variations with citations are displayed in browser tracks/ returned by default
    my $vdisplay_upt_sth  = $dbh->prepare(qq[ update variation set display =? where  variation_id =?  ]);
    my $vfdisplay_upt_sth = $dbh->prepare(qq[ update variation_feature set display =? where  variation_id =?  ]);

    my $varfeat_ext_sth   = $dbh->prepare(qq[ select variation_feature_id  
                                              from variation_feature
                                              where  variation_id =? ]);
 
    my $tvdisplay_upt_sth = $dbh->prepare(qq[ update transcript_variation 
                                              set display =? 
                                              where  variation_feature_id =?
                                           ]);

    my @var_objects ;

    if(  defined $var->[0]){ 
        @var_objects = @{$var} ;
    }
    elsif( defined $pub->{variants}->[0]){ 
        @var_objects = @{$pub->{variants}} ;
    }
    else{ 
        throw("No variants to link to PMID". $pub->{pmid} ); 
        return;
    }


    foreach my $var_obj ( @var_objects  ){

        if( ! $var_obj->isa("Bio::EnsEMBL::Variation::Variation")) {
            throw("Bio::EnsEMBL::Variation::Variation object expected");
        }

        ## avoid duplicates
        $citation_ext_sth->execute( $var_obj->dbID(),  $pub->{dbID} );
        my $count = $citation_ext_sth->fetchall_arrayref();
        next unless $count->[0]->[0] ==0;

        $citation_ins_sth->execute( $var_obj->dbID(), $pub->{dbID});

	## set cited variants to be displayable
	$vdisplay_upt_sth->execute( 1,  $var_obj->dbID());
        $vfdisplay_upt_sth->execute( 1,  $var_obj->dbID());
        ## don't join TV & VF in update
        $varfeat_ext_sth->execute($var_obj->dbID());
        my $vf_ids = $varfeat_ext_sth->fetchall_arrayref();
        foreach my $vf_id (@{$vf_ids}){
            $tvdisplay_upt_sth->execute( 1,  $vf_id->[0]);
        }
    }
}

sub update_ucsc_id {

    my $self    = shift;       
    my $pub     = shift;
    my $ucsc_id = shift;   

    throw("No UCSC external id defined") unless defined $ucsc_id;
    throw("No publication defined")      unless defined $pub->{dbID};

    $pub->{ucsc_id} = $ucsc_id;;

    my $dbh = $self->dbc->db_handle;
    
    my $publication_updt_sth = $dbh->prepare(qq[ update publication set ucsc_id = ? where publication_id =?  ]);   

    $publication_updt_sth->execute( $ucsc_id, $pub->{dbID});

}


1;
