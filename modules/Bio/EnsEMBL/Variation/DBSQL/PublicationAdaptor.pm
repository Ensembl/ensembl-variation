=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/legal/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

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

=head2 fetch_all_by_variant

  Arg [1]    : listref $list
  Example    : $publication = $publication_adaptor->fetch_all_by_variant( $var_object]);
  Description: Retrieves a listref of publication objects via a list of internal
               dbID identifiers
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
  return qw(p.publication_id p.title p.authors p.pmid p.pmcid );
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

  my ($publication_id,$title,$authors,$pmid,$pmcid);
  
  $sth->bind_columns(\$publication_id, \$title, \$authors, \$pmid, \$pmcid);
  
  while($sth->fetch()) {
		
		
      push @publication, Bio::EnsEMBL::Variation::Publication->new
	  (-dbID       => $publication_id,
	   -ADAPTOR    => $self,
	   -TITLE      => $title,
	   -AUTHORS    => $authors,
	   -PMID       => $pmid,
	   -PMCID      => $pmcid,
	  );
  }

  return \@publication;
}


sub store{

    my $self = shift;
    my $pub  = shift;
    
    my $dbh = $self->dbc->db_handle;
    
    my $pub_ins_sth = $dbh->prepare(qq[ insert into publication( title, authors, pmid, pmcid) values ( ?,?,?,?) ]);   
    
    $pub_ins_sth->execute( $pub->{title},
			   $pub->{authors},
			   $pub->{pmid}  || undef,
			   $pub->{pmcid} || undef );


    $pub->{dbID} = $dbh->last_insert_id(undef, undef, 'publication', 'publication_id');

    ## check for cited variants to link 
    if( $pub->{variants}->[0] ){
	
	$self->update_variant_citation($pub);
    }
}

sub update_variant_citation {

    my $self = shift;       
    my $pub  = shift;
    my $var  = shift;   

    my $dbh = $self->dbc->db_handle;

    my $citation_ext_sth = $dbh->prepare(qq[ select count(*) from variation_citation where variation_id =? and publication_id =?  ]);   
    
    my $citation_ins_sth = $dbh->prepare(qq[ insert into variation_citation( variation_id, publication_id) values ( ?,?) ]);   
    
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
	
    }    
}


1;
