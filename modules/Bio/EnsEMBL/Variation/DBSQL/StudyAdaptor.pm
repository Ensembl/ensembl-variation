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

#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::StudyAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::StudyAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $sta = $reg->get_adaptor("human","variation","study");

  # fetch a study by its name
  $study = $sta->fetch_by_name('estd1'); 

  # fetch all study for a source
  $sta = $reg->get_adaptor("human","variation","study");
  $st = $sta->fetch_all_by_source('NHGRI_GWAS_catalog');
  foreach $study (@{$sta->fetch_all_by_source('NHGRI_GWAS_catalog')}){
  	print $study->dbID, " - ", $study->external_reference ,"\n"; 
  }
  

=head1 DESCRIPTION

This adaptor provides database connectivity for Study objects.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::StudyAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Study;

use base qw{Bio::EnsEMBL::DBSQL::BaseAdaptor};

my %cache;

=head2 fetch_by_name

  Arg [1]    : string $name
  Example    : $study = $study_adaptor->fetch_by_name('estd1');
  Description: Retrieves a study object via its name
  Returntype : Bio::EnsEMBL::Variation::Study
  Exceptions : throw if name argument is not defined
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_name {
  my $self = shift;
  my $name = shift;

  throw('name argument expected') if(!defined($name));

  my $result = $self->generic_fetch("st.name='$name'");

  return ($result ? $result->[0] : undef);
}


=head2 fetch_by_dbID

  Arg [1]    : int $dbID
  Example    : $study = $study_adaptor->fetch_by_dbID(254);
  Description: Retrieves a Study object via its internal identifier.
               If no such study exists undef is returned.
  Returntype : Bio::EnsEMBL::Variation::Study
  Exceptions : throw if dbID arg is not defined
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  throw('dbID argument expected') if(!defined($dbID));

  if (exists($cache{$dbID})) {
    return $cache{$dbID};
  }
  my $result = $self->generic_fetch("st.study_id=$dbID");

  if ($result) {
    $cache{$dbID} = $result->[0];
  }

  return ($result ? $result->[0] : undef);
}

	
=head2 fetch_all_by_dbID_list

  Arg [1]    : listref $list
  Example    : $study = $study_adaptor->fetch_all_by_dbID_list([907,1132]);
  Description: Retrieves a listref of study objects via a list of internal
               dbID identifiers
  Returntype : listref of Bio::EnsEMBL::Variation::Study objects
  Exceptions : throw if list argument is not defined
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_dbID_list {
  my $self = shift;
  my $list = shift;

  if(!defined($list) || ref($list) ne 'ARRAY') {
    throw("list reference argument is required");
  }
  
  return undef if (scalar(@$list)==0);
  
  my $id_str = (@$list > 1)  ? " IN (".join(',',@$list).")"   :   ' = \''.$list->[0].'\'';
	
	my $result = $self->generic_fetch("st.study_id $id_str");

  return ($result ? $result : undef);
}


=head2 fetch_all_by_source

  Arg [1]     : string $source_name
  Example     : my $study = $study_adaptor->fetch_by_name('EGAS00000000001');
  Description : Retrieves all Study objects associated with a source.
	Returntype : listref of Bio::EnsEMBL::Variation::Study
  Exceptions : thrown if source_name not provided
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_source{
	my $self = shift;
	my $source_name = shift;

	throw('source_name argument expected') if(!defined($source_name));
		
	my $result = $self->generic_fetch("s.name='$source_name'");

	return ($result ? $result : undef);
}


sub _fetch_all_associate_study_id {
		my $self = shift;
    my $study_id = shift;

	  my $a_study;
		my @study_list;
		
		my $sth = $self->prepare(qq{(SELECT DISTINCT study1_id FROM associate_study WHERE study2_id=?)
																 UNION
																(SELECT DISTINCT study2_id FROM associate_study WHERE study1_id=?)});
		$sth->bind_param(1,$study_id,SQL_INTEGER);
		$sth->bind_param(2,$study_id,SQL_INTEGER);
    $sth->execute();
		$sth->bind_columns(\$a_study);

  	while($sth->fetch()) {
			push(@study_list,$a_study);
		}
		return \@study_list;
}

=head2 fetch_all_by_external_reference

  Arg [1]    : string $external_reference
  Example    : my $study = $study_adaptor->fetch_by_external_reference('PMID:25917933');
  Description: Retrieves the Study object associated with an external reference.
  Returntype : listref of Bio::EnsEMBL::Variation::Study
  Exceptions : thrown if external reference not provided
  Caller     : general
  Status     : Experimental

=cut
sub fetch_all_by_external_reference{

  my $self    = shift;
  my $ext_ref = shift;

  throw('external reference argument expected') if !defined $ext_ref ;

  my $result = $self->generic_fetch("st.external_reference='$ext_ref'");

  ## Could be multiple studies if 2 groups extracted data from same publication
  return ($result ? $result : undef);
}

sub _columns {
  return qw(st.study_id st.name st.description st.url st.external_reference st.study_type st.source_id);
}

sub _tables { return (['study', 'st'],['source', 's']); }

sub _default_where_clause {
  my $self = shift;
  return 'st.source_id = s.source_id';
}

#
# private method, creates study objects from an executed statement handle
# ordering of columns must be consistant
#
sub _objs_from_sth {
  my $self = shift;
  my $sth  = shift;

  my @study;

  my ($study_id,$study_name,$study_description,$study_url,$external_reference,$study_type,$source_id,$associate);

  $sth->bind_columns(\$study_id, \$study_name, \$study_description, \$study_url, 
	                   \$external_reference, \$study_type, \$source_id);

  while($sth->fetch()) {
		
		$associate = $self->_fetch_all_associate_study_id($study_id);
		
    push @study, Bio::EnsEMBL::Variation::Study->new
      (-dbID => $study_id,
       -ADAPTOR => $self,
       -NAME => $study_name,
       -DESCRIPTION => $study_description,
			 -URL => $study_url,
			 -EXTERNAL_REFERENCE => $external_reference,
			 -TYPE => $study_type,
			 -_SOURCE_ID => $source_id,
			 -ASSOCIATE => $associate);
  }

  return \@study;
}

sub store {

  my ($self, $study ) = @_;
     
  my $dbh = $self->dbc->db_handle;
     
  my $sth = $dbh->prepare(q{
         INSERT INTO study (
             source_id,
             name,
             description,
             url,
             external_reference,
             study_type
         ) VALUES (?,?,?,?,?,?)
     });
     
  $sth->execute(
               $study->source->dbID(),
               $study->name,
               $study->description || undef,
               $study->url || undef,
               $study->external_reference || undef,
               $study->type || undef
              );
     
  $sth->finish;
     
# get dbID
  my $dbID = $dbh->last_insert_id(undef, undef, 'study', 'study_id');

  $study->{dbID}    = $dbID;

  $study->{adaptor} = $self;

}


1;
