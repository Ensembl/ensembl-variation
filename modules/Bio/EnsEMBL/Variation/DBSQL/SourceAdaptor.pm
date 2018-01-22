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
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::SourceAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::SourceAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $sta = $reg->get_adaptor("human","variation","source");

  # fetch a source by its name
  $source = $sta->fetch_by_name('dbSNP'); 
  

=head1 DESCRIPTION

This adaptor provides database connectivity for Source objects.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::SourceAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Source;

use DBI qw(:sql_types);

use base qw{Bio::EnsEMBL::DBSQL::BaseAdaptor};


=head2 fetch_by_name

  Arg [1]    : string $name
  Example    : $source = $source_adaptor->fetch_by_name('dbSNP');
  Description: Retrieves a source object via its name
  Returntype : Bio::EnsEMBL::Variation::Source
  Exceptions : throw if name argument is not defined
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_name {
  my $self = shift;
  my $name = shift;

  throw('name argument expected') if(!defined($name));

  my $result = $self->generic_fetch("s.name='$name'");

  return ($result ? $result->[0] : undef);
}


=head2 fetch_by_dbID

  Arg [1]    : int $dbID
  Example    : $source = $source_adaptor->fetch_by_dbID(1);
  Description: Retrieves a Source object via its internal identifier.
               If no such source exists undef is returned.
  Returntype : Bio::EnsEMBL::Variation::Source
  Exceptions : throw if dbID arg is not defined
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  throw('dbID argument expected') if(!defined($dbID));

  my $result = $self->generic_fetch("s.source_id=$dbID");

  return ($result ? $result->[0] : undef);
}

	
=head2 fetch_all_by_dbID_list

  Arg [1]    : listref $list
  Example    : $source = $source_adaptor->fetch_all_by_dbID_list([1,2]);
  Description: Retrieves a listref of source objects via a list of internal
               dbID identifiers
  Returntype : listref of Bio::EnsEMBL::Variation::Source objects
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
	
	my $result = $self->generic_fetch("s.source_id $id_str");

  return ($result ? $result : undef);
}


sub _columns {
  return qw(s.source_id s.name s.version s.description s.url s.type s.somatic_status s.data_types);
}

sub _tables { return (['source', 's']); }


#
# private method, creates source objects from an executed statement handle
# ordering of columns must be consistant
#
sub _objs_from_sth {
  my $self = shift;
  my $sth  = shift;

  my @source;

  my ($source_id,$source_name,$source_version,$source_description,$source_url,$source_type,$source_somatic_status,$source_data_types);

  $sth->bind_columns(\$source_id, \$source_name, \$source_version, \$source_description, \$source_url, \$source_type, \$source_somatic_status, \$source_data_types);

  while($sth->fetch()) {
		
    my @data_types = (defined($source_data_types)) ? split(/,/,$source_data_types) : undef;
		
    push @source, Bio::EnsEMBL::Variation::Source->new
      (-dbID => $source_id,
       -ADAPTOR => $self,
       -NAME => $source_name,
       -VERSION => $source_version,
       -DESCRIPTION => $source_description,
			 -URL => $source_url,
			 -TYPE => $source_type,
			 -SOMATIC_STATUS => $source_somatic_status,
			 -DATA_TYPES => \@data_types);
  }

  return \@source;
}


=head2 get_source_version

  Arg[1]      : string $name
  Example     : $version = $sa->get_source_version('dbSNP');
  Description : Retrieves from the database the version for the source given as an argument
  ReturnType  : int
  Exceptions  : none
  Caller      : general
  Status      : Stable

=cut

sub get_source_version{
    my $self = shift;
    my $name = shift;
    my $version;
    my $sth = $self->prepare(qq{SELECT version from source where name = ? });

    $sth->bind_param(1,$name,SQL_VARCHAR);

    $sth->execute();
    $sth->bind_columns(\$version);
    $sth->fetch();
    $sth->finish();

    return $version;
}


sub store {
  my ($self, $source) = @_;
    
	my $dbh = $self->dbc->db_handle;
    
  my $sth = $dbh->prepare(q{
        INSERT INTO source (
            name,
            version,
            description,
            url,
            type,
            somatic_status,
            data_types
        ) VALUES (?,?,?,?,?,?,?)
    });
    
    $sth->execute(
        $source->name,
        $source->version || undef,
        $source->description || undef,
        $source->url || undef,
        $source->type || undef,
        $source->somatic_status || 'germline',
        (join ",", @{$source->get_all_data_types}) || undef
    );
    
    $sth->finish;
    
    # get dbID
	  my $dbID = $dbh->last_insert_id(undef, undef, 'source', 'source_id');
    $source->{dbID}    = $dbID;
    $source->{adaptor} = $self;
}

=head2 update_version

  Arg[1]      : source object
  Example     : $sa->update_version( $source_object);
  Description : Update the version for the source
  ReturnType  : none
  Exceptions  : none
  Caller      : internal pipelines
  Status      : Experimental

=cut
sub update_version{

  my ($self, $source) = @_;

  ## don't over-write to null
  return unless $source->version();

  my $dbh = $self->dbc->db_handle;
  my $sth = $dbh->prepare(q[ update source set version = ? where source_id =?]);
  
  $sth->execute( $source->version(), $source->dbID() );

}

1;
