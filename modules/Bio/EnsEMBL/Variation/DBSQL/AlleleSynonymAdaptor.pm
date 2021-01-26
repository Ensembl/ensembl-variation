=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Variation::DBSQL::AlleleSynonymAdaptor

=head1 SYNOPSIS

  $reg = 'Bio::EnsEMBL::Registry';
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org', -user => 'anonymous');
  
  $asa = $reg->get_adaptor('homo_sapiens', 'variation', 'allelesynonym');
  $va = $reg->get_adaptor('homo_sapiens', 'variation', 'variation');

  # Get an AlleleSynonym by its internal identfier
  $allele_synonym = $asa->fetch_by_dbID(3674251);

  # Fetch all allele synonyms by name
  $allele_synonyms = $asa->fetch_all_by_name('CA75242691');

  # Fetch all allele synonyms for a particular variation
  $v = $va->fetch_by_name('rs916798809');
  foreach my $as (@{$asa->fetch_all_by_Variation($v)}) {
    print join("\t", $v->name, $as->name(), $as->hgvs_genomic()), "\n";
  }
  

=head1 DESCRIPTION

This adaptor provides database connectivity for AlleleSynonym objects.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::AlleleSynonymAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::AlleleSynonym;

use base qw{Bio::EnsEMBL::DBSQL::BaseAdaptor};

=head2 fetch_by_dbID

  Arg [1]    : int $dbID
  Example    : $allele_synonym = $as_adaptor->fetch_by_dbID(3674251);
  Description: Retrieves an AlleleSynonym object via its internal identifier.
               If no such allele synonym exists undef is returned.
  Returntype : Bio::EnsEMBL::Variation::AlleleSynonym
  Exceptions : throw if dbID arg is not defined
  Caller     : general
  Status     : At Risk

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  throw('dbID argument expected') if(!defined($dbID));

  my $result = $self->generic_fetch("vas.allele_synonym_id=$dbID");
 
  return ($result ? $result->[0] : undef);
}

=head2 fetch_all_by_Variation

  Arg [1]    : Bio::EnsEMBL::Variation::Variation $var
  Example    : my @allele_synonyms = @{$as_adaptor->fetch_all_by_Variation($v)
  Description: Retrieves all AlleleSynonyms for a given variation object.
  Returntype : arrayref of Bio::EnsEMBL::Variation::AlleleSynonym objects
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
    
    my @vas;
    my $allele_synonym_id;

    my $sth = $self->prepare(qq{SELECT allele_synonym_id from allele_synonym where variation_id = ?  });
    $sth->execute($var_obj->dbID);
    $sth->bind_columns(\$allele_synonym_id);
    while ($sth->fetch()){
        push @vas, $self->fetch_by_dbID($allele_synonym_id)
    }
    $sth->finish;
    
    return \@vas;
}

=head2 fetch_all_by_name

  Arg [1]    : string $name 
               The name of the allele_synonym to retrieve.
  Example    : my @allele_synonyms = @{$as_adaptor->fetch_all_by_name('CA75242691')};
  Description: Retrieves all allele synonyms with a specified name.  Allele synonyms
               names may be non-unique which is why this method returns a
               reference to a list.
  Returntype : arrayref of Bio::EnsEMBL::Variation::AlleleSynonym objects
  Exceptions : throw if no argument passed
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_name {
  my $self = shift;
  my $name = shift;

  throw('name argument expected') unless ($name);

  my $sth = $self->prepare(q{
    SELECT allele_synonym_id, variation_id, hgvs_genomic, name
    FROM   allele_synonym 
    WHERE  name = ?
  });

  $sth->bind_param(1, $name, SQL_VARCHAR);
  $sth->execute();
  my $vas =  $self->_objs_from_sth($sth);
  $sth->finish();
  return $vas;
}


sub _columns {
  return qw(vas.allele_synonym_id  vas.variation_id vas.hgvs_genomic vas.name);
}

sub _tables { return (['allele_synonym', 'vas']); }

sub _default_where_clause {
  my $self = shift;
  return ;
}

#
# private method, creates allele synonym objects from an executed statement handle
# ordering of columns must be consistent
#
sub _objs_from_sth {
  my $self = shift;
  my $sth  = shift;

  my @allele_synonyms;

  my ($allele_synonym_id, $variation_id, $hgvs_genomic, $name);
  my @results;
  
  $sth->bind_columns(\$allele_synonym_id, \$variation_id, \$hgvs_genomic, \$name);
  
  while($sth->fetch()) {
    push @results,
         Bio::EnsEMBL::Variation::AlleleSynonym->new
            (-dbID          => $allele_synonym_id,
             -ADAPTOR       => $self,
             -_VARIATION_ID => $variation_id,
             -HGVS_GENOMIC  => $hgvs_genomic,
             -NAME          => $name
            );
  }
  return \@results; 
}

1;
