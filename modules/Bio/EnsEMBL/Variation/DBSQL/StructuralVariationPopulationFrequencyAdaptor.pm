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


# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::StructuralVariationPopulationFrequencyAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::StructuralVariationPopulationFrequencyAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $sva   = $reg->get_adaptor("human","variation","structuralvariation");
  $svpf_adaptor = $reg->get_adaptor("human","variation","structuralvariationpopulationfrequency");
  
  # Get a StructuralVariation by its internal identifier
  $sv = $sva->fetch_by_dbID(145);

  # Get a StructuralVariation by its name
  $sv = $sva->fetch_by_name('esv3631253');
  
  # Get the StructuralVariationPopulationFrequency object from the StructuralVariation object
  my $svpfs = $svpf_adaptor->fetch_all_by_StructuralVariation($sv);
  
  foreach my $svpf (@$svpfs) {
    my $pop_name = $svpf->name;
    my $samples_count = 0;
    
    # Global frequency
    foreach my $SO_term (keys(%{$svpf->{samples_class}})) {
      $samples_count += scalar(keys %{$svpf->{samples_class}->{$SO_term}});
    }
	  print $pop_name.">> Global frequency: ".sprintf("%.4f",$svpf->frequency)." (Samples: $samples_count | Pop size: ".$svpf->size.")\n";
	  
	  # Allele class frequency
	  my $freqs_by_SO_term = $svpf->frequencies_by_class_SO_term;
	  foreach my $SO_term (keys(%$freqs_by_SO_term)) {
	    print "> $SO_term: ".sprintf("%.4f",$freqs_by_SO_term->{$SO_term})."\n";
	  }
  }
 
=head1 DESCRIPTION

This adaptor provides database connectivity for StructuralVariationPopulationFrequency objects.
This StructuralVariationPopulationFrequency object provides the population allele class frequencies 
for a given structural variant.


=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::StructuralVariationPopulationFrequencyAdaptor;

use Bio::EnsEMBL::Variation::StructuralVariationPopulationFrequency;
use Bio::EnsEMBL::Variation::StructuralVariation;
use Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use DBI qw(:sql_types);

our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor');


# method used by superclass to construct SQL
sub _tables { 
  my $self = shift;
  my @tables = (
    ['structural_variation', 'sv'],
    ['structural_variation_association', 'sva'],
    ['structural_variation_sample', 'svs'],
    ['structural_variation_feature', 'svf'],
    ['seq_region', 'seq'],
    ['sample_population', 'sp']
  );
  
  # If we are excluding failed_structural_variations, add that table
  push(@tables,['failed_structural_variation', 'fsv']) unless ($self->db->include_failed_variations());
  
  return @tables;
}

sub _columns {
  return qw( sv.structural_variation_id sp.population_id sv.class_attrib_id sp.sample_id svs.zygosity seq.name );
}

# Add a left join to the failed_structural_variation table
sub _left_join {
  my $self = shift;
  
  # If we are including failed structural variations, skip the left join
  return () if ($self->db->include_failed_variations());
  return (['failed_structural_variation', 'fsv.structural_variation_id=sva.structural_variation_id']);
}

sub _default_where_clause {
  my $self = shift;
  return 'sva.supporting_structural_variation_id=sv.structural_variation_id AND svs.sample_id=sp.sample_id AND svs.structural_variation_id=sva.supporting_structural_variation_id AND sv.structural_variation_id=svf.structural_variation_id AND svf.seq_region_id=seq.seq_region_id';
}

sub _objs_from_sth {
  my ($self, $sth) = @_;
   
  my %row;
  # Create the row hash using column names as keys
  $sth->bind_columns( \( @row{ @{$sth->{NAME_lc} } } ));

  while ($sth->fetch) {

      # we don't actually store the returned object because
      # the _obj_from_row method stores them in a temporary
      # hash _temp_objs in $self 
      $self->_obj_from_row(\%row);
  }

  # Get the created objects from the temporary hash
  my @objs = sort { $a->population->dbID <=> $b->population->dbID} values %{ $self->{_temp_objs} };
  delete $self->{_temp_objs};
 
  # Return the created objects 
  return \@objs;
}


sub _obj_from_row {

  my ($self, $row, $mapper, $dest_slice) = @_;

  my $aa  = $self->db->get_AttributeAdaptor;

  my $obj = $self->{_temp_objs}{$row->{population_id}}; 

  unless (defined($obj)) {
 
    my $pa = $self->db()->get_PopulationAdaptor();
    my $pop = $pa->fetch_by_dbID($row->{population_id});
    my $pop_name = $pop->name;
    my $pop_desc = $pop->description;
    my $pop_size = $pop->size;
 
    $obj = Bio::EnsEMBL::Variation::StructuralVariationPopulationFrequency->new(
         -adaptor                  => $self,
         -_structural_variation_id => $row->{structural_variation_id},
         -_population_id           => $row->{population_id},
         -name                     => $pop_name,
         -description              => $pop_desc,
         -size                     => $pop_size,
         -region_name              => $row->{name}
    );

    $self->{_temp_objs}{$row->{population_id}} = $obj;
  }

  my $class_SO_term = $aa->attrib_value_for_id($row->{class_attrib_id});

  ## Add sample IDs by SO terms (with zygosity)
  
  # Gets only the unique ontology accessions in a hash
  $obj->{samples_class_hash}{$class_SO_term}{$row->{sample_id}} = $row->{zygosity};
  
  # Then we transform the hash into an array
  foreach my $SO_term (keys(%{$obj->{samples_class_hash}})) {
    my $samples = $obj->{samples_class_hash}{$SO_term};
    $obj->{samples_class}{$SO_term} = $samples;
  }
}


=head2 fetch_all_by_StructuralVariation

  Arg [1]    : Bio::EnsEMBL:Variation::StructuralVariation $sv
  Example    : my @svpfs = @{$svpfa->fetch_all_by_StructuralVariation($sv)};
  Description: Retrieves all structural variation features for a given structural variation. Most
               structural variations should only hit the genome once and only a return
               a single structural variation feature.
  Returntype : reference to list Bio::EnsEMBL::Variation::StructuralVariationPopulationFrequency
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut


sub fetch_all_by_StructuralVariation {
  my $self = shift;
  my $var  = shift;

  if(!ref($var) || !$var->isa('Bio::EnsEMBL::Variation::StructuralVariation')) {
    throw('Bio::EnsEMBL::Variation::StructuralVariation arg expected');
  }

  if(!defined($var->dbID())) {
    throw("StructuralVariation arg must have defined dbID");
  }
  my $constraint = $self->_internal_exclude_failed_constraint("sva.structural_variation_id = ".$var->dbID());
  
  return $self->generic_fetch($constraint);
}


# Exclude the constraint for failed structural variant
sub _internal_exclude_failed_constraint {
  my $self = shift;
  my $constraint = shift;
  $constraint .= " AND " . $self->db->_exclude_failed_structural_variations_constraint();
  
  return $constraint;
}




1;
