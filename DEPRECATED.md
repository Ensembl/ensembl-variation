Ensembl Variation Deprecated Methods
===================

This file contains the list of methods deprecated in the Ensembl Variation API. A method is deprecated when it is not functional anymore (schema/data change) or has been replaced by a better one. Backwards compatibility is provided whenever possible. When a method is deprecated, a deprecation warning is thrown whenever the method is used. The warning also contains instructions on replacing the deprecated method and when it will be removed. A year after deprecation (4 Ensembl releases), the method is removed from the API.

### Removed in Ensembl Release 87 (deprecated from Ensembl Release 83) ###

#### Deprecated methods related to the validation_status data:####
- Bio::EnsEMBL::Variation::**Variation**::*get_all_validation_states()*
- Bio::EnsEMBL::Variation::**Variation**::*add_validation_state()*
- Bio::EnsEMBL::Variation::**VariationFeature**::*get_all_validation_states()*
- Bio::EnsEMBL::Variation::**VariationFeature**::*add_validation_state()*
- Bio::EnsEMBL::Variation::Utils::**Sequence**::*get_all_validation_states()*
- Bio::EnsEMBL::Variation::Utils::**Sequence**::*get_validation_code()*
- Bio::EnsEMBL::Variation::Utils::**Sequence**::*add_validation_state()*

#### Deprecated method source_object() in various objects:####
- Bio::EnsEMBL::Variation::**Variation**::*source_object()*
- Bio::EnsEMBL::Variation::**VariationFeature**::*source_object()*
- Bio::EnsEMBL::Variation::**PhenotypeFeature**::*source_object()*
- Bio::EnsEMBL::Variation::**BaseStructuralVariation**::*source_object()*
- Bio::EnsEMBL::Variation::**StructuralVariationFeature**::*source_object()*
- Bio::EnsEMBL::Variation::**Study**::*source_object()*


