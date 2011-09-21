package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::Constants;

use base qw(Exporter);

our @EXPORT_OK = qw(FULL UPDATE NONE);

use constant FULL   => 1;
use constant UPDATE => 2;
use constant NONE   => 3;

1;

