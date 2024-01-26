# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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




=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

# Read a variation config module and generate a perl module
# defining various constants used by the API which should
# be stored in Bio/EnsEMBL/Variation/Utils/Constants.pm and
# then checked into the variation CVS repository

use strict;
use warnings;

use Getopt::Long;

use Data::Dumper;

$Data::Dumper::Terse = 1;
$Data::Dumper::Sortkeys = 1;

my $config;
my $help;

GetOptions(
    "config=s"  => \$config,
    "help|h"    => \$help,
);

if ($help) {
    print "Usage: $0 --config <module> --help > Constants.pm\n";
    exit(0);
}

# pull in our config module

$config ||= 'Bio::EnsEMBL::Variation::Utils::Config';

eval qq{require $config};

die "Failed to require config module '$config':\n$@" if $@;

# and import the variables we need

our @ATTRIB_TYPES;
our @VARIATION_CLASSES;
our @OVERLAP_CONSEQUENCES; 
our @FEATURE_TYPES;
our $OVERLAP_CONSEQUENCE_CLASS;
our %SO_ACC_MAPPER;

eval {
    $config->import(qw(
        @ATTRIB_TYPES
        @VARIATION_CLASSES
        @OVERLAP_CONSEQUENCES 
        @FEATURE_TYPES
        %SO_ACC_MAPPER
        $OVERLAP_CONSEQUENCE_CLASS
    ));
};

die "Failed to import required data structures from config module '$config':\n$@" if $@;

my $hdr = q{package Bio::EnsEMBL::Variation::Utils::Constants;

#####################################################################
# NB: THIS FILE HAS BEEN AUTOMATICALLY GENERATED, EDIT WITH CAUTION #
#####################################################################

use strict;
use warnings;

use base qw(Exporter);

};

my $code;

my $to_export;

for my $attrib_type (@ATTRIB_TYPES) {
   
    my $type_code = delete $attrib_type->{code} or die "code required for attrib_type";
    
    my $const = "ATTRIB_TYPE_".uc($type_code);

    next if $to_export->{attrib_types}->{$const};
    
    $to_export->{attrib_types}->{$const} = 1;
   
    $code .= "use constant $const => '$type_code';\n";
}

$code .= "\nuse constant ATTRIB_TRAIT => 'trait';\n\n";

my $class_code = "our %VARIATION_CLASSES = (\n";

for my $class_set (@VARIATION_CLASSES) {

    my $term = $class_set->{SO_term};

    my $const = "SO_TERM_".uc($term);

    next if $to_export->{SO_class_terms}->{$const};

    $to_export->{SO_class_terms}->{$const} = 1;

    $code .= "use constant $const => '$term';\n";

    delete $class_set->{SO_term};

    $class_set->{display_term} ||= $term;

    my $display_term = $class_set->{display_term};
    if ($display_term =~ /^[A-Z][a-z]+/) {
      $display_term = lcfirst($display_term);
    }
    $class_set->{somatic_display_term} ||= 'somatic '.$display_term;
    
    $class_code .= "'$term' => ".Dumper($class_set).",\n";
}

$class_code .= ");\n";

my $cons_code = "our \%OVERLAP_CONSEQUENCES = (\n";

my $default_consequence_code;

for my $cons_set (@OVERLAP_CONSEQUENCES) {

    my $term = $cons_set->{SO_term};

    # set a default rank
    $cons_set->{rank} ||= '999';
    
    my $const = "SO_TERM_".uc($term);

    $code .= "use constant $const => '$term';\n" unless $to_export->{SO_consequence_terms}->{$const} || $to_export->{SO_class_terms}->{$const};

    $to_export->{SO_consequence_terms}->{$const} = 1;
    
    if ($cons_set->{is_default}) {
        die "Cannot have more than one default OverlapConsequence object!" if $default_consequence_code;
        $default_consequence_code = "our \$DEFAULT_OVERLAP_CONSEQUENCE = $OVERLAP_CONSEQUENCE_CLASS->new_fast(".Dumper($cons_set).");\n";
        $cons_code .= "'$term' => \$DEFAULT_OVERLAP_CONSEQUENCE,\n";
    }
    else {
        $cons_code .= "'$term' => $OVERLAP_CONSEQUENCE_CLASS->new_fast(".Dumper($cons_set)."),\n";
    }
}

$cons_code = $default_consequence_code. "\n". "\n" . $cons_code .");\n";



my $feat_so_acc = "our \$SO_ACC_MAPPER = " . Dumper(\%SO_ACC_MAPPER) . ";\n";


# avoid any duplicate exports by putting all constants to export into a single hash

my $all_to_export;

for my $type (keys %$to_export) {
    for my $export (keys %{ $to_export->{$type} }) {
        $all_to_export->{$export} = 1;
    }
}

my @header = qw(
  %OVERLAP_CONSEQUENCES
  %VARIATION_CLASSES
  $DEFAULT_OVERLAP_CONSEQUENCE
  $SO_ACC_MAPPER
  ATTRIB_TRAIT
);

$hdr .= 'our @EXPORT_OK = qw(' .
 "\n  " . (join "\n  ", @header) .
 "\n  " . (join "\n  ", sort keys %$all_to_export) . "\n);";

$hdr .= "\n\n";

$hdr .= 'our %EXPORT_TAGS = (';

for my $key (sort keys %$to_export) {
    $hdr .= "\n  $key => [qw(\n    " .
            (join "\n    ", sort keys %{ $to_export->{$key} }) .
            "\n  )],";
}

$hdr .= "\n);\n\n";

$hdr .= "use $OVERLAP_CONSEQUENCE_CLASS;\n";

print $hdr, "\n", $code, "\n", $class_code, "\n", $cons_code, "\n", $feat_so_acc, "\n1;\n";
