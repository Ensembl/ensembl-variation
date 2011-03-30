use strict;
use warnings;

use Data::Dumper;

$Data::Dumper::Terse = 1;

use Bio::EnsEMBL::Variation::Utils::Config qw(
    @ATTRIB_TYPES
    @VARIATION_CLASSES 
    @OVERLAP_CONSEQUENCES 
    @FEATURE_TYPES
    $OVERLAP_CONSEQUENCE_CLASS
);


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

$code .= "\n";

for my $class_set (@VARIATION_CLASSES) {

    my $term = $class_set->{SO_term};

    my $const = "SO_TERM_".uc($term);

    next if $to_export->{SO_class_terms}->{$const};

    $to_export->{SO_class_terms}->{$const} = 1;

    $code .= "use constant $const => '$term';\n";
}

$code .= "\n";

my $cons_code = "our \@OVERLAP_CONSEQUENCES = (\n";

my $default_consequence_code;

for my $cons_set (@OVERLAP_CONSEQUENCES) {

    my $term = $cons_set->{SO_term};
    
    my $const = "SO_TERM_".uc($term);

    $code .= "use constant $const => '$term';\n" unless $to_export->{SO_consequence_terms}->{$const};

    $to_export->{SO_consequence_terms}->{$const} = 1;
    
    if ($cons_set->{is_default}) {
        die "Cannot have more than one default OverlapConsequence object!" if $default_consequence_code;
        $default_consequence_code = "our \$DEFAULT_OVERLAP_CONSEQUENCE = $OVERLAP_CONSEQUENCE_CLASS->new_fast(".Dumper($cons_set).");\n";
        $cons_code .= "\$DEFAULT_OVERLAP_CONSEQUENCE,\n";
    }
    else {
        $cons_code .= "$OVERLAP_CONSEQUENCE_CLASS->new_fast(".Dumper($cons_set)."),\n";
    }
}

$cons_code = $default_consequence_code. "\n". $cons_code .");\n";

$hdr .= 'our @EXPORT_OK = qw(@OVERLAP_CONSEQUENCES $DEFAULT_OVERLAP_CONSEQUENCE '.(join ' ', map { keys %{ $to_export->{$_} } } keys %$to_export).');';

$hdr .= "\n\n";

$hdr .= 'our %EXPORT_TAGS = ( ';

for my $key (keys %$to_export) {
    $hdr .= "$key => [qw(".(join ' ', keys %{ $to_export->{$key} }).")], "
}

$hdr .= " );\n\n";

$hdr .= "use $OVERLAP_CONSEQUENCE_CLASS;\n";

print $hdr, "\n", $code, "\n", $cons_code, "\n1;\n";
