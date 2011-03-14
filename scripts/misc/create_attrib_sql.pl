use Bio::EnsEMBL::Variation::Utils::Config qw(
    @ATTRIB_TYPES
    @VARIATION_CLASSES 
    @OVERLAP_CONSEQUENCES 
    @FEATURE_TYPES
);

my %attrib_ids;
my $attrib_type_ids;

my $last_attrib_type_id = 1;
my $last_attrib_id      = 1;
my $last_attrib_set_id  = 1;

my $attrib_type_fmt = 
    q{INSERT INTO attrib_type (attrib_type_id, code, name, description) VALUES (%d, %s, %s, %s);};
my $attrib_fmt = 
    q{INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (%d, %d, '%s');};
my $set_fmt = 
    q{INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (%d, %d);};

my $SQL;

# first define the attrib type entries

for my $attrib_type (@ATTRIB_TYPES) {
    
    my $attrib_type_id = $last_attrib_type_id++;

    my $code        = delete $attrib_type->{code} or die "code required for attrib_type";
    my $name        = delete $attrib_type->{name};
    my $description = delete $attrib_type->{description};

    die "Unexpected entries in attrib_type definition: ".(join ',', keys %$attrib_type)
        if keys %$attrib_type;

    $SQL .= sprintf($attrib_type_fmt, 
        $attrib_type_id, 
        "'$code'",
        ($name ? "'$name'" : "''"),
        ($description ? "'$description'" : 'NULL'),
    )."\n";

    $attrib_type_ids->{$code} = $attrib_type_id;
}

# now define the SO term mappings

my $seen_SO_mappings;

$SQL .= "\n";

for my $class_set (@VARIATION_CLASSES, @OVERLAP_CONSEQUENCES, @FEATURE_TYPES) {

    my $set_sql;

    my @SO_mappings;

    for my $type (qw(SO_accession SO_term display_term)) {
    
        my $attrib_type_id = $attrib_type_ids->{$type} or die "Unrecognised attrib_type: $type";

        my $value = $class_set->{$type};
        
        next unless $value;
            
        push @SO_mappings, $value;

        my $attrib_id_key = $type.'_'.$value;
        
        my $attrib_id = $attrib_ids{$attrib_id_key};
        
        unless ($attrib_id) {
            $attrib_id = $last_attrib_id++;
            $SQL .= sprintf($attrib_fmt, $attrib_id, $attrib_type_id, $value)."\n";
            $attrib_ids{$attrib_id_key} = $attrib_id;
        }

        $set_sql .= (sprintf $set_fmt, $last_attrib_set_id, $attrib_id)."\n";
    }

    $SQL .= $set_sql."\n" unless $seen_SO_mappings->{join '_', @SO_mappings}++;

    $last_attrib_set_id++;
}

print $SQL;
