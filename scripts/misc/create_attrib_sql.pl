use Bio::EnsEMBL::Variation::Utils::Config qw(
    @ATTRIB_TYPES
    @ATTRIB_SETS
    %ATTRIBS
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

# second, take the entries from the ATTRIBS and add them as single-element hashes to the @ATTRIB_SETS array
while (my ($type,$values) = each(%ATTRIBS)) {
    
    map {push(@ATTRIB_SETS,{$type => $_})} @{$values};
} 

# third, loop over the ATTRIB_SETS array and add attribs and assign them to attrib_sets as necessary
for my $set (@ATTRIB_SETS) {
    
    #ÊKeep the attrib_ids
    my @attr_ids;
    
    # Iterate over the type => value entries in the set
    while (my ($type,$value) = each(%{$set})) {
        
        # Lookup the attrib_type
        my $attrib_type_id = $attrib_type_ids->{$type} or next;
        
        # insert a new attrib if we haven't seen it before
        my $attrib_id = $attrib_ids{$attrib_type_ids . "_" . $value};
        
        unless (defined($attrib_id)) {
            
            $attrib_id = $last_attrib_id++;
            $SQL .= sprintf($attrib_fmt, $attrib_id, $attrib_type_id, $value)."\n"; 
            $attrib_ids{$attrib_type_ids . "_" . $value} = $attrib_id;
            
        }
        
        push(@attr_ids,$attrib_id);   
    }
    
    # If the set had more than one attribute, group them into a set
    if (scalar(@attr_ids) > 1) {
        
        my $attrib_set_id = $last_attrib_set_id++;
        map {$SQL .= sprintf($set_fmt, $attrib_set_id, $_)."\n"} @attr_ids;
        
    }
}

print $SQL . "\n";
