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

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Iterator;
use Bio::EnsEMBL::Variation::Utils::EnsEMBL2GFF3;
use FileHandle;
package BaseGVFDumper;

sub new {
	my $caller = shift;
	my $class = ref($caller) || $caller;
	my %params = @_;
	my $self = {};
	bless $self, $class;
	foreach my $attrib (keys %params) {
		$self->{$attrib} = $params{$attrib};
	}
	return $self;
}

my $vfa;
my $svfa;
my $slice_adaptor;
my $slices;

my $working_slice;
my $sub_slice_start;
my $sub_slice_end;
my $slice_end;

my $chunk_size;

my $set;

my $id_count = 0;

my $failed_ids;
my $failed_descs;
my $structural_failed_ids;

my $fh;

sub init {
	my $self = shift;

    $fh = FileHandle->new($self->{'output'}, "w");
    
    if ($self->{'include_failed'} || $self->{'just_failed'}) {
        $self->include_failed_variations();
    }
    if ($self->{'include_svs'} || $self->{'just_failed'}) {
        $svfa = $self->{'vdba'}->get_StructuralVariationFeatureAdaptor;
    }
    if ($self->{'set_name'}) {
        $self->init_set_name();
    }

    $vfa           = $self->{'vdba'}->get_VariationFeatureAdaptor;
    $slice_adaptor = $self->{'cdba'}->get_SliceAdaptor;
    $chunk_size    = $self->{'chunk_size'};

    if (@{$self->{'seq_region_names'}}) {
        my @seq_region_names = @{$self->{'seq_region_names'}};
        for my $name (@seq_region_names) {
            my $slice = $slice_adaptor->fetch_by_region('toplevel', $name)
                or die "Failed to get a slice for seq region name: $name\n";
            push @$slices, $slice;
        }
    } else {
        $slices = $slice_adaptor->fetch_all('toplevel');
        die "Didn't find any toplevel slices.\n" unless @$slices;
    }
   
    # print header before slices are used

	$working_slice = shift @$slices;
	$sub_slice_start = $working_slice->start; 
	$sub_slice_end = $sub_slice_start + $chunk_size - 1;
	$slice_end = $working_slice->end;
}

sub include_failed_variations {
    my $self = shift;
    # include failed variations, but flag them with an attribute
    my $vdba = $self->{'vdba'};
    $vdba->include_failed_variations(1);
    $self->{'vdba'} = $vdba;
    # joining to variation and failed variation is too slow, so
    # we cache the entire failed_variation table in memory, just 
    # storing the failed_description_id to save memory. We also 
    # store the failed_description table in memory to map these 
    # ids to the description
    my $dbh = $self->{'vdba'}->dbc->db_handle;
    my $sth = $dbh->prepare(qq{
        SELECT  failed_description_id, description
        FROM    failed_description
    });
    $sth->execute;
    while (my ($id, $desc) = $sth->fetchrow_array) {
        $failed_descs->{$id} = $desc;
    }

    $sth = $dbh->prepare(qq{
        SELECT  variation_id, failed_description_id
        FROM    failed_variation
    });
    $sth->execute;
    my $v_id;
    my $desc_id;
    $sth->bind_columns(\$v_id, \$desc_id);
    while ($sth->fetch) {
        $failed_ids->{$v_id} = $desc_id;
    }

    $sth = $dbh->prepare(qq{
        SELECT  structural_variation_id, failed_description_id
        FROM    failed_structural_variation
    });
    $sth->execute;
    $sth->bind_columns(\$v_id, \$desc_id);
    while ($sth->fetch) {
        $structural_failed_ids->{$v_id} = $desc_id;
    }
}

sub init_set_name {
    my $self = shift;
    my $vsa = $self->{'vdba'}->get_VariationSetAdaptor;
    $set = $vsa->fetch_by_name($self->{'set_name'})
        or die "Didn't find set $self->{'set_name'}";
}

sub has_next_slice {
	return ($sub_slice_end < $slice_end || (scalar @$slices > 0));	
}

sub next_slice {
    my $self = shift;
    return unless($self->has_next_slice);
	if ($sub_slice_start > $slice_end) {
		$working_slice = shift @$slices;
		$sub_slice_start = $working_slice->start; 
		$slice_end = $working_slice->end;
	}
	$sub_slice_end = $sub_slice_start + $chunk_size - 1;
	$sub_slice_end = $working_slice->end if ($sub_slice_end > $working_slice->end);
	my $sub_slice = $working_slice->sub_Slice($sub_slice_start, $sub_slice_end);
	$sub_slice_start = $sub_slice_end + 1;
	return $sub_slice;
}

sub vfs_on_slice {
	my $self = shift;
	my $slice = shift;
	my $vfs;
	if ($self->{'use_iterator'}) { 
        $vfs = $vfa->fetch_Iterator_by_Slice($slice);
        return $vfs;
    }
    elsif ($self->{'somatic'}) {
		$vfs = $vfa->fetch_all_somatic_by_Slice_constraint_with_TranscriptVariations($slice);
    }
    elsif ($self->{'set'}) {
		$vfs = $vfa->fetch_all_by_Slice_VariationSet($slice, $set);
	}
    elsif ($self->{'include_consequences'}) {
		$vfs = $vfa->fetch_all_by_Slice_constraint_with_TranscriptVariations($slice);
    }
    else {
    	$vfs = $vfa->fetch_all_by_Slice($slice);
    }
    return $self->_generic_Iterator($vfs);	
}

sub svfs_on_slice {
    my $self  = shift;
    my $slice = shift;
    my $svfs = $svfa->fetch_all_by_Slice($slice, 1);
    return $self->_generic_Iterator($svfs);
}

sub _generic_Iterator {
    my ($self, $vfs) = @_;
    my $iterator =  Bio::EnsEMBL::Utils::Iterator->new($vfs);   
    return $iterator;
}

sub print_gvf_header {
    my $self = shift;
    my $header = shift;
    print $fh $working_slice->gvf_header(%$header);

    # print a sequence-region line in the GVF file for each slice @risk
    print $fh '##sequence-region ', $working_slice->seq_region_name, ' ',$working_slice->start, ' ', $working_slice->end, "\n"; 
    for my $slice (@$slices) {
        print $fh '##sequence-region ', $slice->seq_region_name, ' ',$slice->start, ' ', $slice->end, "\n"; 
    }
}

sub print_variation {
    my $self  = shift;
    my $vf    = shift;
    my $attrs = shift;
    
    if ($self->{'include_failed'} || $self->{'just_failed'}) {
        my $desc = $self->failure_reason($vf->{_variation_id});
        if ($desc) {
            $attrs->{ensembl_failure_reason} = $desc;
        } else {
            return if $self->{'just_failed'};
        }
    }

    # if we get here, we want this vf included in the dump, so convert
    # it to GVF including any extra attributes defined above
    $attrs->{ID} = ++$id_count;
    my $gvf_line = $vf->to_gvf(
                        extra_attrs               => $attrs, 
                        include_consequences      => $self->{'include_consequences'},
                        include_coding_details    => $self->{'include_coding_details'},
                        include_global_maf        => $self->{'include_global_maf'},
                        include_validation_states => $self->{'include_validation_states'},
                        include_clinical_significance => $self->{'include_clinical_significance'},
                    );
    print $fh $gvf_line,  "\n" if $gvf_line;
}

my $prev_svs;

sub print_structural_variation {
    my $self  = shift;
    my $svf   = shift;
    my $attrs = shift; 
    return if $svf->var_class eq 'CNV_PROBE';

   if ($self->{'include_failed'} || $self->{'just_failed'}) {
        my $desc = $self->failure_reason($svf->{structural_variation_id});
        if ($desc) {
            $attrs->{ensembl_failure_reason} = $desc;
        } else {
            return if $self->{'just_failed'};
        }
    }
    my $coords = join '-', $svf->seq_region_name, $svf->seq_region_start, $svf->seq_region_end; 
                
    if (my $prev_coords = $prev_svs->{$svf->variation_name}) {
        warn "repeated SV: ".$svf->variation_name." coords 1: $prev_coords, coords 2: $coords\n" if $prev_coords ne $coords;
        return;
    }
    
    # we can now have SVs that map to multiple locations so we can't use the
    # feature's own identifier and we have to use a file-wide count as for
    # normal variations

    $attrs->{ID} = ++$id_count;
    my $gvf_line = $svf->to_gvf(extra_attrs => $attrs);     
    print $fh "$gvf_line\n" if $gvf_line;
    $prev_svs->{$svf->variation_name} = $coords if $gvf_line;
}

sub failure_reason {
    my $self         = shift;
    my $variation_id = shift;

    my $desc_id      = $failed_ids->{$variation_id};
    unless ($desc_id) {
        $desc_id = $structural_failed_ids->{$variation_id};
    }
    if ($desc_id) {
        return $failed_descs->{$desc_id};
    } else {
        return 0;
    }
}

1;
