use strict;
use warnings;
#object that contains the specific methods to dump data when there are chromosome coordinates from dbSNP (not contigs, as usual). 
#So far, this is the case for rat and chicken
package dbSNP::GenericChromosome;

use dbSNP::GenericContig;
use vars qw(@ISA);
use ImportUtils qw(debug load dumpSQL create_and_load);

@ISA = ('dbSNP::GenericContig');

sub variation_feature{
    my $self = shift;

    ### TBD not sure if variations with map_weight > 1 or 2 should be
    ### imported.
    
    debug("Dumping seq_region data");
    
    #only take chromosome coordinates
    dumpSQL($self->{'dbCore'}, qq{SELECT sr.seq_region_id, sr.name
				      FROM   seq_region sr
				      WHERE coord_system_id = 1});

    debug("Loading seq_region data");
    create_and_load($self->{'dbVariation'}, "tmp_seq_region", "seq_region_id", "name *");
    
    debug("Dumping SNPLoc data");
    
    my $tablename = $self->{'species_prefix'} . 'SNPContigLoc';
    dumpSQL($self->{'dbSNP'}, qq{SELECT snp_id, contig_chr, phys_pos_from, phys_pos,
				 IF(orientation, -1, 1)
				 FROM $tablename
				 $self->{'limit'}});
    
    
    debug("Loading SNPLoc data");
    
    create_and_load($self->{'dbVariation'}, "tmp_contig_loc", "snp_id i*", "chr *", "start i", 
		    "end i", "strand i");
    
    debug("Creating variation_feature data");
    
    $self->{'dbVariation'}->do(qq{INSERT INTO variation_feature 
				      (variation_id, seq_region_id,
				       seq_region_start, seq_region_end, seq_region_strand,
				       variation_name)
				      SELECT v.variation_id, ts.seq_region_id, tcl.start, tcl.end,
				      tcl.strand, v.name
				      FROM   variation v, tmp_contig_loc tcl, tmp_seq_region ts
				      WHERE  v.snp_id = tcl.snp_id
				      AND    tcl.chr = ts.name});
    
    $self->{'dbVariation'}->do("DROP TABLE tmp_contig_loc");
    $self->{'dbVariation'}->do("DROP TABLE tmp_seq_region");
}

1;
