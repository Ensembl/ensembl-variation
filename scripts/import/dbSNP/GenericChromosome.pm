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
    dumpSQL($self->{'dbSNP'}, qq{SELECT snp_id, contig_chr, phys_pos_from, 
				 IF(loc_type = 3,  phys_pos_from-1,
				    IF (loc_type = 1, substring_index(phys_pos, '..', -1),phys_pos)),
				 IF(orientation, -1, 1)
				 FROM $tablename
				 $self->{'limit'}});
    
    
    debug("Loading SNPLoc data");
    
    create_and_load($self->{'dbVariation'}, "tmp_contig_loc", "snp_id i*", "chr *", "start i", 
		    "end i", "strand i");

    debug("Creating genotyped variations");
    #creating the temporary table with the genotyped variations
    $self->{'dbVariation'}->do(qq{CREATE TABLE tmp_genotyped_var SELECT DISTINCT variation_id FROM individual_genotype_single_bp});
    $self->{'dbVariation'}->do(qq{INSERT IGNORE INTO tmp_genotyped_var SELECT DISTINCT variation_id FROM individual_genotype_multiple_bp});
    $self->{'dbVariation'}->do(qq{CREATE INDEX variation_idx ON tmp_genotyped_var (variation_id)});
 
    debug("Creating tmp_variation_feature data");
    
    dumpSQL($self->{'dbVariation'},qq{SELECT v.variation_id, ts.seq_region_id, tcl.start, tcl.end, tcl.strand, v.name, v.source_id, v.validation_status
					  FROM variation v, tmp_contig_loc tcl, tmp_seq_region ts
					  WHERE v.snp_id = tcl.snp_id
					  AND tcl.chr = ts.name });

    create_and_load($self->{'dbVariation'},'tmp_variation_feature',"variation_id *","seq_region_id", "seq_region_start", "seq_region_end", "seq_region_strand", "variation_name", "source_id", "validation_status");
    
    debug("Dumping data into variation_feature table");
    $self->{'dbVariation'}->do(qq{INSERT INTO variation_feature (variation_id, seq_region_id,seq_region_start, seq_region_end, seq_region_strand,
								 variation_name, flags, source_id, validation_status)
				      SELECT tvf.variation_id, tvf.seq_region_id, tvf.seq_region_start, tvf.seq_region_end, tvf.seq_region_strand,
				      tvf.variation_name,IF(tgv.variation_id,'genotyped',NULL), tvf.source_id, tvf.validation_status
				      FROM tmp_variation_feature tvf LEFT JOIN tmp_genotyped_var tgv ON tvf.variation_id = tgv.variation_id
				  });
    $self->{'dbVariation'}->do("DROP TABLE tmp_contig_loc");
    $self->{'dbVariation'}->do("DROP TABLE tmp_seq_region");
    $self->{'dbVariation'}->do("DROP TABLE tmp_genotyped_var");
    $self->{'dbVariation'}->do("DROP TABLE tmp_variation_feature");
    #for the chicken, delete 13,000 SNPs that cannot be mapped to EnsEMBL coordinate
    if ($self->{'species_prefix'} eq 'gga'){
	$self->{'dbVariation'}->do("DELETE FROM variation_feature WHERE seq_region_end = -1");
    }
}

1;
