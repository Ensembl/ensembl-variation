use strict;
use warnings;
#object that contains the specific methods to dump data when there specie is a mosquito (not contigs, as usual). 
package dbSNP::Mosquito;

use dbSNP::GenericContig;
use vars qw(@ISA);
use ImportUtils qw(debug load dumpSQL create_and_load);

@ISA = ('dbSNP::GenericContig');


sub variation_feature{
    my $self = shift;

    my %scaff; #contains relation AAAB -> NW
    open(IN,"gzip -dc dbSNP/scaff_NW_chr_CRA_AAAB.gz |") || die "Could not get file with scaff mapping data: $!\n";
    my @line;
    while (<IN>){
	chomp;
	@line = split /\t/;
	$scaff{$line[3]} = $line[0];
    }
    ### TBD not sure if variations with map_weight > 1 or 2 should be
    ### imported.
    debug("Dumping seq_region data");
    my $sth = $self->{'dbCore'}->dbc()->prepare(qq{SELECT sr.seq_region_id, sr.name
						 FROM   seq_region sr
						 WHERE coord_system_id = 3});
    $sth->execute();
    open(FH,">" . $self->{'tmpdir'} . '/' . $self->{'tmpfile'}); #open the file with the data dump from the core database
    my $row;
    while($row = $sth->fetchrow_arrayref()) {
	my @row = @$row;
	$row[1] = $scaff{$row[1]} if (defined $scaff{$row[1]});
	@row = map {(defined($_)) ? $_ : '\N'} @row;  # convert undefined to NULL;
	print FH join("\t", @row), "\n";
    }

    $sth->finish();
    close(FH);

    debug("Loading seq_region data");
    create_and_load($self->{'dbVariation'}, "tmp_seq_region", "seq_region_id", "name *");
    
    debug("Dumping SNPLoc data");
    
    my $tablename = $self->{'species_prefix'} . 'SNPContigLoc';
    dumpSQL($self->{'dbSNP'}, qq{SELECT snp_id, contig_acc, asn_from,
				 IF(loc_type = 3,  asn_from - 1, asn_to), # 3 = between
				 IF(orientation, -1, 1)
				 FROM $tablename
				 $self->{'limit'}});
    
    
    debug("Loading SNPLoc data");
    
    create_and_load($self->{'dbVariation'}, "tmp_contig_loc", "snp_id i*", "chr *", "start i", 
		    "end i", "strand i");
    
    #creating the temporary table with the genotyped variations
    dumpSQL($self->{'dbVariation'},qq{SELECT DISTINCT variation_id
					  FROM individual_genotype
				      });
    create_and_load($self->{'dbVariation'},'tmp_genotyped_var',"variation_id *");

    debug("Creating variation_feature data");
    
    $self->{'dbVariation'}->do(qq{INSERT INTO variation_feature 
				      (variation_id, seq_region_id,
				       seq_region_start, seq_region_end, seq_region_strand,
				       variation_name, flags, source_id, validation_status)
				      SELECT v.variation_id, ts.seq_region_id, tcl.start, tcl.end,
				      tcl.strand, v.name, IF(tgv.variation_id,'genotyped',NULL), v.source_id, v.validation_status
				      FROM   variation v LEFT JOIN tmp_genotyped_var tgv ON v.variation_id = tgv.variation_id, tmp_contig_loc tcl, tmp_seq_region ts
				      WHERE  v.snp_id = tcl.snp_id
				      AND    tcl.chr = ts.name});
    
    $self->{'dbVariation'}->do("DROP TABLE tmp_contig_loc");
    $self->{'dbVariation'}->do("DROP TABLE tmp_seq_region");
    $self->{'dbVariation'}->do("DROP TABLE tmp_genotyped_var");
}

1;
