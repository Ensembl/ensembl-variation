use strict;
use warnings;
#object that contains the specific methods to dump data when the specie is a HUMAN (adds HGVbase and TSC information). 
package dbSNP::Human;

use dbSNP::GenericContig;
use vars qw(@ISA);
use ImportUtils qw(debug load dumpSQL create_and_load);

@ISA = ('dbSNP::GenericContig');

sub dump_dbSNP{
    my $self = shift;
    #first, dump all dbSNP data as usual
    $self->SUPER::dump_dbSNP();
    #then, get HGVbase IDs from Yuans file
    $self->dump_HGVbaseIDs();
    #and finally, get TSC data from dbSNP
    $self->dump_TSCIDs();    
}

#specific function to get the HGVbase IDs from a file provided by Yuan and add them to the variation_synonym table
sub dump_HGVbaseIDs{
    my $self = shift;
    
    #copy the file with the rs-> HGVbaseID information to the temp folder
    system('gunzip -c','dbSNP/rs_hgvbase.txt.gz',$self->{'tmpdir'} . '/' . $self->{'tmpfile'});
    debug("Loading HGVbase data");
    create_and_load($self->{'dbVariation'},"tmp_rs_hgvbase","rsID","HGVbaseID");
    #add a new source to the Source table
    $self->{'dbVariation'}->do(qq{INSERT INTO source (name,version) values ('HGVbase',15)
				  });
    debug("Adding HGVbaseIDs to synonym table");
    my $source_id = $self->{'dbVariation'}->dbh()->{'mysql_insertid'}; #get the last autoinc id from the database (the one from the HGVbase source)
    #add the HGVbaseIDs to the Variation_Synonym table
    $self->{'dbVariation'}->do(qq{INSERT INTO variation_synonym (variation_id,source_id,name)
				      SELECT v.variation_id, $source_id, trh.HGVbaseID
				      FROM variation v, tmp_rs_hgvbase trh
				      WHERE v.name = trh.rsID
				  });
    #and finally, remove the temporary table
    $self->{'dbVariation'}->do(qq{DROP TABLE tmp_rs_hgvbase
				  });
}

#specific function to get the TSC IDs from dbSNP and add them to the variation_synonym table
sub dump_TSCIDs{
    my $self = shift;
    
    #add the TSC source to the table
    $self->{'Variation'}->do(qq{INSERT INTO source (name,version) values ('TSC',1)
    });
    my $source_id = $self->{'dbVariation'}->dbh()->{'mysql_insertid'}; #get the last autoinc id in the database (the one from the TSC source)
    #and finally add the TSC ids to the synonyms table
    debug("Dumping TSC information from dbSNP");
    dumpSQL($self->{'dbSNP'}, qq{SELECT concat('rs',ss.snp_id), $source_id, s.loc_snp_id 
				     FROM SubSNP s, SNPSubSNPLink ss 
				     WHERE ss.subsnp_id = s.subsnp_id 
				     AND s.loc_snp_id like 'TSC%'
				 }
	    );
    debug("Loading TSC ids into temporary table");
    create_and_load($self->{'dbVariation'},"tmp_rs_TSC","rsID","sourceID","TSCid");
    $self->{'Variation'}->do(qq{ INSERT INTO variation_synonym (variation_id, source_id, name)
				     SELECT v.variation_id, trt.source_id, trt.TSCid 
				     FROM Variation v, tmp_rs_TSC trt
				     WHERE v.name = trt.rsID 
				}
			     );
}
