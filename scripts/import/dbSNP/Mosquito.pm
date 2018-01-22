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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

use warnings;
#object that contains the specific methods to dump data when there specie is a mosquito (not contigs, as usual). 
package dbSNP::Mosquito;

use dbSNP::GenericContig;
use vars qw(@ISA);
use ImportUtils qw(debug load dumpSQL create_and_load);

@ISA = ('dbSNP::GenericContig');

# will be used to store the 2 different strains and the description
our %STRAINS = ( 'PEST' => 'SNPs detected by comparing mosquito pest and mopti reads against Celera Anopheles Gambiae build 2 Genomic Reference Sequence using ssahaSNP',
		 'MOPTI' => 'SNPs detected by comparing mosquito pest and mopti reads against Celera Anopheles Gambiae build 2 Genomic Reference Sequence using ssahaSNP');

sub dump_dbSNP{
    my $self = shift;
    #first, dump all dbSNP data as usual
    $self->SUPER::dump_dbSNP();
    #then, get strain information
    $self->add_strains();    
}



sub add_strains{
    my $self = shift;

    my $pop_MOPTI; #population_id for the mopti strain
    my $pop_PEST;  #population_id for the pest strain
    #first of all, add to the population table the 2 strains: PEST and MOPTI
    foreach my $population_name (keys %STRAINS){
      $self->{'dbVariation'}->do(qq{INSERT INTO sample (name,description)
				    VALUES ('$population_name', '$STRAINS{$population_name}')});
      #store in the same hash the id for the population
      $STRAINS{$population_name} = $self->{'dbVariation'}->{'mysql_insertid'};
    }
    #and copy the data from the sample to the Population table
    debug("Loading population table with data from sample");
    $self->{'dbVariation'}->do(qq{INSERT INTO population (sample_id) 
				  SELECT sample_id
				  FROM sample
				 });
    

    #then, get from dbSNP the relation between ssId and mopti/pest strain
    debug("Dumping strain information from dbSNP");

    my $stmt = qq{
                  SELECT 
                    subsnp_id, 
                    loc_snp_id
	          FROM 
	            SubSNP
                 };
    dumpSQL($self->{'dbSNP'},$stmt);
    
    create_and_load($self->{'dbVariation'},'tmp_strain', "subsnp_id i*", "strain");

    debug("Updating allele table with strain alleles");
    #first, update the MOPTI alleles
    $self->{'dbVariation'}->do(qq{UPDATE allele , variation_synonym, tmp_strain
				      SET allele.sample_id = $STRAINS{'MOPTI'}, allele.frequency = 1
				      WHERE allele.variation_id = variation_synonym.variation_id
				      AND SUBSTRING(variation_synonym.name,3) = tmp_strain.subsnp_id
				      AND tmp_strain.strain like '%mopti%'
				  });

    #and do the same for the PEST alleles
    $self->{'dbVariation'}->do(qq{UPDATE allele , variation_synonym, tmp_strain
				      SET allele.sample_id = $STRAINS{'PEST'}, allele.frequency = 1
				      WHERE allele.variation_id = variation_synonym.variation_id
				      AND SUBSTRING(variation_synonym.name,3) = tmp_strain.subsnp_id
				      AND tmp_strain.strain like '%pest%'
				  });
    
    #and finally, drop the temporary table
    $self->{'dbVariation'}->do(qq{DROP TABLE tmp_strain});
}


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
						 FROM   seq_region sr, coord_system cs
						 WHERE sr.coord_system_id = cs.coord_system_id
						 AND cs.name = "scaffold"});
    $sth->execute();
    open(FH,">" . $self->{'tmpdir'} . '/' . $self->{'tmpfile'}); #open the file with the data dump from the core database

    my ($tablename1,$tablename2,$row);

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

    
    my ($assembly_version) =  $self->{'assembly_version'} =~ /^[a-zA-Z]+(\d+)\.*.*$/; $assembly_version=2;
    print "assembly_version again is $assembly_version\n";


    my $stmt = qq{
                  SELECT 
                    name 
                  FROM 
                    $self->{'snp_dbname'}..sysobjects 
                  WHERE 
                    name LIKE '$self->{'dbSNP_version'}\_SNPContigLoc\_$assembly_version\__'
               };                 
    my $sth1 = $self->{'dbSNP'}->prepare($stmt);
    $sth1->execute();

    while($row = $sth1->fetchrow_arrayref()) {
      $tablename1 = $row->[0];
    }
    
    $stmt = qq{
               SELECT 
                 name 
               FROM 
                 $self->{'snp_dbname'}..sysobjects 
               WHERE 
                 name LIKE '$self->{'dbSNP_version'}\_ContigInfo\_$assembly_version\__'
            };
    my $sth2 = $self->{'dbSNP'}->prepare($stmt);
    $sth2->execute();

    while($row = $sth2->fetchrow_arrayref()) {
      $tablename2 = $row->[0];
    }
    print "table_name1 is $tablename1 table_name2 is $tablename2\n";

    $stmt = "SELECT ";
    if ($self->{'limit'}) {
      $stmt .= "TOP $self->{'limit'} ";
    }
    $stmt .= qq{
                  t1.snp_id AS sorting_id, 
                  t2.contig_acc,
		  t1.lc_ngbr+2,t1.rc_ngbr,
		  CASE WHEN
		    orientation
		  THEN
		    -1
		  ELSE
		    1
		  END
		FROM   
		  $tablename1 t1, 
		  $tablename2 t2
		WHERE 
		  t1.ctg_id=t2.ctg_id
	        ORDER BY
		  sorting_id ASC  
               };
    dumpSQL($self->{'dbSNP'},$stmt);
    
    debug("Loading SNPLoc data");
    
    create_and_load($self->{'dbVariation'}, "tmp_contig_loc", "snp_id i*", "contig *", "start i", 
		    "end i", "strand i");
    
    #creating the temporary table with the genotyped variations

     $self->{'dbVariation'}->do(qq{CREATE TABLE tmp_genotyped_var SELECT DISTINCT variation_id FROM tmp_individual_genotype_single_bp});
     $self->{'dbVariation'}->do(qq{CREATE UNIQUE INDEX variation_idx ON tmp_genotyped_var (variation_id)});
     $self->{'dbVariation'}->do(qq{INSERT IGNORE INTO  tmp_genotyped_var SELECT DISTINCT variation_id FROM individual_genotype_multiple_bp});



    debug("Creating variation_feature data");
    
    $self->{'dbVariation'}->do(qq{INSERT INTO variation_feature 
				      (variation_id, seq_region_id,
				       seq_region_start, seq_region_end, seq_region_strand,
				       variation_name, flags, source_id, validation_status)
				      SELECT v.variation_id, ts.seq_region_id, tcl.start, tcl.end,
				      tcl.strand, v.name, IF(tgv.variation_id,'genotyped',NULL), v.source_id, v.validation_status
				      FROM   variation v LEFT JOIN tmp_genotyped_var tgv ON v.variation_id = tgv.variation_id, tmp_contig_loc tcl, tmp_seq_region ts
				      WHERE  v.snp_id = tcl.snp_id
				      AND    tcl.contig = ts.name});
    
    $self->{'dbVariation'}->do("DROP TABLE tmp_contig_loc");
    $self->{'dbVariation'}->do("DROP TABLE tmp_seq_region");
    $self->{'dbVariation'}->do("DROP TABLE tmp_genotyped_var");
}

1;
