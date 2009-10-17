use strict;
use warnings;
#object that contains the specific methods to dump data when the specie is a HUMAN (adds HGVbase and TSC information). 
package dbSNP::EnsemblIds;

#use dbSNP::GenericContig;
use dbSNP::GenericChromosome;

use vars qw(@ISA);
use ImportUtils qw(debug load dumpSQL create_and_load);

#@ISA = ('dbSNP::GenericContig');
@ISA = ('dbSNP::GenericChromosome');

sub dump_dbSNP{
    my $self = shift;
    #first, dump all dbSNP data as usual
    #$self->SUPER::dump_dbSNP();
    #then, get ENS IDs from dbSNP
    $self->dump_ENSIDs() if $self->{'dbCore'}->species =~ /hum|rat|mouse|platypus|tetraodon/i;
    $self->dump_AFFYIDs() if $self->{'dbCore'}->species =~ /hum/i;
}

#get ENSIDs from dbSNP
sub dump_ENSIDs{
  my $self = shift;

  my %ens_batch_name = ('rat' => ['RAT_COMPUTATIONAL_CELERA','RAT_COMPUTATIONAL_STAR'],
			'platypus' => ['2008FEB_platypus_assembly','2008FEB_platypus_reads'],
		       );

  #add the ENS source to the table
  foreach my $species (keys %ens_batch_name) {
    if ($self->{'dbCore'}->species =~ /$species/i) {
      foreach my $batch_name (@{$ens_batch_name{$species}}) {
	debug("Processing $batch_name...");
	$self->{'dbVar'}->do(qq{INSERT INTO source (name,version,description) values ("ENSEMBL_$batch_name",1,"Variation called by ENSEMBL")});
	my $source_id = $self->{'dbVar'}->db_handle->{'mysql_insertid'}; #get the last autoinc id in the database (the one from the ENS source)
	#and finally add the ENS ids to the synonyms table
	debug("Dumping ENS information from dbSNP");
	dumpSQL($self->{'dbSNP'}, qq{SELECT concat('rs',s.snp_id), s.subsnp_id, $source_id, s.loc_snp_id
				     FROM SubSNP s, Batch b
				     WHERE s.batch_id = b.batch_id
                                     AND b.handle = 'ENSEMBL'
                                     AND b.loc_batch_id = "$batch_name"
				 }
	       );
	debug("Loading ENS ids into temporary table");
	create_and_load($self->{'dbVar'},"tmp_rs_ENS_$batch_name","rsID *","ss_id i","source_id i","ENSid");
	debug("Loading ENS ids into variation_synonym table");
	$self->{'dbVar'}->do(qq{ INSERT IGNORE INTO variation_synonym (variation_id, subsnp_id, source_id, name)
				     SELECT v.variation_id, trt.ss_id, trt.source_id, trt.ENSid 
				     FROM variation v, tmp_rs_ENS_$batch_name trt
				     WHERE v.name = trt.rsID 
				}
			    );
	#and finally, remove the temporary table
	#$self->{'dbVar'}->do(qq{DROP TABLE tmp_rs_ENS_$batch_name});
      }
    }
  }
}


sub dump_AFFYIDs{

  my $self = shift;
  my ($source_name,$set_name);

  debug("Dumping AFFY information from dbSNP");
=head
  dumpSQL($self->{'dbSNP'}, qq{SELECT concat('rs',s.snp_id), s.loc_snp_id, loc_batch_id
				     FROM SubSNP s, Batch b
                                     WHERE s.batch_id = b.batch_id
                                     AND b.handle = "AFFY"
				 }
	 );
  debug("Loading  ids into temporary table");
  create_and_load($self->{'dbVar'},"tmp_rs_AFFY","rsID *","AFFYid", "affy_name");
=cut
  foreach my $table ("yuan_aff_100k_var_46","yuan_aff_500k_var_46","yuan_aff_genome6_var_47") {
    if ($table =~ /100k/i) {
      $source_name = "Affy GeneChip 100K Array";
      $set_name = "Mapping50K";
    }
    elsif ($table =~ /500k/i) {
      $source_name = "Affy GeneChip 500K Array";
      $set_name = "Mapping250K";
    }
    elsif ($table =~ /genome6/i) {
      $source_name = "Affy GenomeWideSNP_6.0";
      $set_name = "6.0";
    }

    debug("Creating name_pair table with source_name $source_name...");

    $self->{'dbVar'}->do(qq{CREATE TABLE $table\_name_pair like $table.name_pair});
    $self->{'dbVar'}->do(qq{insert into $table\_name_pair select * from $table.name_pair});
    $self->{'dbVar'}->do(qq{insert ignore into $table\_name_pair
                            select a.AFFYid as affy_name,a.rsID as rs_name 
                            from tmp_rs_AFFY a, $table.name_pair c 
                            where a.AFFYid=c.affy_name});

    $self->{'dbVar'}->do(qq{UPDATE yuan_aff_genome6_var_47_name_pair
                            SET affy_name = REPLACE(affy_name,'AFFY_6_1M_','')});

    my $source_id_ref = $self->{'dbVar'}->db_handle->selectall_arrayref(qq{
                     SELECT source_id from source where name = "$source_name"});
    my $source_id = $source_id_ref->[0][0];

    if (!$source_id) {
      $self->{'dbVar'}->do(qq{insert into source (name) values("$source_name")});
      $source_id = $self->{'dbVar'}->db_handle->{'mysql_insertid'};
    }

    debug("Inserting in variation_synonym table from $table\_name_pair...");
    $self->{'dbVar'}->do(qq{ INSERT IGNORE INTO variation_synonym_test (variation_id, source_id, name)
                                     SELECT v.variation_id, $source_id as source_id, t.affy_name as name 
                                     FROM variation v, $table\_name_pair t
                                     WHERE v.name = t.rs_name
                                }
			);

    #update rs_ID to rsCurrent from rsHigh
    #$self->{'dbVar'}->do(qq{update tmp_rs_AFFY_test t, rsHist h set t.rsID=h.rsCurrent where t.rsID=h.rsHigh});
    debug("Inserting in variation_synonym table from table  tmp_rs_AFFY with $source_name...");
    $self->{'dbVar'}->do(qq{ INSERT IGNORE INTO variation_synonym_test (variation_id, source_id, name)
				     SELECT v.variation_id, $source_id as source_id, t.AFFYid as name 
				     FROM variation v, tmp_rs_AFFY_test t
				     WHERE v.name = t.rsID
                                     AND t.affy_name like "%$set_name%"
				}
			);

   #and finally, remove the temporary table
   #$self->{'dbVar'}->do(qq{DROP TABLE $table\_name_pair});

  }
}

1;
