use strict;
use warnings;
#object that loading mapping data into the variation_feature table

package dbSNP::MappingChromosome;

use dbSNP::GenericContig;
use vars qw(@ISA);
use ImportUtils qw(debug load create_and_load dumpSQL);

@ISA = ('dbSNP::GenericContig');


sub variation_feature{
  my $self = shift;
    
  my (%rec,%source,%status,%rec_pos,%rec_seq_region);

  debug("Dumping Variation");

  # Create hashes of variation_id's, source_id's etc keyed by variation name
  my $sth = $self->{'dbVariation'}->prepare (qq{
SELECT variation_id, name, source_id, validation_status
FROM variation});
  $sth->execute();
  while( my ($variation_id, $name, $source_id, $validation_status) 
         = $sth->fetchrow_array()) {
    $rec{$name} = $variation_id;
    $source{$name} = $source_id;
    $status{$name} = $validation_status;
  }
  $sth->finish();

  # Create hash of top-level seq_region_id keyed by seq_region_name
  my ($seq_region_id,$seq_region_name);
  my $sth1 = $self->{'dbCore'}->dbc->prepare (qq{
select sr.seq_region_id, sr.name 
from   seq_region sr, seq_region_attrib sra, attrib_type at 
where  sr.seq_region_id=sra.seq_region_id 
and    sra.attrib_type_id = at.attrib_type_id 
and    (at.code="toplevel" or at.code="non_ref") } );
  $sth1->execute();
  $sth1->bind_columns(\$seq_region_id,\$seq_region_name);
  while ($sth1->fetch) {
    $rec_seq_region{$seq_region_name} = $seq_region_id;
  }
  $sth1->finish();


  open (FH, ">" . $self->{'tmpdir'} . "/" . $self->{'tmpfile'} );

  # Concatenate all mapping_file_N files into a single result file
  my $mapping_dir = $self->{'mapping_file'};
  system("cat $mapping_dir/mapping_file* > $mapping_dir/all_mapping_file");
  open (IN, "$mapping_dir/all_mapping_file") 
      or die "can't open mapping_file:$!";
  
  # Process results file
  while (<IN>) {
    chomp;
    next if /^more |^PARSING/;
    s/^MORE_HITS//;
    my ($ref_id, $slice_name, $start, $end, $strand, $ratio) = split;
    
    # Skip mappings where %ID < 50%
    next if $ratio <0.5;
      
    my ($coord_sys,$assembly,$seq_region_name,$seq_region_start,$seq_region_end,$version);

    # Parse the slice name 
    if( $slice_name =~ /-/ ){
      # Old style; <SRname>-<SRstart>-<SRend>
      ($seq_region_name,
       $seq_region_start,
       $seq_region_end) 
          = split /\-/, $slice_name;
    }
    elsif( $slice_name =~ /\:/ ){
      # New style: uses Bio::EnsEMBL::Slice->name
      ($coord_sys,
       $assembly,
       $seq_region_name,
       $seq_region_start,
       $seq_region_end,
       $version) 
          = split /\:/, $slice_name
    }
    else{
      # Asume that we have a seq_region_name only
      $seq_region_name  = $slice_name;
      $seq_region_start = 1;
      $seq_region_end   = '';
    }

    # Get the coreDB seq_region_id from its name
    my $seq_region_id = $rec_seq_region{$seq_region_name}; # Find in cache
    if (! $seq_region_id) {
      # Not in cache - query database 
      my $sth = $self->{'dbCore'}->dbc->prepare (qq{SELECT seq_region_id from seq_region where name = ?});
      $sth->execute("$seq_region_name");
      $seq_region_id = $sth->fetchrow_array();
    }
    unless( $seq_region_id ){
      debug( "Seq region $seq_region_name not found in DB" );
      next;
    }

    # Calculate the seq_region_start for the variation_feature
    my $new_seq_region_start = $start;
    my $new_seq_region_end   = $end;
    if( $seq_region_start ){ # Apply offset if available
      $new_seq_region_start = $seq_region_start + $start -1;
      $new_seq_region_end   = $seq_region_start + $end   -1;
    }
    
    next if ( $rec_pos{$ref_id}{$seq_region_id}{$new_seq_region_start}{$new_seq_region_end}); # Skip if variation_feature already processed

    # Print to output file
    print FH join
        ( "\t",
          $seq_region_id,           # variation_feature.seq_region_id
          $new_seq_region_start,    # variation_feature.seq_region_start
          $new_seq_region_end,      # variation_feature.seq_region_end
          $strand,                  # variation_feature.seq_region_strand
          $rec{$ref_id},            # variation_feature.variation_id
          $ref_id,                  # variation_feature.variation_name
          $source{$ref_id},         # variation_feature.source_id
          $status{$ref_id}||'NULL', # variation_feature.validation_status
          ) . "\n";

    # Flag as processed
    $rec_pos{$ref_id}{$seq_region_id}{$new_seq_region_start}{$new_seq_region_end}=1;
  }
    
  close IN;
  close FH;
  
    
  debug("Creating genotyped variations");

  # Creating temporary variation_feature table
  create_and_load( $self->{'dbVariation'}, 
                   "tmp_variation_feature",
                   "seq_region_id",
                   "seq_region_start",
                   "seq_region_end",
                   "seq_region_strand",
                   "variation_id *",
                   "variation_name", 
                   "source_id", 
                   "validation_status" );

  #creating the temporary table with the genotyped variations
   $self->{'dbVariation'}->do(qq{
CREATE TABLE tmp_genotyped_var 
SELECT DISTINCT variation_id FROM tmp_individual_genotype_single_bp});
   $self->{'dbVariation'}->do(qq{
CREATE UNIQUE INDEX variation_idx ON tmp_genotyped_var (variation_id)});
   $self->{'dbVariation'}->do(qq{
INSERT IGNORE INTO  tmp_genotyped_var 
SELECT DISTINCT variation_id FROM individual_genotype_multiple_bp});
  
  # Copy from tmp variation_feature to final variation_feature
  $self->{'dbVariation'}->do(qq{
INSERT INTO variation_feature
      (variation_id,
       seq_region_id, 
       seq_region_start, 
       seq_region_end, 
       seq_region_strand,
       variation_name, 
       flags, 
       source_id, 
       validation_status)
SELECT tv.variation_id, 
       tv.seq_region_id, 
       tv.seq_region_start, 
       tv.seq_region_end,
       tv.seq_region_strand, 
       tv.variation_name, 
       IF(tgv.variation_id,'genotyped',NULL), 
       tv.source_id,
       tv.validation_status
FROM   tmp_variation_feature tv LEFT JOIN 
       tmp_genotyped_var tgv ON tv.variation_id = tgv.variation_id });

  $self->{'dbVariation'}->do(qq{DROP TABLE tmp_variation_feature});
  $self->{'dbVariation'}->do(qq{DROP TABLE tmp_genotyped_var});
}

1;
