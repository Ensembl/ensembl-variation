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
#object that loading mapping data into the variation_feature table

package dbSNP::MappingChromosome;

use dbSNP::GenericContig;
use vars qw(@ISA);
use ImportUtils qw(debug load create_and_load dumpSQL);

@ISA = ('dbSNP::GenericContig');


sub variation_feature{
  my $self = shift;

  my (%rec,%source,%status,%rec_pos,%rec_seq_region,$three_mappings);


  debug("Dumping Variation");

  my ($variation_id,$name,$source_id,$validation_status,$variation_feature_id,$allele_string);

  debug("Move current variation_feature table to variation_feature_before_mapping");
  $self->{'dbVar'}->do(qq{RENAME TABLE variation_feature TO variation_feature_before_mapping});
  $self->{'dbVar'}->do(qq{CREATE TABLE variation_feature LIKE variation_feature_before_mapping});

  #we need to create a seq_region table
  $self->{'dbVar'}->do(qq{CREATE TABLE IF NOT EXISTS seq_region(
                               seq_region_id int(10),
                               name varchar(40),
                               primary key (seq_region_id),
                               unique key (name))
 });

  # Create hash of top-level seq_region_id keyed by seq_region_name
  my ($seq_region_id,$seq_region_name);
  my $sth1 = $self->{'dbCore'}->dbc->prepare (qq{
select sr.seq_region_id, sr.name 
from   seq_region sr, seq_region_attrib sra, attrib_type at 
where  sr.seq_region_id=sra.seq_region_id 
and    sra.attrib_type_id = at.attrib_type_id 
and    at.code="toplevel" } );
  $sth1->execute();
  $sth1->bind_columns(\$seq_region_id,\$seq_region_name);
  while ($sth1->fetch) {
    $rec_seq_region{$seq_region_name} = $seq_region_id;
    $self->{'dbVar'}->do(qq{INSERT IGNORE INTO seq_region(seq_region_id,name)values($seq_region_id,"$seq_region_name")});
  }
  $sth1->finish();


  debug("Reading Mapping file");
#  open (FH, ">" . $self->{'tmpdir'} . "/" . $self->{'tmpfile'} );

  my $mapping_file_dir = $self->{'mapping_file_dir'};

  debug("Mapping file is $mapping_file_dir");
  # If mapping_file_dir is a file, read from it
  if (-f "$mapping_file_dir") {
    open (IN, "$mapping_file_dir") or die "can't open mapping_file:$!";
    my @ar = split /\//, $mapping_file_dir;
    $ar[-1] = "THREE_MAPPINGS";
    $three_mappings = join "/", @ar;
    print "three_mappings is $three_mappings\n";
  }
  # Concatenate all mapping_file_N files into a single result file
  elsif (-d "$mapping_file_dir") {
    system("cat $mapping_file_dir/mapping_file* > $mapping_file_dir/all_mapping_file");
    system("cat $mapping_file_dir/THREE_MAPPINGS > $mapping_file_dir/THREE_MAPPINGS");
    $three_mappings = "$mapping_file_dir/THREE_MAPPINGS";
    open (IN, "$mapping_file_dir/all_mapping_file") or die "can't open mapping_file:$!";
  }


  #create table three_mappings and failed_variation
  $self->{'dbVar'}->do(qq{CREATE TABLE IF NOT EXISTS THREE_MAPPINGS (name varchar(15), unique key name(name))});
  system("cp  $three_mappings " . $self->{'tmpdir'} . "/" . $self->{'tmpfile'});
  load($self->{'dbVar'}, "THREE_MAPPINGS", "name");
  
  $self->{'dbVar'}->do(qq{insert ignore into failed_variation (variation_id,failed_description_id) select v.variation_id,1 from variation v, THREE_MAPPINGS t where v.name=t.name});
  open (FH, ">" . $self->{'tmpdir'} . "/" . $self->{'tmpfile'} );

  # Process results file
  while (<IN>) {
    chomp;
    next if /^more |^PARSING/;
    s/^MORE_HITS//;
    my ($ref_id, $slice_name, $start, $end, $strand, $ratio) =split;

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
    $strand = ($strand eq "+") ? 1 : -1;
    #print "$ref_id\n";
    #next if ( $rec_pos{$ref_id}{$seq_region_id}{$new_seq_region_start}{$new_seq_region_end}); # Skip if variation_feature already processed


    # Print to output file
    print FH join
        ( "\t",
          $seq_region_id,           # variation_feature.seq_region_id
          $new_seq_region_start,    # variation_feature.seq_region_start
          $new_seq_region_end,      # variation_feature.seq_region_end
          $strand,                  # variation_feature.seq_region_strand
          $ref_id,                  # variation_feature.variation_name
          ) . "\n";

    # Flag as processed
    #$rec_pos{$ref_id}{$seq_region_id}{$new_seq_region_start}{$new_seq_region_end}=1; needs too much memory to hashing
  }

  close IN;
  close FH;

  debug("Creating tmp variation_features");

  # Creating temporary variation_feature table
  create_and_load( $self->{'dbVar'}, 
                   "tmp_variation_feature",
                   "seq_region_id i*",
                   "seq_region_start i*",
                   "seq_region_end",
                   "seq_region_strand",
                   "variation_name *", 
                   );

  # Copy from tmp variation_feature to final variation_feature, think about multi-mapping


  $self->{'dbVar'}->do(qq{
CREATE TABLE tmp_genotyped_var 
SELECT DISTINCT variation_id FROM tmp_individual_genotype_single_bp});
   $self->{'dbVar'}->do(qq{
CREATE UNIQUE INDEX variation_idx ON tmp_genotyped_var (variation_id)});
   $self->{'dbVar'}->do(qq{
INSERT IGNORE INTO  tmp_genotyped_var 
SELECT DISTINCT variation_id FROM individual_genotype_multiple_bp});



  # Copy from tmp variation_feature to final variation_feature

  $self->{'dbVar'}->do(qq{
INSERT INTO variation_feature
      (variation_id,
       seq_region_id, 
       seq_region_start, 
       seq_region_end, 
       seq_region_strand,
       variation_name,
       source_id, 
       validation_status)
SELECT v.variation_id, 
       tv.seq_region_id, 
       tv.seq_region_start, 
       tv.seq_region_end,
       tv.seq_region_strand, 
       tv.variation_name, 
       v.source_id,
       v.validation_status
FROM   variation v, tmp_variation_feature tv
WHERE v.name = tv.variation_name
AND   v.variation_id=vf.variation_id
});


  debug("copy allele_string,flags from table variation_feature_before_mapping if it's null");

  $self->{'dbVar'}->do(qq{UPDATE variation_feature vf, variation_feature_before_mapping vfm 
                          SET vf.allele_string = vfm.allele_string
                          WHERE vf.variation_id=vfm.variation_id});

  $self->{'dbVar'}->do(qq{UPDATE variation_feature vf, tmp_genotyped_var tgv
                          SET vf.flags = 'genotyped'
                          WHERE vf.variation_id=tgv.variation_id});

  #$self->{'dbVar'}->do(qq{DROP TABLE tmp_variation_feature});
  #$self->{'dbVar'}->do(qq{DROP TABLE tmp_genotyped_var});
}

1;
