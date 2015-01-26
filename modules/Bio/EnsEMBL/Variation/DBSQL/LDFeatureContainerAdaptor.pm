=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut

#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::LDFeatureContainerAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::LDFeatureContainerAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $sa = $reg->get_adaptor("human","core","slice");
  $lda = $reg->get_adaptor("human","variation","ldfeaturecontainer");
  $vfa = $reg->get_adaptor("human","variation","variationfeature");

  # Get a LDFeatureContainer in a region
  $slice = $sa->fetch_by_region('chromosome', 'X', 1e6, 2e6);

  $ldContainer = $lda->fetch_by_Slice($slice);

  print "Name of the ldContainer is: ", $ldContainer->name();

  # fetch ld featureContainer for a particular variation feature

  $vf = $vfa->fetch_by_dbID(145);

  $ldContainer = $lda->fetch_by_VariationFeature($vf);

  print "Name of the ldContainer: ", $ldContainer->name();


=head1 DESCRIPTION

This adaptor provides database connectivity for LDFeature objects.
LD Features may be retrieved from the Ensembl variation database by
several means using this module.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::LDFeatureContainerAdaptor;

use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Variation::LDFeatureContainer;
use vars qw(@ISA);
use Data::Dumper;

use POSIX;
use FileHandle;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use constant MAX_SNP_DISTANCE => 100_000;

use base qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

our $MAX_SNP_DISTANCE = 100000;
our $BINARY_FILE      = '';
our $TMP_PATH         = '';

sub executable {
  my $self = shift;
  $BINARY_FILE = shift if @_;
  unless( $BINARY_FILE ) {
    my $binary_name = 'calc_genotypes';
    ($BINARY_FILE) = grep {-e $_} map {"$_/calc_genotypes"} split /:/,$ENV{'PATH'};
  }
  return $BINARY_FILE; 
}

sub temp_path {
  my $self = shift;
  $TMP_PATH = shift if @_;
  $TMP_PATH ||= '/tmp'; 
  return $TMP_PATH;
}

=head2 fetch_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice to fetch genes on. Assuming it is always correct (in the top level)
  Arg [2]    : (optional) Bio::EnsEMBL::Variation::Population $population. Population where 
                we want to select the LD information
  Example    : $ldFeatureContainer = $ldfeaturecontainer_adaptor->fetch_by_Slice($slice);
  Description: Overwrites superclass method to add the name of the slice to the LDFeatureContainer.
  Returntype : Bio::EnsEMBL::Variation::LDFeatureContainer
  Exceptions : thrown on bad argument
  Caller     : general
  Status     : Stable

=cut
sub fetch_by_Slice {
    my $self = shift;
    my $slice = shift;
    my $population = shift;

    if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
	throw('Bio::EnsEMBL::Slice arg expected');
    }
    my $sth;
    my $in_str;
    my $siblings = {};
    #when there is no population selected, return LD in the HapMap and PerlEgen populations
    $in_str = $self->_get_LD_populations($siblings);
    #if a population is passed as an argument, select the LD in the region with the population
    if ($population){
	if(!ref($population) || !$population->isa('Bio::EnsEMBL::Variation::Population')) {
	    throw('Bio::EnsEMBL::Variation::Population arg expected');
	}
	my $population_id = $population->dbID;
	$in_str = " = $population_id";
#	if ($in_str =~ /$population_id/){
#	    $in_str = "IN ($population_id)";
#'	}
#	else{
#	    warning("Not possible to calculate LD for a non HapMap or PerlEgen population: $population_id");
#	    return {};
#	}
    }

    if ($in_str eq ''){
	#there is no population, not a human specie or not passed as an argument, return the empy container
	my $t = Bio::EnsEMBL::Variation::LDFeatureContainer->new(
								 '-ldContainer'=> {},
								 '-name' => $slice->name,
								 '-variationFeatures' => {}
								 );
	return $t
    }

    $sth = $self->prepare(qq{SELECT c.individual_id,c.seq_region_id,c.seq_region_start,c.seq_region_end,c.genotypes,ip.population_id
				 FROM compressed_genotype_region c, individual_population ip
				 WHERE  ip.individual_id = c.individual_id
				 AND   ip.population_id $in_str
				 AND   c.seq_region_id = ?
				 AND   c.seq_region_start >= ? and c.seq_region_start <= ?
				 AND   c.seq_region_end >= ?
				 ORDER BY c.seq_region_id, c.seq_region_start},{mysql_use_result => 1});

    $sth->bind_param(1,$slice->get_seq_region_id,SQL_INTEGER);
    $sth->bind_param(2,$slice->start - MAX_SNP_DISTANCE,SQL_INTEGER) if ($slice->start - MAX_SNP_DISTANCE >= 1);
    $sth->bind_param(2,1,SQL_INTEGER) if ($slice->start - MAX_SNP_DISTANCE < 1);
    $sth->bind_param(3,$slice->end,SQL_INTEGER);
    $sth->bind_param(4,($slice->start < 1 ? 1 : $slice->start),SQL_INTEGER);

    $sth->execute();
    
    my $ldFeatureContainer = $self->_objs_from_sth($sth,$slice,$siblings);

    $sth->finish();
    #and store the name of the slice in the Container
    $ldFeatureContainer->name($slice->name());
    return $ldFeatureContainer;
}

=head2 fetch_by_VariationFeature

  Arg [1]    : Bio::EnsEMBL:Variation::VariationFeature $vf
  Arg [2]    : (optional) Bio::EnsEMBL::Variation::Population $pop
  Example    : my $ldFeatureContainer = $ldFetureContainerAdaptor->fetch_by_VariationFeature($vf);  Description: Retrieves LDFeatureContainer for a given variation feature.  If optional population is supplied, values are only returned for that population.
  Returntype : reference to Bio::EnsEMBL::Variation::LDFeatureContainer
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_VariationFeature {
  my $self = shift;
  my $vf  = shift;
  my $pop = shift;

  if(!ref($vf) || !$vf->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
    throw('Bio::EnsEMBL::Variation::VariationFeature arg expected');
  }

  if(!defined($vf->dbID())) {
    throw("VariationFeature arg must have defined dbID");
  }
  
  # cache the position so objs_from_sth picks it up later to filter
  $self->{_vf_pos} = $vf->seq_region_start;
  
  # fetch by slice using expanded feature slice
  my $ldFeatureContainer = $self->fetch_by_Slice($vf->feature_Slice->expand(MAX_SNP_DISTANCE,MAX_SNP_DISTANCE),$pop);
  
  # delete the cached pos
  delete $self->{_vf_pos};
  
  $ldFeatureContainer->name($vf->dbID);
  
  return $ldFeatureContainer;
  
 
}


sub get_populations_by_Slice{
  my $self = shift;
  
  my $population_hash = $self->get_populations_hash_by_Slice(@_);
  return [values(%{$population_hash})];
}

sub get_populations_hash_by_Slice{
  my $self = shift;
  my $slice = shift;
  
  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Bio::EnsEMBL::Slice arg expected');
  }
  
  my $pop_list = $self->_get_LD_populations();
  
  my ($sr, $slice_start, $slice_end) = ($slice->get_seq_region_id, $slice->start, $slice->end);
  
  my %results;
  my $pop_threshold = 20;	# number of individuals required
  my $gen_threshold = 3;	# number of genotypes per individual required

  # just get the population list if it's really too long
  if($slice->length > 10000000) {
	
	my $sth = $self->prepare(qq{SELECT population_id, name FROM population WHERE population_id $pop_list;});
	$sth->execute;
	
	
	%results = map {$_->[0] => $_->[1]} @{$sth->fetchall_arrayref()};
  }

  # do a guesstimate for long slices, otherwise it takes too long
  elsif($slice->length > 5000000) {
	
	my $sth = $self->prepare(qq{
	  SELECT distinct(c.individual_id), p.name
	  FROM compressed_genotype_region c, individual_population ip, population p, individual i
	  WHERE c.individual_id = ip.individual_id
	  AND ip.population_id = p.population_id
	  AND c.individual_id = i.individual_id
	  AND c.seq_region_id = ?
	  AND c.seq_region_start >= ? and c.seq_region_start <= ?
	  AND c.seq_region_end >= ?
	  AND i.father_individual_id is NULL AND i.mother_individual_id is NULL
	  AND (p.population_id $pop_list)
	});
	
	$sth->execute($slice->get_seq_region_id, $slice->start, $slice->end, $slice->start);
	
	my %counts = ();
	
	while(my $row = $sth->fetchrow_arrayref()) {
	  my ($ind_id, $pop_id, $pop_name) = @$row;
	  $results{$pop_id} = $pop_name;
	  $counts{$pop_id}++;
	}
	
	#Delete the populations that doesn't have enough individuals
	delete @results{ grep {$counts{$_} <= $pop_threshold} keys(%counts)};
  }
  
  else {

	my $sth = $self->prepare(qq{
	  SELECT p.population_id, p.name, c.individual_id, c.seq_region_start, c.seq_region_end, c.genotypes 
	  FROM compressed_genotype_region c, individual_population ip, population p, individual i
	  WHERE c.individual_id = ip.individual_id
	  AND ip.population_id = p.population_id
	  AND c.individual_id = i.individual_id
	  AND c.seq_region_id = ?
	  AND c.seq_region_start >= ? and c.seq_region_start <= ?
	  AND c.seq_region_end >= ?
	  AND i.father_individual_id is NULL AND i.mother_individual_id is NULL
	  AND (p.population_id $pop_list)
	});
	
	$sth->execute($sr, $slice_start, $slice_end, $slice_start);
	
	my (%enough, %counts, %sample_pop, %counts_pop);
	
	my $row_count = 0;
	
	while(my $row = $sth->fetchrow_arrayref()) {
	  my ($population_id, $population_name, $individual_id, $start, $end, $genotypes) = @$row;
	  
	  $row_count++;
	  
	  next if $enough{$individual_id};
	  
	  
	  $results{$population_id} = $population_name;
	  $sample_pop{$population_id} = $population_name;
	  
	  # if the row is only partially within the slice
	  if($start < $slice_start || $end > $slice_end) {
		
		my @genotypes = unpack("(www)*", $genotypes);
		my $snp_start = $start;
		
		while( my( $variation_id, $gt_code, $gap ) = splice @genotypes, 0, 3 ) {
		  if(
			($snp_start >= $slice_start) &&
			($snp_start <= $slice_end)
		  ) {		
			$counts{$individual_id}++;
		  }
		  
		  $snp_start += $gap + 1 if defined $gap;
		  last if $snp_start > $slice_end;
		}
	  }
	  
	  # if the row is fully within the slice
	  else {
		$counts{$individual_id} += (((length($genotypes) - 2) / 4) + 1);
	  }
	  
	  $enough{$individual_id} = 1 if $counts{$individual_id} >= $gen_threshold;
	  $counts_pop{$population_id}++ if $counts{$individual_id} >= $gen_threshold;
	}
	
	delete @results{grep {$counts_pop{$_} <= $pop_threshold} keys %counts_pop};
  }
  
  return \%results;
}


#for a given population, gets all individuals that are children (have father or mother)
sub _get_siblings{
    my $self = shift;
    my $population_id = shift;
    my $siblings = shift;

    my $sth_individual = $self->db->dbc->prepare(qq{SELECT i.individual_id
							     FROM individual i, individual_population ip
							     WHERE ip.individual_id = i.individual_id
							     AND ip.population_id = ? 
							     AND i.father_individual_id IS NOT NULL
							     AND i.mother_individual_id IS NOT NULL
							 });
    my ($individual_id);
    $sth_individual->execute($population_id);
    $sth_individual->bind_columns(\$individual_id);
    while ($sth_individual->fetch){
	$siblings->{$population_id.'-'.$individual_id}++; #necessary to have in the key the population, since some individuals are shared between
	                                                   #populations
    }
}

#reads one line from the compress_genotypes table, uncompress the data, and writes it to the different hashes: one containing the number of bases for the variation and the other with the actual genotype information we need to print in the file
sub _store_genotype{
	my $self = shift;
	my $individual_information = shift;
	my $alleles_variation = shift;
	my $individual_id = shift;
	my $seq_region_start = shift;
	my $genotype = shift;
	my $population_id = shift;
	my $slice = shift;
	
	my $snp_start = $seq_region_start;
	my ($slice_start, $slice_end) = ($slice->start, $slice->end);
	
	my @genotypes = unpack("(www)*", $genotype);
	while( my( $variation_id, $gt_code, $gap ) = splice @genotypes, 0, 3 ) {
		my $gt = $self->_genotype_from_code($gt_code);
		
		if(
			defined $gt->[0] &&
			($snp_start >= $slice_start) &&
			($snp_start <= $slice_end) &&
			$gt->[0] =~ /^[ACGT]$/ &&
			(defined $gt->[1] && $gt->[1] =~ /[ACGT]$/)
		) {		
			$alleles_variation->{$snp_start}->{$population_id}->{$gt->[0]}++;
			$alleles_variation->{$snp_start}->{$population_id}->{$gt->[1]}++ if defined $gt->[1];
			
			$individual_information->{$population_id}->{$snp_start}->{$individual_id}->{allele_1} = $gt->[0];
			$individual_information->{$population_id}->{$snp_start}->{$individual_id}->{allele_2} = $gt->[1] if defined $gt->[1];
		}
		
		$snp_start += $gap + 1  if defined $gap;
	}
}

sub _genotype_from_code {
	my $self = shift;
	my $code = shift;
	
	return $self->{'_genotype_codes'}->{$code} if defined($self->{'_genotype_codes'}) && defined($self->{'_genotype_codes'}->{$code});
	
	if(!defined($self->{'_gtca'})) {
		$self->{'_gtca'} = $self->db->get_GenotypeCodeAdaptor;
	}
	my $gtca = $self->{'_gtca'};
	
	my $gtc = $gtca->fetch_by_dbID($code);
	$self->{'_genotype_codes'}->{$code} = $gtc ? $gtc->genotype : undef;
	
	return $self->{'_genotype_codes'}->{$code};
}

#
# Converts the genotype into the required format for the calculation of the pairwise_ld value: AA, Aa or aa
# From the Allele table, will select the alleles and compare to the alleles in the genotype
#

sub _convert_genotype{
    my $self = shift;
    my $alleles_variation = shift; #reference to the hash containing the alleles for the variation present in the genotypes
    my $individual_information = shift; #reference to a hash containing the values to be written to the file
    my @alleles_ordered; #the array will contain the alleles ordered by apparitions in the genotypes (only 2 values possible)
    
    @alleles_ordered = sort({$alleles_variation->{$b} <=> $alleles_variation->{$a}} keys %{$alleles_variation});
    
    #let's convert the allele_1 allele_2 to a genotype in the AA, Aa or aa format, where A corresponds to the major allele and a to the minor
    foreach my $individual_id (keys %{$individual_information}){
	#if both alleles are different, this is the Aa genotype
	if ($individual_information->{$individual_id}{allele_1} ne $individual_information->{$individual_id}{allele_2}){
	    $individual_information->{$individual_id}{genotype} = 'Aa';
	}
	#when they are the same, must find out which is the major
	else{	    
	    if ($alleles_ordered[0] eq $individual_information->{$individual_id}{allele_1}){
		#it is the major allele
		$individual_information->{$individual_id}{genotype} = 'AA';
	    }
	    else{
		$individual_information->{$individual_id}{genotype} = 'aa';
	    }
	    
	}
    }
}

sub _get_LD_populations{
    my $self = shift;
    my $siblings = shift;
    my ($pop_id,$population_name);
#    my $sth = $self->db->dbc->prepare(qq{SELECT s.sample_id, s.name
#				     FROM population p, sample s
#				     WHERE (s.name like 'PERLEGEN:AFD%'
#				     OR s.name like 'CSHL-HAPMAP%')
#				     AND s.sample_id = p.sample_id});
	my $sth = $self->db->dbc->prepare(qq{SELECT population_id, name FROM population WHERE display = 'LD'});

    $sth->execute();
    $sth->bind_columns(\$pop_id,\$population_name);
    #get all the children that we do not want in the genotypes
    my @pops;
    while($sth->fetch){
	if($population_name =~ /CEU|YRI|MEX/){
	    $self->_get_siblings($pop_id,$siblings);
	}
	push @pops, $pop_id;
    }
    
    my $in_str = " IN (" . join(',', @pops). ")";
	
    return $in_str if (defined $pops[0]);
    return '' if (!defined $pops[0]);

}


sub _objs_from_sth {
  my $self = shift;
  my $sth  = shift;
  my $slice = shift;
  my $siblings = shift;

  my ($ld_region_id,$ld_region_start,$ld_region_end,$d_prime,$r2,$sample_count);
  my ($vf_id1,$vf_id2);

  my %feature_container = ();
  my %vf_objects = ();

  #get all Variation Features in Slice
  my $vfa = $self->db->get_VariationFeatureAdaptor();
  my $variations = $vfa->fetch_all_by_Slice($slice); #retrieve all variation features
  #create a hash that maps the position->vf_id
  my %pos_vf = ();
  my $region_Slice = $slice->seq_region_Slice();
  map {$pos_vf{$_->seq_region_start} = $_->transfer($region_Slice)} sort {($a->source_name eq 'dbSNP') <=> ($b->source_name eq 'dbSNP')} @{$variations};

  my %alleles_variation = (); #will contain a record of the alleles in the variation. A will be the major, and a the minor. When more than 2 alleles
  #, the genotypes for that variation will be discarded
  my %individual_information = (); #to store the relation of snps->individuals
  my %regions; #will contain all the regions in the population and the number of genotypes in each one
  my $previous_seq_region_id = 0;

  my %_pop_ids;

  my ($individual_id, $seq_region_id, $seq_region_start,$seq_region_end,$genotypes, $population_id);
  my @cmd = qw(calc_genotypes);

  #open the pipe between processes if the binary file exists in the PATH
  my $bin = $self->executable;
  #open the pipe between processes if the binary file exists in the PATH
  if( ! $bin ) {
    warning("Binary file calc_genotypes not found. Please, read the ensembl-variation/C_code/README.txt file if you want to use LD calculation\n");
    goto OUT;
  }
	
  #my $file= $self->temp_path."/".sprintf( "ld%08x%08x%08x", $$, time, rand( 0x7fffffff) );
  #open IN, ">$file.in";
  $sth->bind_columns(\$individual_id, \$seq_region_id, \$seq_region_start, \$seq_region_end, \$genotypes, \$population_id);
  while($sth->fetch()) {
    #only print genotypes without parents genotyped
    if (!exists $siblings->{$population_id . '-' . $individual_id}){ #necessary to use the population_id
      $self->_store_genotype(\%individual_information,\%alleles_variation, $individual_id, $seq_region_start, $genotypes, $population_id, $slice);
      $previous_seq_region_id = $seq_region_id;
    }
  }
  $sth->finish();
  
  #we have to print the variations
  my (%in_files, %in_file_names);
  
  foreach my $snp_start (sort{$a<=>$b} keys %alleles_variation){
	foreach my $population (keys %{$alleles_variation{$snp_start}}){
	  
	  my $fh;
	  
	  # create file handles in hash
	  if(!defined($in_files{$population})) {
		$fh = FileHandle->new;
		my $f_name = $self->temp_path."/".sprintf( "ld%08x%08x%08x%08x", $population, $$, time, rand( 0x7fffffff));		
		$in_file_names{$population} = $f_name;
		$fh->open(">".$f_name.".in");
		$in_files{$population} = $fh;
	  }
	  
	  else {
		$fh = $in_files{$population};
	  }
	  
	  #if the variation has 2 alleles, print all the genotypes to the file
	  if (keys %{$alleles_variation{$snp_start}{$population}} == 2){
		$self->_convert_genotype($alleles_variation{$snp_start}{$population},$individual_information{$population}{$snp_start});
		foreach my $individual_id (keys %{$individual_information{$population}{$snp_start}}){
		  print $fh join("\t",$previous_seq_region_id,$snp_start, $snp_start,
		  $population, $individual_id,
		  $individual_information{$population}{$snp_start}{$individual_id}{genotype})."\n" || warn $!;
		}
	  }
	}
  }
  
  # close file handles and check file sizes
  foreach my $key(keys %in_files) {
	my $f = $in_files{$key};
	
	my @stats = stat $f;
	$f->close;
	if($stats[7] == 0) {
	  unlink($in_file_names{$key}.'.in');
	  delete $in_file_names{$key};
	  delete $in_files{$key};
	}
  }
  
  # run LD binary
  `$bin <$_\.in >$_\.out` for values %in_file_names;
  
  
  foreach my $file(values %in_file_names) {
	open OUT, "$file.out";
	while(<OUT>){
	  my %ld_values = ();
	  
	  #     936	965891	164284	166818	0.628094	0.999996	120 
	  #get the ouput into the hashes
	  chomp;
	  
	  ($population_id,$ld_region_id,$ld_region_start,$ld_region_end,$r2,$d_prime,$sample_count) = split /\s/;
	  
	  # skip entries unrelated to selected vf if doing fetch_all_by_VariationFeature
	  if(defined($self->{_vf_pos})) {
		next unless $ld_region_start == $self->{_vf_pos} || $ld_region_end == $self->{_vf_pos};
	  }
	  
	  $ld_values{'d_prime'} = $d_prime;
	  $ld_values{'r2'} = $r2;
	  $ld_values{'sample_count'} = $sample_count;
	  
	  if (!defined $pos_vf{$ld_region_start} || !defined $pos_vf{$ld_region_end}){
		next; #problem to fix in the compressed genotype table: some of the positions seem to be wrong
	  }
	  $vf_id1 = $pos_vf{$ld_region_start}->dbID();
	  $vf_id2 = $pos_vf{$ld_region_end}->dbID();
	  
	  $feature_container{$vf_id1 . '-' . $vf_id2}->{$population_id} = \%ld_values;
	  $vf_objects{$vf_id1} = $pos_vf{$ld_region_start};
	  $vf_objects{$vf_id2} = $pos_vf{$ld_region_end};
	  
	  $_pop_ids{$population_id} = 1;	  
	}
	close OUT || die "Could not close filehandle: $!\n";
	
	unlink( "$file.in" );
	unlink( "$file.out" );
  }
  
OUT:
  my $t = Bio::EnsEMBL::Variation::LDFeatureContainer->new(
	'-ldContainer'=> \%feature_container,
	'-name' => '',
	'-variationFeatures' => \%vf_objects
  );

  $t->{'_pop_ids'} =\%_pop_ids;

  return $t;      
}

1;
