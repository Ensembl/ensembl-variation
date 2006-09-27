#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor

=head1 SYNOPSIS

  $db = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(...);

  $va = $db->get_VariationAdaptor();
  $vga = $db->get_VariationGroupAdaptor();
  $pa = $db->get_PopulationAdaptor();

  # Get a Variation by its internal identifier
  $var = $va->fetch_by_dbID(145);

  # fetch a variation by its name
  $var = $va->fetch_by_name('rs100');


  # fetch all variations from a population
  $pop = $pa->fetch_by_name('PACIFIC');
  @vars = {$va->fetch_all_by_Population($pop)};

=head1 DESCRIPTION

This adaptor provides database connectivity for Variation objects.
Variations (SNPs, etc.) may be retrieved from the Ensembl variation database by
several means using this module.

=head1 AUTHOR - Graham McVicker

=head1 CONTACT

Post questions to the Ensembl development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Variation::Allele;

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');


=head2 fetch_by_dbID

  Arg [1]    : int $dbID
  Example    : $var = $var_adaptor->fetch_by_dbID(5526);
  Description: Retrieves a Variation object via its internal identifier.
               If no such variation exists undef is returned.
  Returntype : Bio::EnsEMBL::Variation::Variation
  Exceptions : throw if dbID arg is not defined
  Caller     : general, IndividualAdaptor
  Status     : Stable

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  throw('dbID argument expected') if(!defined($dbID));

  my $sth = $self->prepare
    (q{SELECT v.variation_id, v.name, v.validation_status, s1.name, v.ancestral_allele,
              a.allele_id, a.allele, a.frequency, a.sample_id, vs.moltype,
              vs.name, s2.name
       FROM   (variation v, source s1, allele a)
                LEFT JOIN variation_synonym vs on v.variation_id = vs.variation_id 
              LEFT JOIN source s2 on  vs.source_id = s2.source_id
       WHERE  v.variation_id = a.variation_id
       AND    v.source_id = s1.source_id
       AND    v.variation_id = ?});
  $sth->bind_param(1,$dbID,SQL_INTEGER);
  $sth->execute();

  my $result = $self->_objs_from_sth($sth);
  $sth->finish();
  return undef if(!@$result);

  return $result->[0];
}



=head2 fetch_by_name

  Arg [1]    : string $name
  Example    : $var = $var_adaptor->fetch_by_name('rs1453','dbSNP');
  Description: Retrieves a population object via its name
  Returntype : Bio::EnsEMBL::Variation::Variation
  Exceptions : throw if name argument is not defined
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_name {
  my $self = shift;
  my $name = shift;
  my $source = shift;

  throw('name argument expected') if(!defined($name));
  my $extra_sql = "";

  if ( defined $source ) {
    $extra_sql = qq(AND    s1.name = ?);
  }

  my $sth = $self->prepare
    (qq{SELECT v.variation_id, v.name, v.validation_status, s1.name, v.ancestral_allele,
              a.allele_id, a.allele, a.frequency, a.sample_id, vs.moltype,
              vs.name, s2.name
#       FROM   variation v, source s1, source s2, allele a, variation_synonym vs
	  FROM   (variation v, source s1, allele a) 
	      LEFT JOIN variation_synonym vs on v.variation_id = vs.variation_id 
              LEFT JOIN source s2 on  vs.source_id = s2.source_id
       WHERE  v.variation_id = a.variation_id
#       AND    v.variation_id = vs.variation_id
       AND    v.source_id = s1.source_id
#       AND    vs.source_id = s2.source_id
       AND    v.name = ?
       $extra_sql  
       ORDER BY a.allele_id});

  $sth->bind_param(1,$name,SQL_VARCHAR);
  $sth->bind_param(2,$source,SQL_VARCHAR) if defined $source;
  $sth->execute();

  my $result = $self->_objs_from_sth($sth);
  $sth->finish();

  if(!@$result) {
    # try again if nothing found, but check synonym table instead
    $sth = $self->prepare
      (qq{SELECT v.variation_id, v.name, v.validation_status, s1.name, v.ancestral_allele,
                a.allele_id, a.allele, a.frequency, a.sample_id, vs1.moltype,
                vs2.name, s2.name
         FROM   variation v, source s1, source s2, allele a,
                variation_synonym vs1, variation_synonym vs2
         WHERE  v.variation_id = a.variation_id
         AND    v.variation_id = vs1.variation_id
         AND    v.variation_id = vs2.variation_id
         AND    v.source_id = s1.source_id
         AND    vs2.source_id = s2.source_id
         AND    vs1.name = ?
                $extra_sql
         ORDER BY a.allele_id});
    $sth->bind_param(1,$name,SQL_VARCHAR);
    $sth->bind_param(2,$source,SQL_VARCHAR) if defined $source;
    $sth->execute();
    $result = $self->_objs_from_sth($sth);

    return undef if(!@$result);

    $sth->finish();
  }
  return $result->[0];
}



=head2 fetch_all_by_dbID_list

  Arg [1]    : reference to list of ints $list
  Example    : @vars = @{$va->fetch_all_by_dbID_list([124, 56, 90])};
  Description: Retrieves a set of variations via their internal identifiers.
               This is faster than repeatedly calling fetch_by_dbID if there
               are a large number of variations to retrieve
  Returntype : reference to list of Bio::EnsEMBL::Variation::Variation objects
  Exceptions : throw on bad argument
  Caller     : general, IndividualGenotypeAdaptor, PopulationGenotypeAdaptor
  Status     : At Risk

=cut


sub fetch_all_by_dbID_list {
  my $self = shift;
  my $list = shift;

  if(!defined($list) || ref($list) ne 'ARRAY') {
    throw("list reference argument is required");
  }

  return [] if(!@$list);

  my @out;

  # mysql is faster and we ensure that we do not exceed the max query size by
  # splitting large queries into smaller queries of 200 ids
  my $max = 200;

  while(@$list) {
    my @ids = (@$list > $max) ? splice(@$list, 0, $max) : splice(@$list, 0);

    my $id_str = (@ids > 1)  ? " IN (".join(',',@ids).")"   :   ' = '.$ids[0];

    my $sth = $self->prepare
      (qq{SELECT v.variation_id, v.name, v.validation_status, s1.name, v.ancestral_allele,
                 a.allele_id, a.allele, a.frequency, a.sample_id, vs.moltype,
                 vs.name, s2.name
	     FROM   (variation v, source s1, allele a)
	      LEFT JOIN variation_synonym vs on v.variation_id = vs.variation_id 
              LEFT JOIN source s2 on  vs.source_id = s2.source_id
          WHERE  v.variation_id = a.variation_id
          AND    v.source_id = s1.source_id
          AND    v.variation_id $id_str});
    $sth->execute();

    my $result = $self->_objs_from_sth($sth);

    $sth->finish();

    push @out, @$result if(@$result);
  }

  return \@out;
}


=head2 get_all_sources

  Args        : none
  Example     : $sources = $va->get_all_sources();
  Description : Retrieves from the database all sources in the Source table
  ReturnType  : array ref of string
  Exceptions  : none
  Caller      : web
  Status      : At Risk

=cut

sub get_all_sources{
    my $self = shift;
    my @sources;

    my $source_name;
    my $sth = $self->prepare(qq{SELECT name from source
				});
    $sth->execute();
    $sth->bind_columns(\$source_name);
    
    while ($sth->fetch()){
	push @sources, $source_name
    }
    $sth->finish();

    return \@sources;
}


=head2 get_default_source

  Args      : none
  Example     : $default_source = $va->get_default_source();
  Description : Retrieves from the database the default source used for display purposes
  ReturnType  : string
  Exceptions  : none
  Caller      : web
  Status      : At Risk

=cut

sub get_default_source{
    my $self = shift;

    my $source_name;
    my $sth = $self->prepare(qq{SELECT meta_value from meta where meta_key = ?
				});
    $sth->bind_param(1,'source.default_source',SQL_VARCHAR);
    $sth->execute();
    $sth->bind_columns(\$source_name);
    $sth->fetch();
    $sth->finish();

    return $source_name;
}


=head2 get_source_version

  Arg[1]      : string $name
  Example     : $version = $va->get_source_version('dbSNP');
  Description : Retrieves from the database the version for the source given as an argument
  ReturnType  : int
  Exceptions  : none
  Caller      : general
  Status      : At Risk

=cut

sub get_source_version{
    my $self = shift;
    my $name = shift;
    my $version;
    my $sth = $self->prepare(qq{SELECT version from source where name = ?
				});

    $sth->bind_param(1,$name,SQL_VARCHAR);

    $sth->execute();
    $sth->bind_columns(\$version);
    $sth->fetch();
    $sth->finish();

    return $version;
}

=head2 get_flanking_sequence

    Arg[1]      : int $variationID
    Example     : $flankinq_sequence = $va->get_flanking_sequence('652');
    Description : Retrieves from the database the appropriate flanking sequence (five,three) for the variation. If the flanking sequence is not in
                  the Flankinq_sequence table, access the core database with the coordinates
    ReturnType  : reference to a list containing (three_flank,five_flank)
    Exceptions  : throw when not possible to obtain sequence
    Caller      : general, Variation
    Status      : Stable

=cut

sub get_flanking_sequence{
  my $self = shift;
  my $variationID = shift;

  my $flanking_sequence; #reference to an array for the three_prime and five_prime seqs
  my ($seq_region_id, $seq_region_strand, $up_seq, $down_seq, $up_seq_region_start, $up_seq_region_end, $down_seq_region_start, $down_seq_region_end);
  
  my $sth = $self->prepare(qq{
			      SELECT seq_region_id, seq_region_strand, up_seq, 
			      down_seq, up_seq_region_start, up_seq_region_end, 
			      down_seq_region_start, down_seq_region_end
			      FROM flanking_sequence
			      WHERE variation_id = ?
			     });

  $sth->bind_param(1,$variationID,SQL_INTEGER);
  $sth->execute(); #retrieve the flank from the variation database
  $sth->bind_columns(\($seq_region_id, $seq_region_strand, $up_seq, $down_seq, $up_seq_region_start, $up_seq_region_end, $down_seq_region_start, $down_seq_region_end));
$sth->fetch();
$sth->finish();

if (!defined $down_seq){
  if( $seq_region_id){
    $down_seq = $self->_get_flank_from_core($seq_region_id, 
                                            $down_seq_region_start, 
                                            $down_seq_region_end, 
                                            $seq_region_strand);
  } else {
    warn( "*****[ERROR]: No seq_region_id for SNP with dbID: $variationID. ".
          "Cannot retrieve flanking region******\n" );    
  }
}
if (!defined $up_seq){
  if( $seq_region_id){
    $up_seq = $self->_get_flank_from_core($seq_region_id, 
                                          $up_seq_region_start, 
                                          $up_seq_region_end, 
                                          $seq_region_strand);
  } else {
    warn( "*****[ERROR]: No seq_region_id for SNP with dbID: $variationID. ".
          "Cannot retrieve flanking region******\n" );
  }
}

push @{$flanking_sequence},$down_seq,$up_seq; #add to the array the 3 and 5 prime sequences

return $flanking_sequence;
}

=head2 fetch_all_by_Population

  Arg [1]    : Bio::EnsEMBL::Variation::Population
  Example    : $pop = $pop_adaptor->fetch_by_dbID(1345);
               @vars = @{$va_adaptor->fetch_all_by_Population($pop)};
  Description: Retrieves all variations which are stored for a specified
               population.
  Returntype : listref of Bio::EnsEMBL::Variation::Variation
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_Population {
  my $self = shift;
  my $pop = shift;

  if(!ref($pop) || !$pop->isa('Bio::EnsEMBL::Variation::Population')) {
    throw('Bio::EnsEMBL::Variation::Population argument expected');
  }

  if(!defined($pop->dbID())) {
    warning("Cannot retrieve genotypes for population without set dbID");
    return [];
  }

  my $sth = $self->prepare
    (q{SELECT v.variation_id, v.name, v.validation_status, s1.name, v.ancestral_allele,
              a.allele_id, a.allele, a.frequency, a.sample_id, vs.moltype,
              vs.name, s2.name
	    FROM   (variation v, source s1, allele a)
	      LEFT JOIN variation_synonym vs on v.variation_id = vs.variation_id 
              LEFT JOIN source s2 on  vs.source_id = s2.source_id
	      WHERE  v.variation_id = a.variation_id
	      AND    v.source_id = s1.source_id
	      AND    a.sample_id = ?});
  $sth->bind_param(1,$pop->dbID,SQL_INTEGER);
  $sth->execute();

  my $results = $self->_objs_from_sth($sth);
  $sth->finish();

  return $results;
}


sub _get_flank_from_core{
    my $self = shift;
    my $seq_region_id = shift;
    my $seq_region_start = shift;
    my $seq_region_end = shift;
    my $seq_region_strand = shift;

    my $flanking_sequence;
    if (defined $self->db()->dnadb()){
	my $slice_adaptor = $self->db()->dnadb()->get_SliceAdaptor();
	my $slice = $slice_adaptor->fetch_by_seq_region_id($seq_region_id);
	if (!$slice){
	    throw("Not possible to obtain slice for seq_region_id \"$seq_region_id\"\n");
	}
	my $flank = $slice->subseq($seq_region_start,$seq_region_end,$seq_region_strand);
	return $slice->subseq($seq_region_start,$seq_region_end,$seq_region_strand);
    }
    return '';
}


sub _objs_from_sth {
  my $self = shift;
  my $sth = shift;

  my ($var_id, $name, $vstatus, $source, $ancestral_allele, $allele_id, $allele, $allele_freq,
      $allele_sample_id, $moltype, $syn_name, $syn_source,
      $cur_allele_id, $cur_var, $cur_var_id);

  $sth->bind_columns(\$var_id, \$name, \$vstatus, \$source, \$ancestral_allele, \$allele_id,
                     \$allele, \$allele_freq, \$allele_sample_id, \$moltype, \$syn_name,
                     \$syn_source);

  my @vars;

  my %seen_syns;
  my %seen_pops;

  my $pa = $self->db()->get_PopulationAdaptor();

  while($sth->fetch()) {
    if(!defined($cur_var) || $cur_var_id != $var_id) {
	$vstatus = 0 if (!defined $vstatus);
      my @states = split(',',$vstatus);
      $cur_var = Bio::EnsEMBL::Variation::Variation->new
        (-dbID   => $var_id,
         -ADAPTOR => $self,
         -NAME   => $name,
         -SOURCE => $source,
	 -ANCESTRAL_ALLELE => $ancestral_allele,
	 -MOLTYPE => $moltype,
         -VALIDATION_STATES => \@states);
      push @vars, $cur_var;
      $cur_var_id = $var_id;
    }

    if(!defined($cur_allele_id) || $cur_allele_id != $allele_id) {
      my $pop;
      if($allele_sample_id) {
        $pop = $seen_pops{$allele_sample_id} ||=
          $pa->fetch_by_dbID($allele_sample_id);
    }
      my $allele = Bio::EnsEMBL::Variation::Allele->new
        (-dbID      => $allele_id,
         -ALLELE    => $allele,
         -FREQUENCY => $allele_freq,
         -POPULATION => $pop);
      $cur_var->add_Allele($allele);

      $cur_allele_id = $allele_id;
    }

    if(defined ($syn_source) && !$seen_syns{"$syn_source:$syn_name"}) {
      $seen_syns{"$syn_source:$syn_name"} = 1;
      $cur_var->add_synonym($syn_source, $syn_name);
    }
}
  return \@vars;
}




1;
