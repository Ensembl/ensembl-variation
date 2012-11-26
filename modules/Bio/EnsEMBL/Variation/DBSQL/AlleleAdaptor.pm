=head1 LICENSE

 Copyright (c) 1999-2012 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::AlleleAdaptor
#
# Copyright (c) 2011 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::AlleleAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $va = $reg->get_adaptor("human","variation","allele");



=head1 DESCRIPTION

This adaptor provides database connectivity for Allele objects.
Alleles may be retrieved from the Ensembl variation database by
several means using this module.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::AlleleAdaptor;

use Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);

use Bio::EnsEMBL::Variation::Allele;
use Bio::EnsEMBL::Utils::Iterator;

use Scalar::Util qw(weaken);

use DBI qw(:sql_types);

our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor');

my $DEFAULT_ITERATOR_CACHE_SIZE = 1000;


sub store {
	my ($self, $allele) = @_;
	
	my $dbh = $self->dbc->db_handle;
	
	# look up or store allele code
	my $allele_code = $self->_allele_code($allele->allele);
	
	my $sth = $dbh->prepare_cached(q{
		INSERT DELAYED INTO allele (
			variation_id,
			subsnp_id,
			allele_code_id,
			sample_id,
			frequency,
			count			
		) VALUES (?,?,?,?,?,?)
	});
	
	$sth->execute(
		$allele->{_variation_id} || $allele->variation->dbID,
		$allele->{subsnp},
		$allele_code,
		$allele->population ? $allele->population->dbID : undef,
		$allele->frequency,
		$allele->count
	);
	
	$sth->finish;
}

sub store_multiple {
	my ($self, $alleles) = @_;
	
	my $dbh = $self->dbc->db_handle;
	
	my $q_string = join ",", map {'(?,?,?,?,?,?)'} @$alleles;
	
	my @args = map {
		$_->{_variation_id} || $_->variation->dbID,
		$_->{subsnp},
		$self->_allele_code($_->allele),
		$_->population ? $_->population->dbID : undef,
		$_->frequency,
		$_->count
	} @$alleles;
	
	my $sth = $dbh->prepare_cached(qq{
		INSERT INTO allele (
			variation_id,
			subsnp_id,
			allele_code_id,
			sample_id,
			frequency,
			count			
		) VALUES $q_string
	});
	
	$sth->execute(@args);
	
	$sth->finish;
}

sub store_to_file_handle {
	my ($self, $allele, $file_handle) = @_;
	
	my $dbh = $self->dbc->db_handle;
	
	print $file_handle join("\t",
		$allele->{_variation_id} || $allele->variation->dbID || '\N',
		$allele->{subsnp} || '\N',
		$self->_allele_code($allele->allele),
		$allele->population ? $allele->population->dbID : '\N',
		defined($allele->frequency) ? $allele->frequency :  '\N',
		defined($allele->count) ? $allele->count : '\N',
	)."\n";
}

=head2 fetch_all

  Description: fetch_all should not be used for Alleles.
  Exceptions : thrown on invocation
  Status     : At risk

=cut

sub fetch_all {
    my $self = shift;

    throw("fetch_all cannot be used for Allele objects");
}

=head2 fetch_all_by_subsnp_id

  Arg [1]    : string $subsnp_id
  Example    : $alleles = $allele_adaptor->fetch_all_by_subsnp_id('ss123');
  Description: Retrieves all allele objects via a component subsnp ID
  Returntype : listref of Bio::EnsEMBL::Variation::Allele objects
  Exceptions : throw if name argument is not defined
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_subsnp_id {
    my $self = shift;
    my $name = shift;
  
    $name =~ s/^ss//gi;

    throw('name argument expected') if(!defined($name));

    # Add the constraint on the subsnp_id column and pass to generic_fetch
    my $constraint = qq{ a.subsnp_id = $name };

    return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_Variation

  Arg [1]    : Bio::EnsEMBL::Variation::Variation
  Arg [2]    : Bio::EnsEMBL::Variation::Population (optional)
  Example    : @alleles = @{$allele_adaptor->fetch_all_by_Variation($var)};
  Description: Retrieves all alleles which are associated with a specified
               variation. If the optional population argument is specified, only retrieve
               alleles for that population.
  Returntype : listref of Bio::EnsEMBL::Variation::Allele
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_Variation {
    my $self = shift;
    my $variation = shift;
    my $population = shift;
    
    # Make sure that we are passed a Variation object
    assert_ref($variation,'Bio::EnsEMBL::Variation::Variation');
    
    # If we got a population argument, make sure that it is a Population object
    assert_ref($population,'Bio::EnsEMBL::Variation::Population') if (defined($population));
    
    # Add a constraint on the variation_id column and pass to generic fetch
    my $variation_id = $variation->dbID();
    my $constraint = qq{ a.variation_id = $variation_id };
    
    # If required, add a constraint on the sample id
    if (defined($population)) {
        my $sample_id = $population->dbID();
        $constraint .= qq{ AND a.sample_id = $sample_id };
    }
    
    # Add the constraint for failed alleles
    $constraint .= " AND " . $self->db->_exclude_failed_alleles_constraint();
  
    my $alleles = $self->generic_fetch($constraint);
    
    # If a population was specified, attach the population to the allele object
    map {$_->population($population)} @{$alleles} if (defined($population));
	
	# add freqs from genotypes for human (1KG data)
	push @$alleles, @{$self->_fetch_all_by_Variation_from_Genotypes($variation, $population)} if $self->db->species =~ /homo_sapiens/i;
	
    # Return the alleles
    return $alleles;
}

sub _fetch_all_by_Variation_from_Genotypes {
	my $self = shift;
	my $variation = shift;
    my $population = shift;
    
    # Make sure that we are passed a Variation object
    assert_ref($variation,'Bio::EnsEMBL::Variation::Variation');
    
    # If we got a population argument, make sure that it is a Population object
    assert_ref($population,'Bio::EnsEMBL::Variation::Population') if (defined($population));
	
	# fetch all genotypes
	my $genotypes = $variation->get_all_IndividualGenotypes();
	
	return [] unless scalar @$genotypes;
	
	# get populations for individuals
	my (@pop_list, %pop_hash);
	
	if(defined($population)) {
		@pop_list = ($population);
		map {$pop_hash{$population->dbID}{$_->dbID} = 1} @{$population->get_all_Individuals};
	}
	else {
		my $pa = $self->db->get_PopulationAdaptor();
		%pop_hash = %{$pa->_get_individual_population_hash([map {$_->individual->dbID} @$genotypes])};
		return [] unless %pop_hash;
		
		@pop_list = @{$pa->fetch_all_by_dbID_list([keys %pop_hash])};
	}
	
	return [] unless @pop_list and %pop_hash;
	
	my %ss_list = map {$_->subsnp || '' => 1} @$genotypes;
	
	my @objs;
	
	foreach my $pop(@pop_list) {
		
		# HACK: only include 1KG phase 1 pops for now
		next unless $pop->name =~ /^1000GENOMES:phase_1_/;
		
		foreach my $ss(keys %ss_list) {
			
			my (%counts, $total, @freqs);
			map {$counts{$_}++}
				map {@{$_->{genotype}}}
				grep {$pop_hash{$pop->dbID}{$_->individual->dbID}}
				grep {$_->subsnp || '' eq $ss}
				@$genotypes;
			
			next unless %counts;
			
			my @alleles = keys %counts;
			$total += $_ for values %counts;
			next unless $total;
			
			@freqs = map {defined($counts{$_}) ? ($counts{$_} / $total) : 0} @alleles;
			
			for my $i(0..$#alleles) {
				push @objs, Bio::EnsEMBL::Variation::Allele->new_fast({
					allele     => $alleles[$i],
					count      => scalar keys %counts ? ($counts{$alleles[$i]} || 0) : undef,
					frequency  => @freqs ? $freqs[$i] : undef,
					population => $pop,
					variation  => $variation,
					adaptor    => $self,
					subsnp     => $ss eq '' ? undef : $ss,
				});
				
				weaken($objs[-1]->{'variation'});
			}
		}
	}
	
	return \@objs;
}

=head2 get_all_failed_descriptions

  Arg[1]      : Bio::EnsEMBL::Variation::Allele
	               The allele object to get the failed descriptions for
  Example     : 
                my $failed_descriptions = $adaptor->get_all_failed_descriptions($allele);
                if (scalar(@{$failed_descriptions})) {
		          print "The allele '" . $allele->allele() . "' has been flagged as failed because '" . join("' and '",@{$failed_descriptions}) . "'\n";
                }
		
  Description : Gets the unique descriptions for the reasons why the supplied allele has failed.
  ReturnType  : reference to a list of strings
  Exceptions  : thrown on incorrect argument
  Caller      : general
  Status      : At Risk

=cut

sub get_all_failed_descriptions {
    my $self = shift;
    my $allele = shift;
    
    # Call the internal get method without any constraints
    my $description = $self->_internal_get_failed_descriptions($allele) || [];
    
    return $description;
}

=head2 get_subsnp_handle

  Arg[1]      : Bio::EnsEMBL::Variation::Allele
	               The allele object to get the observation submitter subsnp handle for
  Arg[2]      : Bio::EnsEMBL::Variation::Population (optional)
                       The population object to get the frequency submitter subsnp handle for                           
  Example     : 
                my $handle = $adaptor->get_subsnp_handle($allele);
		print "The allele '" . $allele->allele() . "' of subsnp 'ss" . $allele->subsnp() . "' was submitted by '$handle'\n";

                my $handle = $adaptor->get_subsnp_handle($allele, $population);
		print "Frequency data for subsnp 'ss" . $allele->subsnp() . "' in population " . $population->name() . " was submitted by '$handle'\n";


		
  Description : Gets the submitter handle for the specified allele
                The initial observation and later allele frequencies may be submitted for the same ss record by different groups
                supplying a population returns the submitter handle for data obtained within this population  
  ReturnType  : string
  Exceptions  : thrown on incorrect argument
  Caller      : general
  Status      : At Risk

=cut

sub get_subsnp_handle {
    my $self = shift;
    my $allele = shift;
    my $population = shift ; ## return submitter for specific allele frequency

    # Assert that the object passed is an Allele
    assert_ref($allele,'Bio::EnsEMBL::Variation::Allele');
    
    # If we got a population argument, make sure that it is a Population object
    assert_ref($population,'Bio::EnsEMBL::Variation::Population') if (defined($population));


    # Get the subsnp id and get rid of any 'ss' prefix
    my $ssid = $allele->subsnp() || "";
    $ssid =~ s/^ss//;
    
    my ($stmt, $sth);
    if(defined $population ){
       $stmt = qq{
        SELECT
            submitter_handle.handle
        FROM
            submitter_handle, allele
        WHERE
            allele.subsnp_id = ?
        AND allele.frequency_submitter_handle = submitter_handle.handle_id
        AND allele.sample_id = ?
        LIMIT 1
        };
       $sth = $self->prepare($stmt);
       $sth->execute($ssid, $population->dbID());       
       
    }
    else{
        $stmt = qq{
          SELECT
              handle
          FROM
              subsnp_handle
          WHERE
              subsnp_id = ?
          LIMIT 1
        };


    $sth = $self->prepare($stmt);
    $sth->execute($ssid);

    }
    my $handle;
    $sth->bind_columns(\$handle);
    $sth->fetch();
    return $handle;
}


# API-internal method for getting failed descriptions for an Allele
sub _internal_get_failed_descriptions {
    my $self = shift;
    my $allele = shift;
    my $constraint = shift;
    
    # Assert that the object passed is an Allele
    assert_ref($allele,'Bio::EnsEMBL::Variation::Allele');
    
    my $stmt = qq{
        SELECT DISTINCT
            fd.description
        FROM
            failed_allele fa JOIN
            failed_description fd ON (
                fd.failed_description_id = fa.failed_description_id
            )
        WHERE
            fa.allele_id = ?
    };
    $stmt .= qq{ AND $constraint } if (defined($constraint));
    
    my $sth = $self->prepare($stmt);
    $sth->execute($allele->dbID());
    my @descriptions;
    my $description;
    $sth->bind_columns(\$description);
    while ($sth->fetch()) {
        push(@descriptions,$description);
    }
    return \@descriptions;
}

sub _objs_from_sth {
    my $self = shift;
    my $sth = shift;

    my ($allele_id, $variation_id, $subsnp_id, $allele, $frequency, $sample_id, $count, $last_allele_id);
    my @alleles;
    
    $sth->bind_columns(\$allele_id, \$variation_id, \$subsnp_id, \$allele, \$frequency, \$sample_id, \$count);
    
    while($sth->fetch()) {
    
        # The left join to failed allele can create duplicate rows, so check that we've got a new Allele before creating the object
        unless (defined($last_allele_id) && $last_allele_id == $allele_id) {
        
            my $obj = Bio::EnsEMBL::Variation::Allele->new(
                -dbID           => $allele_id,
                -VARIATION_ID   => $variation_id,
                -SUBSNP         => $subsnp_id,
                -ALLELE         => $allele,
                -FREQUENCY      => $frequency,
                -POPULATION_ID  => $sample_id,
                -COUNT          => $count,
                -ADAPTOR        => $self
            );
              
            push(@alleles,$obj);
            $last_allele_id = $allele_id;
        }
        
    }
    
    return \@alleles;
}

# method used by superclass to construct SQL
sub _tables { 
    my $self = shift;
    
    my @tables = (
        ['allele', 'a'], ['allele_code', 'ac']
    );
    
	# If we are excluding failed_alleles, add that table
	push(@tables,['failed_allele', 'fa']) unless ($self->db->include_failed_variations());
	
	return @tables;
}

# Add a left join to the failed_variation table
sub _left_join { 
    my $self = shift;
    
    # If we are including failed variations, skip the left join
    return () if ($self->db->include_failed_variations());
    return ([ 'failed_allele', 'fa.allele_id = a.allele_id']); 
}

sub _columns {
  return qw( a.allele_id a.variation_id a.subsnp_id ac.allele a.frequency a.sample_id a.count );
}

sub _write_columns {
	return qw(variation_id subsnp_id allele_code_id sample_id frequency count);
}

sub _default_where_clause  {
  return 'a.allele_code_id = ac.allele_code_id';
}

sub _cache_allele_codes {
	my $self = shift;
	
	my $dbh = $self->dbc->db_handle;
	
	my $sth = $dbh->prepare(qq{SELECT allele_code_id, allele FROM allele_code});
	$sth->execute;
	
	my ($code, $allele);
	$sth->bind_columns(\$code, \$allele);
	my %allele_codes;
	$allele_codes{defined($allele) ? $allele : ''} = $code while $sth->fetch;
	$sth->finish();
	
	$self->db->{_allele_codes} = \%allele_codes;
	
	return $self->db->{_allele_codes};
}

sub _allele_code {
	my ($self, $allele) = @_;
	
	# check if cache is loaded
	my $just_loaded = 0;
	
	if(!exists($self->db->{_allele_codes})) {
		$self->_cache_allele_codes;
		$just_loaded = 1;
	}
	
	if(!exists($self->db->{_allele_codes}->{$allele})) {
		
		# check another process hasn't created it by reloading the cache
		$self->_cache_allele_codes unless $just_loaded;
		
		# if it still doesn't exist
		if(!exists($self->db->{_allele_codes}->{$allele})) {
			my $dbh = $self->dbc->db_handle;
			
			my $sth = $dbh->prepare(q{
				INSERT INTO allele_code (
					allele
				)
				VALUES (?)
			});
			eval {
				$sth->execute($allele);
			};
			$sth->finish;
			
			my $allele_code;
			
			# insert failed, another process did it maybe?
			if($@) {
				my $sth2 = $dbh->prepare(qq{
					SELECT allele_code_id
					FROM allele_code
					WHERE allele = ?
				});
				$sth2->execute($allele);
				$sth2->bind_columns(\$allele_code);
				$sth2->fetch();
				
				throw("ERROR: Failed to insert allele '$allele' into allele_code") unless defined($allele_code);
				
				$sth2->finish();
			}
			else {
				$allele_code = $dbh->last_insert_id(undef, undef, 'allele_code', 'allele_code_id');
			}
			
			$self->db->{_allele_codes}->{$allele} = $allele_code;
		}
	}
	
	return $self->db->{_allele_codes}->{$allele};
}
1;
