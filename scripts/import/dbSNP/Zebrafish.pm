use strict;
use warnings;
#object that contains the specific methods to dump data from dbSNP. For zebra-fish, there are no coordinates, so Yuan maps them to our data, and 
#adds the file into the variation_feature table
package dbSNP::Zebrafish;

use dbSNP::GenericContig;
use vars qw(@ISA);
use ImportUtils qw(debug load);

@ISA = ('dbSNP::GenericContig');


sub variation_feature{
    my $self = shift;
    
    my %rec;
    debug("Dumping Variation");
    
    my $sth = $self->{'dbVariation'}->prepare (qq{SELECT variation_id, name
						    FROM variation});
    $sth->execute();
    
    while(my ($variation_id, $name) = $sth->fetchrow_array()) {
	$rec{$name} = $variation_id;
    }
    
    $sth->finish();
    
    open (IN, $self->{'alldiff'});
    open (FH, ">" . $self->{'tmpdir'} . "/" . $self->{'tmpfile'} );
    
    while (<IN>) {
	s/^MORE_HITS//;
	my ($ref_id, $slice_name, $ver, $start, $end, $exact, $snp_type, $strand, $score, $ratio) =split;
	next if $ratio <0.5;
	
	##make start < end if it is a between type
	if ($exact =~ /between/i) {
	    $end = $start-1;
	}
	
	my ($coord_sys,$assembly,$seq_region_name,$seq_region_start,$seq_region_end,$version) = split /\:/, $slice_name
	    if $slice_name =~ /\:/;
	
	my $sth = $self->{'dbCore'}->dbc->prepare (qq{SELECT seq_region_id from seq_region where name = ?});
	$sth->execute("$seq_region_name");
	
	my $seq_region_id = $sth->fetchrow_array();
	
	my $new_seq_region_start = $seq_region_start + $start -1 if ($seq_region_start);
	my $new_seq_region_end = $seq_region_start + $end -1 if ($seq_region_start);
	
	$ref_id = "rs$ref_id";
	print FH "$seq_region_id\t$new_seq_region_start\t$new_seq_region_end\t$strand\t$rec{$ref_id}\t$ref_id\n";
    }
    
    $sth->finish();
    
    close IN;
    close FH;
    
    load($self->{'dbVariation'}, "variation_feature","seq_region_id","seq_region_start","seq_region_end",
	 "seq_region_strand","variation_id","variation_name");
}

1;
