#! perl -w

# Script to dump out a table of variation sets that can be used in the documentation

use strict;
use warnings;
use Data::Dumper;
use Bio::EnsEMBL::Variation::Utils::Config qw(@OVERLAP_CONSEQUENCES);



my %colour = (
	'Essential splice site'  => 'coral',
	'Stop gained'            => '#F00',
	'Stop lost'              => '#F00', 
	'Complex in/del'         => 'mediumspringgreen',
	'Frameshift coding'      => 'hotpink', 
	'Non synonymous coding'  => 'gold',
	'Splice site'            => 'coral',
	'Partial codon'          => 'magenta',
	'Synonymous coding'      => '#76EE00',
	'Coding unknown'         => '#458B00',
	'Within mature miRNA'    => '#458B00',
	'5prime UTR'             => '#7AC5CD',
	'3prime UTR'             => '#7AC5CD',
	'Intronic'               => '#02599C',
	'NMD transcript'         => 'orangered',
	'Within non coding gene' => 'limegreen',
	'Upstream'               => '#A2B5CD',
	'Downstream'             => '#A2B5CD',
	'Regulatory region'      => '#4DFEB8',
	'Intergenic'             => '#636363',
);

my $SO_BASE_LINK = 'http://www.sequenceontology.org/miso/current_release/term';


my %cons_rows;
my %consequences;
my %consequences_rank;

for my $cons_set (@OVERLAP_CONSEQUENCES) {

		my $display_term = $cons_set->{display_term};
    my $so_term      = $cons_set->{SO_term};
    my $so_acc       = $cons_set->{SO_accession};
    my $ens_label    = $cons_set->{label};
    my $so_desc      = $cons_set->{description};
    my $ncbi_term    = $cons_set->{NCBI_term} || '-';
    my $rank         = $cons_set->{rank};

		$display_term = $ens_label if (!defined($display_term));
		$display_term = ucfirst(lc($display_term));
	  $display_term =~ s/_/ /g;
	  $display_term =~ s/Nmd /NMD /g;
	  $display_term =~ s/ utr/ UTR/g;
	  $display_term =~ s/rna/RNA/g;
		
    $so_acc = qq{<a rel="external" href="$SO_BASE_LINK/$so_acc">$so_acc</a>};

    my $row = "$so_term|$so_desc|$so_acc|$ncbi_term";

    $cons_rows{$row} = $rank;
		
		push(@{$consequences{$display_term}},$row);
		
		if ($consequences_rank{$display_term}) {
			$consequences_rank{$display_term} = $rank if ($consequences_rank{$display_term} > $rank);
		}
		else {
			$consequences_rank{$display_term} = $rank;
		}	
}


my $cons_table = 
    qq{<table id="consequence_type_table" class="ss">\n<tr>\n\t<th style="width:5px">*</th>\n\t<th>}.
    (join qq{</th>\n\t<th>}, 'Ensembl term', 'SO term', 'SO Description', 'SO accession', 'NCBI term').
    qq{</th>\n</tr>\n};

my $bg;

for my $d_term (sort {$consequences_rank{$a} <=> $consequences_rank{$b}} keys(%consequences)) {

	my $cons_list = $consequences{$d_term};
	my $count = scalar @$cons_list;
	my $rspan = ($count > 1) ? qq{ rowspan="$count"} : '';
	
	my $c = $colour{$d_term};
	my $first_col = (defined($c)) ? qq{\t<td $rspan style="padding:0px;margin:0px"><div style="padding:0px;margin:0px;height:20px;background-color:$c"></div></td>} : qq{<td$rspan></td>};
	
	$cons_table .= qq{<tr$bg>\n$first_col\n\t<td$rspan>$d_term</td>\n};
	
	my $line = 1;
	
	for my $row (sort {$cons_rows{$a} <=> $cons_rows{$b}} @$cons_list) {
		$row =~ s/\|/<\/td>\n\t<td>/g;
		
		$cons_table .= qq{</tr>\n<tr$bg>\n} if ($line !=1 );
		$cons_table .= qq{\t<td>$row</td>\n};
		$line ++;
	}
	$cons_table .= qq{</tr>\n};
	
	$bg = ($bg eq '') ? qq{ class="bg2"} : '';
	
}

$cons_table .= qq{</table>\n};
$cons_table .= qq{<p>* Corresponding colours in the Genome browser.<p>\n};

print $cons_table;

