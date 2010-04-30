#!perl -w

use strict;

use Getopt::Long;
use List::Util qw(min max);
use LRG::LRG;
use LRG::LRGImport;

my $LRG_EXTERNAL_XML = q{ftp://ftp.ebi.ac.uk/pub/databases/lrgex/};
my $LRG_EXTERNAL_ADDR = q{http://www.lrg-sequence.org/};
my $ASSEMBLY = 'GRCh37';

my $xmlfile;
my $outfile;
my $lrgid;

GetOptions(
    'xml=s'         => \$xmlfile,
    'out=s'         => \$outfile,
    'lrg=s'         => \$lrgid,
    'assembly=s'    => \$ASSEMBLY
);

die("An LRG XML input file or LRG identifier is required") unless (defined($xmlfile) || defined($lrgid));

# If no xml file was specified, attempt to fetch it from the website
if (!defined($xmlfile) && defined($lrgid)) {
    
    # Look for XML records in published and pending directories
    my $urls = [
        $LRG_EXTERNAL_XML,
        $LRG_EXTERNAL_XML . 'pending/'
    ];
    
    print STDOUT localtime() . "\tNo input XML file specified for $lrgid, attempting to get it from the LRG server\n";
    
    my $result = LRGImport::fetch_remote_lrg($lrgid,$urls);
    if ($result->{'success'}) {
        $xmlfile = $result->{'xmlfile'};
        print STDOUT localtime() . "\tSuccessfully downloaded XML file for $lrgid and stored it in $xmlfile\n";
    }
    else {
        my $message = "Could not fetch XML file for $lrgid from server. Server says:";
        foreach my $url (@{$urls}) {
            $message .= "\n\t" . $result->{$url}{'message'};
        }
        die($message . "\n");
    }
}

# If no output file has been specified, use the xmlfile and add a .gff suffix
$outfile = $xmlfile . '.gff' if (!defined($outfile));

#ÊGet a LRG object from the xmlfile
my $root = LRG::LRG::newFromFile($xmlfile) or die("Could not create LRG object from $xmlfile");

# In order to get the correct coordinates relative to the assembly, we need to get the mapping
my $mapping = $root->findNodeArray('updatable_annotation/annotation_set/mapping',{'assembly' => $ASSEMBLY}) or die("Could not get mapping between $lrgid and ASSEMBLY");
# Get the first mapping (should only be one anyway)
$mapping = $mapping->[0];

# For now, we will ignore any mapping spans etc. and just assume that the mapping is linear.
# FIXME: Do a finer mapping, this could be done by (temporarily) inserting the LRG into the core db via the LRGMapping.pm module
my $chr_name = $mapping->data()->{'chr_name'};
my $chr_start = $mapping->data()->{'chr_start'};
my $chr_end = $mapping->data()->{'chr_end'};
my $strand;
my $lrg_start = 1e11;
my $lrg_end = -1;
# Get the extreme points of the LRG mapping
my $spans = $mapping->findNodeArray('mapping_span');
foreach my $span (@{$spans}) {
    $strand = $span->data()->{'strand'};
    $lrg_start = min($lrg_start,$span->data()->{'lrg_start'});
    $lrg_end = max($lrg_end,$span->data()->{'lrg_end'});
}
my $strand_ch = ($strand > 0 ? '+' : '-');

my @output;
# Add a GFF entry for the entire LRG
my @row = (
    $chr_name,
    $LRG_EXTERNAL_ADDR,
    'LRG',
    $chr_start,
    $chr_end,
    '.',
    $strand_ch,
    '.',
    $lrgid
);
push(@output,join("\t",@row));

#ÊGet the coordinates for the gene which will be the extreme points of the transcripts
my $track_line;
my $gene_start = 1e11;
my $gene_end = -1;
my $transcripts = $root->findNodeArray('fixed_annotation/transcript') or die("Could not get transcripts from $lrgid");
foreach my $transcript (@{$transcripts}) {
    $gene_start = min($gene_start,$transcript->data()->{'start'});
    $gene_end = max($gene_end,$transcript->data()->{'end'});
    
    my $tr_name = $transcript->data()->{'name'};
    @row[2] = 'exon';
    @row[8] = $lrgid . '_' . $tr_name;
    
    # Loop over the exons and print them to the GFF
    my $nodes = $transcript->{'nodes'};
    
    # Add a track line for this transcript
    $track_line = "track name=\"$lrgid\_$tr_name\" description=\"Transcript $tr_name for gene of $lrgid\" color=128,128,255";
    push(@output,$track_line);
    for (my $i=0; $i<scalar(@{$nodes}); $i++) {
        next if ($nodes->[$i]->name() ne 'exon');
        
        my $phase;
        # Get the start phase of the exon from the last phase of the intron
        if ($i > 0 && $nodes->[$i-1]->name() eq 'intron') {
            $phase = $nodes->[$i-1]->data()->{'phase'};
        }
        else {
            $phase = '.';
        }
        
        # Get the exon information
        if ($nodes->[$i]->name() eq 'exon') {
            my $exon_start = $nodes->[$i]->findNode('lrg_coords')->data()->{'start'};
            my $exon_end = $nodes->[$i]->findNode('lrg_coords')->data()->{'end'};
            @row[3] = ($strand > 0 ? ($chr_start + $exon_start - 1) : ($chr_end - $exon_end + 1));
            @row[4] = ($strand > 0 ? ($chr_start + $exon_end - 1) : ($chr_end - $exon_start + 1));
            @row[7] = $phase;
            push(@output,join("\t",@row));
        }
    }
}
#ÊShift off the first element (the LRG region)
my $lrgrow = shift(@output);
@row[2] = 'gene';
@row[3] = ($strand > 0 ? ($chr_start + $gene_start - 1) : ($chr_end - $gene_end + 1));
@row[4] = ($strand > 0 ? ($chr_start + $gene_end - 1) : ($chr_end - $gene_start + 1));
@row[8] = $lrgid . '_g1';
# Add the gene entry to the top of the array
unshift(@output,join("\t",@row));
# Add a track line for the LRG gene entry
$track_line = "track name=\"$lrgid" . "_g1" . "\" description=\"Genomic region spanned by gene of $lrgid\" color=128,128,255";
unshift(@output,$track_line);

#ÊAdd the LRG region entry to the top
unshift(@output,$lrgrow);
# Add a track line for the LRG region
$track_line = "track name=\"$lrgid\" description=\"Genomic region spanned by $lrgid\" color=128,128,255";
unshift(@output,$track_line);

# Add browser line
my $browser_line = "browser position $chr_name\:$chr_start\-$chr_end";
unshift(@output,$browser_line);

#ÊOpen the outfile for writing
open(GFF,'>',$outfile) or die("Could not open $outfile for writing");
print GFF join("\n",@output);
close(GFF);
