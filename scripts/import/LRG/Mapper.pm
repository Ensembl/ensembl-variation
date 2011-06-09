#! perl -w

###
# A framework for mapping a sequence to something using different methods
# 

use strict;
use warnings;

package Mapper;

use LRG::SMALT;
use LRG::Samtools;
use LRG::Revcomp;
use LRG::Exonerate;
use LRG::bsub;
use Bio::EnsEMBL::Utils::Scalar qw(wrap_array);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use List::Util qw (min max);

#ÊUse autoload for get/set methods
our $AUTOLOAD;

sub DESTROY { }

# Some default values
sub defaults {
    return {
        'tmpdir' => '/tmp/',
        'endsize' => 500,
        '_suffix' => time()
    };
}

sub permitted {
    return [
        'query',
        'target',
        'hash',
        'tmpdir',
        'pid',
        'endsize',
        'logprefix',
        'mappings',
        '_suffix'
    ];
}

sub AUTOLOAD {
    my $self = shift;
    my $type = ref($self) or die("$self is not an object");

    #ÊThe name of the called subroutine
    my $name = $AUTOLOAD;
    # Strip away the pre-pended package info
    $name =~ s/.*://;

    # Check that the subroutine should exist for this module
    unless (grep {/^$name$/} @{$self->_get_set('_permitted')}) {
        die "Can't access `$name' field in class $type";
    }

    # Call the get/set method
    return $self->_get_set('_' . $name,@_);
}

#ÊConstructor
sub new {
    
    my $class = shift;
    
    my $self = bless({},$class);
    $self->initialize(@_);
    
    return $self;
    
}

# Initialize some variables
sub initialize {
    my $self = shift;
    my %vals = @_;
    
    #ÊThe get/set methods that exist in this module
    $self->_get_set('_permitted',$self->permitted());

    # Set the pid of this process
    $self->pid($$);

    #ÊSet the default values
    my $DEFAULTS = $self->defaults();
    map { $self->$_($DEFAULTS->{$_}) } keys(%{$DEFAULTS});
    
    #ÊSet any field passed in via the parameter hash
    map { $self->$_($vals{$_}) } keys(%vals);
     
}

# Internal get/set method
sub _get_set {
    my $self = shift;
    my $attribute = shift;
    my $value = shift;
    
    if (defined($value)) {
        $self->{$attribute} = $value;
    }
    
    return $self->{$attribute};
}

# The procedure for mapping a large sequence to the genome looks roughly like:
#   1. Do a paired read mapping to the genome using SMALT with the 5' and 3' ends of the query
#   2. Extract the genomic sequence where the pair has been "anchored" using faidx of samtools
#   3. Do a detailed alignment using exonerate
#ÊThe Mapper module wraps this all up in a pipeline. The input needed to the module is the query and the target database.
# Additional parameters can be set
sub do_mapping {
    my $self = shift;
    
    # Check that a query has been specified
    die ("A query sequence or file has to be specified") unless(defined($self->query()));
    # Check that a target has been specified
    die ("A target database has to be specified") unless(defined($self->target()));
    
    my $anchors = 5;
    my $anchor_length = 200;
    my $min_anchor_spacing = 2000;
    my $max_anchor_spacing;
    my $max_length_tolerance;
    my %alignments;
    
    # Create wrappers for the pipeline components
    my $smalt = SMALT->new('hash' => $self->hash());
    # Set additional parameters for SMALT
    $smalt->extra_parameters(['-f','cigar']);
    
    my $samtools = Samtools->new('program' => 'faidx');
    my $exonerate = Exonerate->new('model' => 'affine:bestfit', 'extra_parameters' => ['--bestn','1','--showalignment','yes','--showvulgar','no','--ryo',q{"=QMATCH=%qi,%qab,%qae,%qS\n=TMATCH=%ti,%tab,%tae,%tS\n=VULGAROUT=%V\n=SCORE=%s\n"},'--refine','region']);
    
    # Process each query seq
    my $seqs = $self->from_fasta($self->query()); 
    
    while (my ($id,$seq) = each(%{$seqs})) {
        
        my ($short_id) = $id =~ m/^(\S+)/;
        printf("Processing \%s...\n",$id);
        
        #ÊExtract a few samples from the sequence and "anchor" them to the genome
        
        #ÊCalculate the spacing between the anchors 
        my $anchor_spacing = int((length($seq) - $anchor_length*$anchors)/($anchors - 1));
        
        #ÊIf the anchor spacing is too small, say 2kb, just use the entire sequence
        if ($anchor_spacing < $min_anchor_spacing) {
            $anchors = 1;
            $anchor_length = length($seq);
        }
          
        #ÊExtract the anchor sequences          
        my @anchorseq;
        my $offset = 0;
        for (my $i=0; $i<$anchors; $i++) {
            push(@anchorseq,substr($seq,$offset,$anchor_length));
            $offset += $anchor_length + $anchor_spacing;            
        }
        
        #ÊWrite the anchors to a fasta file
        my $anchorfile = $self->to_fasta(\@anchorseq);
        
        # Use the anchorfile as input for SMALT
        $smalt->inputfile($anchorfile);
        
        #ÊSet an output file for SMALT
        my $outfile = $self->tmpfile() . ".smalt";
        $smalt->outputfile($outfile);
        
        # Wrap the SMALT process up in a bsub and execute
        printf("\tAnchoring sequence \%s to reference genome...\n",$id);
        my $logprefix = $self->logprefix() || $self->tmpdir() . "/mapper_" . $self->pid() . "_" . $short_id;
        my $bsub = bsub->new('logout' => qq{$logprefix\.out}, 'logerr' => qq{$logprefix\.err}, 'jobname' => 'smalt_' . $self->pid() . "_$short_id", 'job' => $smalt);
        $bsub->execute();
        
        #ÊParse the SMALT output and determine the most promising genomic region where the anchors map in a sensible way
        printf("\tFinding most promising genomic region where anchors map...\n");
        my $best_span = parse_anchors($outfile,$anchors,$anchor_length,$anchor_spacing) or die ("Could not identify the genomic region where the query sequence, $id, should be anchored");
        
        #ÊExtract the genomic region, including some wiggle room on the flanks
        printf("\tDetermined the best location to be \%s:\%d-\%d:\%d, extracting sequence, including \%d bp flanking...\n",$best_span->[0],$best_span->[1],$best_span->[2],$best_span->[3],$anchor_spacing);
        my $genomic_start = $best_span->[1] - $anchor_spacing;
        my $genomic_end = $best_span->[2] + $anchor_spacing;
        $samtools->extra_parameters([$self->target(),$best_span->[0] . ":" . $genomic_start . "-" . $genomic_end]);
        my $genomic_seq = $samtools->execute();
        # Dump the sequence to a file
        my $genomicfile = $self->to_fasta($genomic_seq);
         
        #ÊDo the full alignment using Exonerate
        # Write the query seq to a tempfile
        my $queryfile = $self->to_fasta($seq);
        $exonerate->query($queryfile);
        $exonerate->target($genomicfile);
        
        printf("\tAligning query sequence to genomic region...\n");
        #ÊBsub the exonerate job
        $bsub->jobname('exonerate_' . $self->pid() . "_$short_id");
        $bsub->job($exonerate);
        $bsub->execute();
        
        unlink($anchorfile);
        unlink($outfile);
        unlink($genomicfile);
        unlink($queryfile);
        
        printf("\tParsing alignment...\n");
        my $alignment = $self->parse_alignment($bsub->logout(),1,$genomic_start);
        $alignments{$id} = $alignment;
    }
    
    $self->mappings(\%alignments);
}

sub parse_anchors {
    my $outfile = shift;
    my $anchors = shift;
    my $anchor_length = shift;
    my $anchor_spacing = shift;
    
    my $max_anchor_spacing = 3*$anchor_spacing;
    my $max_length_tolerance = $anchor_spacing;
    
    # Determine the location of the anchors
    my @locations;
    open (FH,"<",$outfile) or die ("Could not open SMALT output file: $outfile for parsing");
    while (<FH>) {
        chomp;
        my ($q,$qs,$qe,$qo,$sn,$ss,$se,$so) = $_ =~ m/^\S+\s+(\S+)\s+(\d+)\s+(\d+)\s+([\-\+])\s+(\S+)\s+(\d+)\s+(\d+)\s+([\-\+])/;
            
        # Replace the orientations with integers
        $qo =~ s/^\+$//;
        $qo .= '1';
        $so =~ s/^\+$//;
        $so .= '1';
        push(@locations,[$q,$qs,$qe,$qo,$sn,$ss,$se,$so]);
    }
        
    # Sort the locations according to chromosome and position
    @locations = sort {$a->[4] cmp $b->[4] || $a->[5] <=> $b->[5]} @locations;

    #ÊParse the anchor locations and try to find a genomic span where they occur together at the expected intervals. 
    # Allow the distance between to mapped anchors to be at most $max_anchor_spacing and require them to be on the same chromosome
    my @span;
    foreach my $location (@locations) {
            
        # If we are within an accepted span, update the end coordinate of the span and add orientation and the contributing anchor count 
        if (scalar(@span) && $span[-1]->[0] eq $location->[4] && ($location->[5] - $span[-1]->[2]) < $max_anchor_spacing) {
            $span[-1]->[2] = $location->[6];
            $span[-1]->[3] += ($location->[3] * $location->[7]);
            $span[-1]->[4]++;
        }
        #ÊElse, start a new span
        else {
            push(@span,[$location->[4],$location->[5],$location->[6],($location->[3] * $location->[7]),1]);
        }
    } 
        
    # Determine the best span by picking the one having the length closest to the length of the query
    my $query_length = $anchors*$anchor_length + ($anchors-1)*$anchor_spacing;
    @span = sort {abs($b->[2] - $b->[1] + 1 - $query_length) <=> abs($a->[2] - $a->[1] + 1 - $query_length) || $a->[4] <=> $b->[4]} @span;
        
    my $best_span;
    foreach my $s (@span) {   
        #ÊPick the first one that fulfills the criteria
        if (abs($s->[2] - $s->[1] + 1 - $query_length) < $max_length_tolerance && $s->[4] >= ($anchors - 1)) {
            $best_span = $s;
            $best_span->[3] /= abs($best_span->[3]);
            last;
        }
    }
    
    return $best_span;
}

# Parse the exonerate output
sub parse_alignment {
    my $self = shift;
    my $alignfile = shift;
    my $qoffset = shift || 1;
    my $toffset = shift || 1;
    my $allow_query_reversed = shift;
    
    my ($tid,$tstart,$tend,$tstrand);
    my ($qid,$qstart,$qend,$qstrand);
    my ($qseq,$tseq,$target);
    my ($vulgar,$score);
    
    #ÊOpen the output file
    open(FH,"<",$alignfile) or die ("Could not open alignment file $alignfile for parsing");
    while (<FH>) {
        chomp;
        
        # Read the aligned sequences
        if (m/^\s*\d+\s+\:\s+(\S+)\s+\:\s+\d+\s*$/) {
            my $seq = $1;
            $seq =~ s/\-//g;
            $qseq .= $seq if (!$target);
            $tseq .= $seq if ($target);
            $target = !$target;
        }
        #ÊParse the query or target coordinates
        elsif (m/^\=([Q|T])MATCH\=(\S+),(\d+),(\d+),([\+\-])$/) {
            ($qid,$qstart,$qend,$qstrand) = ($2,$3,$4,($5 eq '+' ? 1 : -1)) if ($1 eq 'Q');
            ($tid,$tstart,$tend,$tstrand) = ($2,$3,$4,($5 eq '+' ? 1 : -1)) if ($1 eq 'T');
        }
        #ÊParse the vulgar output
        elsif (m/^\=VULGAROUT\=(.+)/) {
            $vulgar = $1;
        }        
        #ÊGet the score
        elsif (m/^\=SCORE\=(\d+)/) {
            $score = $1;
        }        
    }
    close(FH);
    
    #ÊParse the vulgar string and store the details of the alignment
    my @pairs;
    my $qaoffset = 0;
    my $taoffset = 0;
    while ($vulgar =~ m/([A-Z]+)\s+(\d+)\s+(\d+)/ig) {
        my ($type,$qlen,$tlen) = ($1,$2,$3);
        
        my $qsegment = substr($qseq,$qaoffset,$qlen);
        my $tsegment = substr($tseq,$taoffset,$tlen);
        
        # If this is a indel
        if ($type eq 'G') {
        
            my $segment = [
                $type,
                $qaoffset,
                $qaoffset + ($qlen - 1),
                $taoffset,
                $taoffset + ($tlen - 1),
                $qstrand * $tstrand
            ];    
            if ($qlen > 0) {
                #$segment->{'type'} = 'QINS';
                push(@{$segment},$qsegment);
                push(@{$segment},'-' x length($qsegment)); 
            }
            else {
                #$segment->{'type'} = 'QDEL';
                push(@{$segment},'-' x length($tsegment));
                push(@{$segment},$tsegment); 
            }
            
            push(@pairs,$segment);
        }
        elsif ($type eq 'M') {
            
            push(@pairs,[
                'DNA',
                $qaoffset,
                $qaoffset + ($qlen - 1),
                $taoffset,
                $taoffset + ($tlen - 1),
                $qstrand * $tstrand
            ]);
            for (my $i=0; $i<min($qlen,$tlen); $i++) {
                if (substr($qsegment,$i,1) ne substr($tsegment,$i,1)) {                    
                    #ÊCreate a segment for the mismatch
                    push(@pairs,[
                        $type,
                        $qaoffset + $i,
                        $qaoffset + $i,
                        $taoffset + $i,
                        $taoffset + $i,
                        $qstrand * $tstrand,
                        substr($qsegment,$i,1),
                        substr($tsegment,$i,1)                            
                    ]);
                }
            }
        }
        else {
            warn ("Unrecognized alignment type in vulgar string: $type");
        }
        
        $qaoffset += $qlen;
        $taoffset += $tlen;
        
    } 
    
    # Update the coordinates to accomodate the offset into the sequence
    foreach my $pair (@pairs) {
        foreach my $ref (['query',$qstart,$qoffset,$qstrand],['target',$tstart,$toffset,$tstrand]) {
            my ($obj,$start,$offset,$strand) = @{$ref};
            my $posoffset = ($obj eq 'query' ? 0 : 2);
            if ($strand > 0) {
                $pair->[$posoffset + 1] += $start + $offset;
                $pair->[$posoffset + 2] += $start + $offset;
            }
            else {
                my $s = $pair->[$posoffset + 1];
                $pair->[$posoffset + 1] = ($start - 1) - $pair->[$posoffset + 2] + $offset;
                $pair->[$posoffset + 2] = ($start - 1) - $s + $offset;
                
            }
        }
        #ÊIf the query is opn the negative strand but we're supposed to use the forward, reverse complement
        unless ($allow_query_reversed || $qstrand > 0) {
            reverse_comp(\$pair->[6]);
            reverse_comp(\$pair->[7]);
        }
    }

    # Store and return the overall mapping in the format that is expected
    # Keep only the chromosome name from the target id
    ($tid) = $tid =~ m/^(\S+)\:\d+\-\d+/; 
     return {
         'lrg_id' => $qid,
         'chr_name' => $tid,
         'lrg_start' => min($qstart,$qend) + $qoffset,
         'lrg_end' => max($qstart,$qend) + $qoffset - 1,
         'chr_start' => min($tstart,$tend) + $toffset,
         'chr_end' => max($tstart,$tend) + $toffset - 1,
         'strand' => ($qstrand * $tstrand),
         'score' => $score,
         'pairs' => \@pairs 
    };
        
}

# Generate a unique temporary filename
sub tmpfile {
    my $self = shift;
    my $tmpfile = $self->tmpdir() . "/mapper_tmpfile_" . $self->pid() . "_" . $self->_suffix();
    $self->_suffix(int($self->_suffix())+1);
    return $tmpfile;
}

#ÊWill dump a sequence to a (possibly temporary) fasta file. 
#ÊCan also be passed a filename in which case it will do nothing
sub to_fasta {
    my $self = shift;
    my $seq = shift;
    my $outfile = shift;
    
    #ÊFirst of all, if the $seq is in fact an existing file, just return the name of that file
    my $seqs = wrap_array($seq);
    return $seq unless ((grep {is_dna($_)} @{$seqs}) || !(-e $seq && -f $seq ));
    
    #ÊIf $outfile is undefined, generate a temporary filename
    $outfile = $self->tmpfile() . ".fa" unless (defined($outfile));
    
    # Open a filehandle to the outfile
    open(FA,">",$outfile) or die ("Could not open output fasta file '$outfile' for writing");
    
    # Write the sequence(s) in fasta format. For simplicity, use the name of the file and an iterator as identifier
    my $n = 1;
    map {
        my $name = "";
        if (substr($_,0,1) ne '>') {
            ($name) = $outfile =~ m/^(?:.*\/)?(.+)$/;
            $name = ">" . $name . "_$n\n";
        }
        printf(FA "%s\%s\n",$name,$_); 
        $n++;
    } @{$seqs};
    
    # Close the file
    close(FA);
    
    #ÊReturn the name of the outfile
    return $outfile;
}

# Gets the sequences from a fasta file (as a hash) or, if the supplied argument is a sequence, return a hash with the sequence
sub from_fasta {
    my $self = shift;
    my $file = shift;
    
    # If the supplied file is in fact a nucleotide sequence, wrap it in the hash and return it
    return {'seq_' . length($file) => $file} if (is_dna($file));
    
    #ÊElse, parse it as a fasta file
    return parsefasta($file);
} 

sub is_dna {
    my $str = shift;
    my $ratio = shift || 0.75;
    
    return 0 unless (defined($str));
    
    my $nt_count = lc($str) =~ tr/[acgt]//;
    my $r = $nt_count/length($str);
    
    return ($r >= $ratio);
}

sub parsefasta {
    my $infile = shift;
    
    my %seqs;
    my $id;
    open (FH,"<",$infile) or die ("Could not open $infile for parsing");
    while (<FH>) {
        chomp;
        if (substr($_,0,1) eq ">") {
      
            $id = substr($_,1);
            $seqs{$id} = "";
        }
        else {
            $_ =~ s/\s*//;
            $seqs{$id} .= $_;
      
        }
    }
    close(FH);

    return \%seqs;
}

1;
