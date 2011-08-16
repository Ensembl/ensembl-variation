=head1 LICENSE

 Copyright (c) 1999-2011 The European Bioinformatics Institute and
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

=head1 NAME

Bio::EnsEMBL::Variation::Utils::ProteinFunctionUtils

=head1 DESCRIPTION

This module defines several subroutines used to encode and decode
the compressed format used to store the SIFT and PolyPhen predictions
in the variation databases. Normally users will not use this module
directly, but should instead use the various methods provided to access
these predictions on a Bio::EnsEMBL::Variation::TranscriptVariationAllele
object.

=cut

package Bio::EnsEMBL::Variation::Utils::ProteinFunctionUtils;

use strict;
use warnings;

use POSIX qw(ceil);
use List::Util qw(max);
use Compress::Zlib;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(Exporter);

my $DEBUG = 0;

our @EXPORT_OK = qw(
    @ALL_AAS
    $AA_LOOKUP
    $NO_PREDICTION
    $HEADER
    prediction_from_string 
    prediction_to_short
    prediction_from_short
    compress_prediction_string
    expand_prediction_string
    compute_offset
);

# user-defined constants

# a header which lets us identify prediction strings and
# to check if they may have been corrupted

our $HEADER     = 'VEP';

# the format we use when 'pack'ing our predictions

my $PACK_FORMAT = 'v';

# a hash mapping qualitative predictions to numerical values

my $PREDICTION_TO_VAL = {
    polyphen => {
        'probably damaging' => 0,
        'possibly damaging' => 1,
        'benign'            => 2,
        'unknown'           => 3,
    },

    sift => {
        'tolerated'     => 0,
        'deleterious'   => 1,
    },
};

# all valid amino acids

our @ALL_AAS = qw(A C D E F G H I K L M N P Q R S T V W Y);

# we use a short with all bits set to represent the lack of a prediction in 
# an (uncompressed) prediction string, we will never observe this value
# as a real prediction even if we set all the (6) prediction bits because we 
# limit the max score to 1000 so the 10 score bits will never all be set

our $NO_PREDICTION = pack($PACK_FORMAT, 0xFFFF);

# constants derived from the the user-defined constants

# the maximum number of distinct qualitative predictions used by any tool

my $MAX_NUM_PREDS = max( map { scalar keys %$_ } values %$PREDICTION_TO_VAL ); 

# the number of bits used to encode the qualitative prediction

my $NUM_PRED_BITS = ceil( log($MAX_NUM_PREDS) / log(2) );

throw("Cannot represent more than ".(2**6-1)." predictions") if $NUM_PRED_BITS > 6;

# a hash mapping back from a numerical value to a qualitative prediction

my $VAL_TO_PREDICTION = {
    map {
        my $tool = $_; 
        $tool => {
            map {
                $PREDICTION_TO_VAL->{$tool}->{$_} => $_ 
            } keys %{ $PREDICTION_TO_VAL->{$tool} }
        }
    } keys %$PREDICTION_TO_VAL
};

# a hash from amino acid single letter code to a numerical value

our $AA_LOOKUP = { map {$ALL_AAS[$_] => $_} 0 .. $#ALL_AAS };

# the number of valid amino acids

our $NUM_AAS = scalar(@ALL_AAS);


=head2 prediction_to_short

  Arg[1]      : string $tool - one of 'sift' or 'polyphen' 
  Arg[2]      : string $pred - one of the possible qualitative predictions of the tool 
  Arg[3]      : string $prob - the prediction score (with 3 d.p.s of precision)
  Description : converts a prediction and corresponding score into a 2-byte short value
  Returntype  : the packed short value
  Exceptions  : throws if either the tool or prediction arguments aren't recognised 
  Status      : At Risk
  
=cut

sub prediction_to_short {
    my ($tool, $pred, $prob) = @_;

    throw("Unrecognised prediction tool '$tool'") 
        unless defined $PREDICTION_TO_VAL->{$tool};

    # we only store 3 d.p. so multiply by 1000 to turn our
    # probability into a number between 0 and 1000. 
    # 2^10 == 1024 so we need 10 bits of our short to store 
    # this value
    
    my $val = $prob * 1000;

    # we store the prediction in the top $NUM_PRED_BITS bits
    # so look up the numerical value for the prediction, 
    # shift this $NUM_PRED_BITS bits to the left and then set
    # the appropriate bits of our short value
    
    $pred = lc($pred);
    
    my $pred_val = $PREDICTION_TO_VAL->{$tool}->{$pred};

    throw("No value defined for prediction: '$pred'?")
        unless defined $pred_val;

    $val |= ($pred_val << (16 - $NUM_PRED_BITS));

    printf("p2s: $pred ($prob) => 0x%04x\n", $val) if $DEBUG;
    
    $val = pack $PACK_FORMAT, $val;

    return $val;
}

=head2 prediction_from_short

  Arg[1]      : string $tool - one of 'sift' or 'polyphen' 
  Arg[2]      : string $pred - the packed short value 
  Description : converts a 2-byte short value back into a prediction and a score
  Returntype  : a list containing 2 values, the prediction and the score
  Exceptions  : throws if the tool argument isn't recognised 
  Status      : At Risk
  
=cut

sub prediction_from_short {
    my ($tool, $val) = @_;

    # check it isn't our special null prediction

    if ($val eq $NO_PREDICTION) {
        print "no pred\n" if $DEBUG;
        return undef;
    }

    # unpack the value as a short

    $val = unpack $PACK_FORMAT, $val;

    # shift the prediction bits down and look up the prediction string

    my $pred = $VAL_TO_PREDICTION->{$tool}->{$val >> (16 - $NUM_PRED_BITS)};

    # mask off the top 6 bits reserved for the prediction and convert back to a 3 d.p. float

    my $prob = ($val & (2**10 - 1)) / 1000;

    printf("pfs: 0x%04x => $pred ($prob)\n", $val) if $DEBUG;

    return ($pred, $prob);
}

=head2 compress_prediction_string

  Arg[1]      : string $string - the prediction string to compress 
  Description : compresses a prediction string with gzip
  Returntype  : the compressed string
  Exceptions  : throws if the string is an unexpected length, or if the compression fails
  Status      : At Risk
  
=cut

sub compress_prediction_string {
    my ($string) = @_;

    return undef unless $string;

    throw("Prediction string is an unexpected length")
         unless (length($string) % $NUM_AAS) == 0;

    # prepend a header, so we can tell if our string has been mangled, or
    # is compressed etc.

    my $gzipped = Compress::Zlib::memGzip($HEADER.$string)
        or throw("Failed to gzip: $gzerrno");

    return $gzipped;
}

=head2 header_ok

  Arg[1]      : string $string - the prediction string to check 
  Description : checks if a prediction string has the expected header
  Returntype  : boolean
  Exceptions  : none
  Status      : At Risk
  
=cut

sub header_ok {
    my ($string) = @_;
    return undef unless $string;
    return substr($string,0,length($HEADER)) eq $HEADER;
}

=head2 expand_prediction_string

  Arg[1]      : string $gzipped - the compressed prediction string 
  Description : uncompresses a compressed prediction string
  Returntype  : the uncompressed string
  Exceptions  : throws if the header is incorrect, or if the decompression fails
  Status      : At Risk
  
=cut

sub expand_prediction_string {
    my ($gzipped) = @_;

    return undef unless $gzipped;

    my $raw = Compress::Zlib::memGunzip($gzipped)
        or throw("Failed to gunzip: $gzerrno");

    throw("Malformed prediction string") unless header_ok($raw);
    
    return $raw;
}

=head2 compute_offset

  Arg[1]      : int $pos - the desired position in the peptide 
  Arg[2]      : string $aa - the desired mutant amino acid
  Description : computes the correct offset into a prediction string for a given
                peptide position and mutant amino acid
  Returntype  : the integer offset
  Exceptions  : none
  Status      : At Risk
  
=cut

sub compute_offset {
    my ($pos, $aa) = @_;

    my $offset = length($HEADER) + ( ( (($pos-1) * $NUM_AAS) + $AA_LOOKUP->{$aa} ) * 2);

    return $offset;
}

=head2 prediction_from_string

  Arg[1]      : string $tool - one of 'sift' or 'polyphen' 
  Arg[2]      : string $string - the (possibly compressed) prediction string
  Arg[3]      : int $pos - the desired position in the peptide 
  Arg[4]      : string $aa - the desired mutant amino acid
  Description : returns the prediction and score from the specified tool for 
                the requested position and mutant amino acid in the given
                prediction string
  Returntype  : a list containing 2 values, the prediction and the score
  Exceptions  : throws if the tool argument isn't recognised, if either the 
                position or amino acid are invalid, or if the prediction string
                looks to be malformed
  Status      : At Risk
  
=cut

sub prediction_from_string {
    my ($tool, $string, $pos, $aa) = @_;

    unless (header_ok($string)) {
        # maybe the string is still compressed so we try to uncompress it, 
        # this will throw if its still not valid
        $string = expand_prediction_string($string);
    }

    $aa = uc($aa) if defined $aa;

    throw("Invalid position: $pos") unless (defined $pos && $pos > 0);
    
    throw("Invalid amino acid: $aa") unless (defined $aa && defined $AA_LOOKUP->{$aa});

    # compute our offset into the prediction string

    my $offset = compute_offset($pos, $aa);
    
    print "offset: $offset\n" if $DEBUG;
    
    if ($offset + 1 > length($string)) {
        warning("Offset outside of prediction string for position $pos?");
        return undef;
    }
    
    my $pred = substr($string, $offset, 2);

    return prediction_from_short($tool, $pred);
}

1;


