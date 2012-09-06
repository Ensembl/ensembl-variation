#!/localsw/bin/perl -w

use strict;

# Read the order file
my $idorder = read_order($ARGV[0]);

# Read the 
open FP,"<$ARGV[1]";

my $prevtid = '';

my %tid_hash;
my $last = 0;
while (<FP>) {
  my ($junk,$junk2,$tid,$rest) = split;

  if ($prevtid ne $tid) {
    if ($last) {
      $tid_hash{$prevtid}->{end} = $last-1;
  #    print "$prevtid " . $tid_hash{$prevtid}->{start} . " " . $tid_hash{$prevtid}->{end} . "\n";
    }
    $tid_hash{$tid}->{start} = $last;
    $prevtid = $tid;
  }
  $last = tell(FP);
}
$tid_hash{$prevtid}->{end} = $last;
#print "$prevtid " . $tid_hash{$prevtid}->{start} . " " . $tid_hash{$prevtid}->{end} . "\n";

my $data;
foreach my $id (@$idorder) {
  if (!exists($tid_hash{$id})) {
    print STDERR "Note: Didn't find $id in $ARGV[1] (may be OK if no TVs for this transcript)\n";
  } else {
    seek(FP,$tid_hash{$id}->{start},0);
    read(FP, $data, $tid_hash{$id}->{end} - $tid_hash{$id}->{start} + 1);
    print $data;
  }
}

sub read_order {
  my ($idfile) = shift;

  open IDFP,"<$idfile";

  my @idorder;
  while (<IDFP>) {
    my ($junk,$tid,$rest) = split;
    push @idorder,$tid;
  }

  close IDFP;

  return \@idorder;
}  

