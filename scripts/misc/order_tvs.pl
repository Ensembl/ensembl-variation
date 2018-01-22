#!/localsw/bin/env/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.



=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut


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

