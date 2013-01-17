use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/legal/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut


while(<>) {
	my $copy = $_;
	
	# get sequence-y bits
	while($_ =~ m/([ACGTN-]+)/g) {
		my $s = $1;
		reverse_comp(\$s);
		
		$copy =~ s/$1/$s/;
	}
	
	print $copy;
}