use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

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