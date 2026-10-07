#!/usr/bin/env perl  -w
#
# Module to use Gibbs program and then get a new TFBS matrix.

use TFBS::DB::FlatFileDir;
use TFBS::PatternGen::SimplePFM;

# subroutine to use the Gibbs algorithm
sub gibbs {
	my $path = shift @_;
	my $in_file = shift @_;
	my $lengths = shift @_;
	my @temp = split /\+/, $lengths;
	$lengths = shift @temp;
	foreach (@temp) {
		$lengths .= ",".$_;
	}
	my $nuc = "";
	my $cutoff = "";
	my $expected = "";
	shift @_ if $_[0] eq "aa";

	if (scalar(@_) >= 3) {
		$nuc = " -d -p2 -".shift @_ if ($_[0]) && ($_[0] eq "n");
		$expected = " ".shift @_ if $_[0] > 1;
		$cutoff = " -C".shift @_ if ($_[0]) && ($_[0] <= 1);
	} elsif (scalar(@_) > 0) {
		$nuc = " -d -p2 -".shift @_ if ($_[0]) && ($_[0] eq "n");
		$expected = " ".shift @_ if $_[0] > 1;
		$cutoff = " -C".shift @_ if ($_[0]) && ($_[0] <= 1);
	}

	#system "purge $in_file 1000 -q";	# to remove repeated sequences; advizable when comparing many sequences
	system "$path $in_file $lengths$expected$nuc$frag -f$cutoff";
}

# subroutine to extract and store the matrix found by gibbs. if no output is specified, the matrix will be sent to STDOUT
sub extract_matrix {
	my $in_file = shift @_;
	my $out = "";
	my $db;
	undef my @sequences;
	if ($_[0]) { $out = $_[0]; }

	open IN, "<$in_file";
	while (<IN>) {
		if (/[ACGTacgt]/) {
			chomp;
			s/\W//g;
			push @sequences, $_;
		}
	}

	my $patterngen = TFBS::PatternGen::SimplePFM->new(-seq_list=>\@sequences);
	my $pfm = $patterngen->pattern(); # $pfm is now a TFBS::Matrix::PFM object

	if ($out =~ /\w+/) {

		if (-e  $out) {
			$db = TFBS::DB::FlatFileDir->connect("$out");
		} else {
			$db = TFBS::DB::FlatFileDir->create("$out");
		}
		$db->store_Matrix($pfm);

	} else {
		my $prettystring = $pfm->prettyprint();
		print $prettystring;
	}
}

1;

