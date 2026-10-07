#!/usr/bin/perl
#
# Module to generate HTML reports with hyperlinks to databases. It also features subroutines to search Medline and Google.

use LWP::Simple;
use LWP::UserAgent;
use HTTP::Request;
use HTTP::Response;
use HTML::TreeBuilder;
use HTML::FormatText;

# This subroutine creates outputs as an HTML file with hyperlinks to GenBank, GeneCards, and Medline.
sub info {
	my $out_file = shift @_;
	my $med_search = shift @_;
	open OUT, ">$out_file";

	print OUT "<HTML>\n<TITLE>Auto-generated file with GenBank, GeneCards, and Medline automated searches</TITLE>\n";
	print OUT '<BODY><div align="center"><h1>GenBank, GeneCards, and Medline links</h1></div>';
	print OUT '<TABLE cellSpacing=1 cellPadding=5 width="90%" border=1>';
	print OUT '<TR><TD><B>Gene Symbol</B></TD><TD><B>GenBank</B></TD><TD><B>GeneCards</B></TD><TD><B>Medline</B></TD></TR>';

	foreach $gene (@_) {
		undef @tmp;
		my @tmp = split / /, $gene;
		if ($tmp[1] =~ /\w+/) {
			$acc  = $tmp[1];
		} else {
			$acc = $gene;
		}
		$gene = $tmp[0];
		print OUT "<TR><TD>$gene</TD>";
		print OUT '<TD>'.'<A HREF="http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=nucleotide&tool=arct&cmd=search&term='.$acc.'"/A>GenBank</A></TD>';
		print OUT '<TD>'.'<A HREF="http://bioinfo.weizmann.ac.il/cards-bin/cardsearch.pl?search='.$gene.'"/A>GeneCards</A></TD>';
		print OUT '<TD>'.'<A HREF="http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=search&tool=arct&db=pubmed&term='."$gene+$med_search".'"/A>PubMed</A></TD>';
		print OUT "</TR>";
	}
	print OUT "</TABLE>\n</body>\n</html>";
}

# Subroutine to check the number of entries in Medline for a given query term
sub medline_count {
	my $query = $_[0];
	$query =~ s/\s+/\+/g;
	my $URL ='http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=search&tool=arct&db=pubmed&dispmax=1&doptcmdl=uilist&term='.$query;
	my $retrieve = get($URL);
	my $tree = HTML::TreeBuilder->new_from_content($retrieve);
	my $formatter = HTML::FormatText->new(leftmargin => 0, rightmargin => 50);
	my $text = $formatter->format($tree);
	my @lines = split/\n/, $text;
	my $count = 0;
	foreach (@lines) {
		if (/Item 1 of (\d+).*/) {
			$count = $1;
		}
	}
	return $count;
}

# For a given $query, google_count returns the number of hits in Google plus the correction suggested by Google, if any.
# Example:
# print google_count("aging"),"\n",google_count("aghing");
sub google_count {
	my $query = $_[0];
	$query =~ s/\s+/\+/g;
	my $URL ='http://www.google.com/custom?q='.$query;
	my $ua = new LWP::UserAgent;
	$ua->agent('Mozilla/4.0');
	my $req = new HTTP::Request GET => "$URL";
	my $res = $ua->request($req);
	my $retrieve = $res->content;
	my $tree = HTML::TreeBuilder->new_from_content($retrieve);
	my $formatter = HTML::FormatText->new(leftmargin => 0, rightmargin => 50);
	my $text = $formatter->format($tree);
	my @lines = split/\n/, $text;
	my $count = 0;
	my $correction = "";
	foreach (@lines) {
		if (/Results 1 - 10 of about (\d+,*\d*)/) {
			$count = $1;
		} elsif (/Results 1 - (\d+).*/) {
			$count = $1;
		} elsif (/Did you mean: (.*)/) {
			$correction = $1;
		}
	}
	$count =~ s/,//g;
	return $count, $correction;
}

# Subroutine to check if a given SWISS-PROT variant is a polymorphism or a disease. Make sure you use only numbers and
# not the "VAR_" characters in $query.
# Example:
# print var_type("018942");
#
sub var_type {
	my ($var) = @_;
	my $type = "Unknown";	# Default

	# Make sure $var has leading zeros
	while ( length($var) < 6 ) {
		$var = "0".$var;
	}

	my $URL ='http://au.expasy.org/cgi-bin/get-sprot-variant.pl?VAR_'.$var;

	my $retrieve = get($URL);

	if ($retrieve =~ /polymorphism or unclassified\)<\/td>\s+<td>(\w+)<\/td>/) {
		$type = $1;
	}
	return $type;
}

1;
