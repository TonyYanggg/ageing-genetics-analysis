 #!/usr/bin/perl
#
# ARCT -- Example 2
#
# This script downloads the promoter of the human Werner's syndrome gene (hWRN) from the first (hWRN has many) transcription to the translation start site. It then searches for TFBS and creates a PNG file with the results. Have in mind that you need to have a local TF database (created with TF_db, see main.pl).

use ARCT::Parser;
use ARCT::TF;

# First we retrieve the promoter into a local file
my $filename = "wrn.fa";
my @query = qw /NM_000553/;

get_promoter($filename, @query);

# Then we find the TFBS and save the results
my $db_name = "/gene/Perl/hTFBS";	# Path where you should have your FlatFileDir with the TF
my $out_dir = "TFtest";
my @names = qw/hWRN/;

# To invoque TF_scan is simple
TF_scan($filename, $db_name, $out_dir, @names);