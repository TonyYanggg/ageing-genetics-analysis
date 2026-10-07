#!/usr/bin/perl
#
# ARCT -- Example 1
#
# This script downloads all human DNA repair genes from NCBI, then BLASTs the CDS in search of homologs in man, mouse, fly, worm, and ecoli and produces the corresponding phylogenetic profile (profile.txt) plus a more detailed report (comparison.txt)

use ARCT::Parser;
use ARCT::Profiles;

# The sequences we want to get (126 DNA repair genes, previously parsed from http://www.cgal.icnet.uk/DNA_Repair_Genes.html
my $query = "NM_003362 NM_014311 NM_003925 NM_003211 NM_002542 NM_012222 NM_002528 NM_002434 NM_024608 NM_001641 NM_014481 NM_013975 NM_006297 NM_001618 NM_005485 NM_002412 NM_000251 NM_002439 NM_000179 NM_002440 NM_002441 NM_000534 NM_000249 NM_000535 NM_014381 D38437 D38500 NM_004628 NM_002874 NM_004344 NM_005053 NM_000380 NM_002945 NM_002946 NM_002947 NM_000122 NM_000400 NM_005316 NM_001515 NM_001516 NM_001517 NM_001799 NM_001239 NM_002431 NM_000123 NM_001983 NM_005236 NM_000234 NM_000082 NM_000124 NM_020196 NM_001923 NM_000107 NM_022362 NM_002875 NM_002877 NM_002876 NM_002878 NM_007068 NM_005431 NM_005432 NM_002879 NM_003579 NM_012415 NM_007295 NM_000059 NM_005732 NM_005590 NM_002485 NM_001469 NM_021141 NM_006904 NM_002312 NM_003401 NM_022487 NM_002452 NM_001948 AB036063 NM_002690 NM_002693 NM_002691 NM_006231 NM_002592 NM_002912 NM_006341 NM_016316 NM_006502 NM_007195 NM_006596 NM_016218 NM_013274 NM_013284 NM_004111 NM_007205 NM_003686 NM_012444 NM_003336 NM_003337 NM_020165 NM_003350 NM_003348 NM_002105 NM_000057 NM_000553 NM_004260 NM_000051 NM_000135 NM_000136 NM_033084 NM_021922 NM_022725 NM_004629 D42045 NM_022836 NM_013347 NM_006020 NM_007254 NM_001184 NM_002853 NM_004584 NM_004507 NM_002873 NM_005657 NM_001274 NM_007194";
# Note: NM_007248 was removed temporarily from RefSeq and so we only used 125 genes

# To retrieve the sequences we can do so one by one
my $start = 1;
my $out_file = "DNArepair.sq";
my $format = "EMBL";		# EMBL works better since it keeps more information that can be useful later but other formats are possible
my $db = 1;

get_single($start, $format, $out_file, $db, $query);

# Then we BLAST everything from the NCBI server according to the following parameters
my $in_file = "DNArepair.sq";
my $out_file = "Blast/test";	# Remember that the Directory Blast must exist in the current directory
my $start = 1;
my $end = 125;
my @query_species = qw /man mouse fly worm ecoli/;

&blast_away ($in_file, $format, $out_file, $start, $end, @query_species);
