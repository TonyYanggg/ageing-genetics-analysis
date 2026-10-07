#!/usr/bin/perl
#
# TF module for researching transcription factors
#
# These are the basic modules for TFBS detection, including the creation of a database of TF matrices, retrieving information from the database and searching for TFBS in nucleotide sequences

use TFBS::DB::TRANSFAC;
use TFBS::DB::FlatFileDir;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Graphics;
use Bio::SeqFeature::Generic;
use lib "$ENV{HOME}/projects/bioperl-live";

# The subroutine TFscan scans the FASTA file you pass on to it and then creates a graphical display of the results
sub TF_scan {
	my $filename = shift @_;
	my $db_name = shift @_;
	my $out_dir = shift @_;
	my $query;
	mkdir $out_dir, 0755 unless (-e $out_dir);

	# Getting the sequence(s) from the file
	my $seqio  = Bio::SeqIO->new( '-format' => 'fasta' , -file => "$filename");
	# Getting the TF for each sequence in the file
	while ($seqobj = $seqio->next_seq()) {
		undef %tfstats;
		undef @tfhits;

		# If the user does not name the sequence then the name becomes the FASTA ID
		if (defined($_[0])) {
			$query = shift @_;
		} else {
			$query = $seqobj->display_id();
			$query = substr $query, 1, 8 unless length($query) < 9;
		}
		my $seq_length = $seqobj->length();

		print "Processing $query\n";

		# The full length sequence
		my $full_length = Bio::SeqFeature::Generic->new(-start=>-"$seq_length",-end=>-1,  -seq_id=>"$query upstream region");

		# Create a panel with the sequences' length
		$panel = Bio::Graphics::Panel->new(-segment=> $full_length,
								-key_style => 'between',
								-width  => 800,
								-pad_left => 10,
								-pad_right => 10
								);

		# We create a track in the panel to add the labels with the start and end of the scale, which is a new track
		my $track = $panel->add_track(-glyph  => 'arrow',
				-bgcolor => 'white',
				-fgcolor => 'white',
				-label  => 1,
				-height  => 0,
				);
		my $feature = Bio::SeqFeature::Generic->new(-seq_id=>"-$seq_length",
							-start=>"-$seq_length",-end=>"-$seq_length");
		$track->add_feature($feature);
		my $feature = Bio::SeqFeature::Generic->new(-seq_id=>"-1",
							-start=>-1,-end=>-1);
		$track->add_feature($feature);

		# Then we add a ruler
		$panel->add_track($full_length,
			-glyph   => 'arrow',
			-fgcolor => 'black',
			-tick => 1,
			-east => 'true',
			-label => 0,
			);

		# Then we add the sequence we're searching
		$panel->add_track($full_length,
				-glyph  => 'generic',
				-bgcolor => 'blue',
				-label  => 1,
				);

		# Then for the TFBS, we can use a different track for each TF, each with a different color
		my @colors = qw(aliceblue antiquewhite aqua aquamarine azure beige bisque black blanchedalmond blue bluevoilet brown burlywood cadetblue chartreuse chocolate coral cornflowerblue cornsilk crimson cyan darkblue darkcyan darkgoldenrod darkgray darkgreen darkkhaki darkmagenta darkolivegreen darkorange darkorchid darkred darksalmon darkseagreen darkslateblue darkslategray darkturquoise darkviolet deeppink deepskyblue dimgray dodgerblue firebrick floralwhite forestgreen fuchsia gainsboro ghostwhite gold goldenrod gray green greenyellow honeydew hotpink indianred indigo ivory khaki lavender lavenderblush lawngreen lemonchiffon lightblue lightcoral lightcyan lightgoldenrodyellow lightgreen lightgrey lightpink lightsalmon lightseagreen lightskyblue lightslategray lightsteelblue lightyellow lime limegreen linen magenta maroon mediumaquamarine mediumblue mediumorchid mediumpurple mediumseagreen mediumslateblue mediumspringgreen mediumturquoise mediumvioletred midnightblue mintcream mistyrose moccasin navajowhite navy oldlace olive olivedrab orange orangered orchid palegoldenrod palegreen paleturquoise palevioletred papayawhip peachpuff peru pink plum powderblue purple red rosybrown royalblue saddlebrown salmon sandybrown seagreen seashell sienna silver skyblue slateblue slategray snow springgreen steelblue tan teal thistle tomato turquoise violet wheat white whitesmoke yellow yellowgreen);
		srand(time|$$);

		my $tmp = 0;

		# You declare which matrices to find in the following line. Yet be sure to check how many matrices exist for each TF
		my @matrices = qw |AP-1ii AP-2 C/EBPi CRE-BP1/c-Jun CREBiii E2Fi E47i Elk-1i HNF-4 HSF1 Max MyoD NF-AT V$NFKB_Q6 SRFi STAT1 STAT3 Sp1i TATA c-Myb c-Myc/Maxi p300 p53i|;

		# Opening the necessary files
		my $db = TFBS::DB::FlatFileDir->connect("$db_name");

		# Now for the loop in which we search each of the matrices
		foreach $matrix (@matrices) {

			# First we get the matrix
			# Getting a matrix by name or ID
			if ($matrix =~ /V\$.*/) {
				$pwm = $db->get_Matrix_by_ID("$matrix",'PWM');
				$matrix = $pwm->name;
			} else {
				$pwm = $db->get_Matrix_by_name("$matrix",'PWM');
			}
			$matrix =~ s/i//g;        # To remove the trailing "i" referring to different versions of the TF

			# The threshold is set in accordance with the length of the binding site
			my $limit = 105 - ($pwm->length)*1.7;
			$limit = $limit."%";

			# The actual search
			my $siteset = $pwm->search_seq(-seqobj=>$seqobj, -threshold=>"$limit");
			# The result is a TFBS::Site object

			# Then we add a track for each TF if we have a TFBS
			unless ($siteset->size() == 0) {
				$track = $panel->add_track(
						-glyph  => 'generic',
						-label  => 1,
						-bgcolor  =>  $colors[int(rand(139))],
						-font2color => 'black',
						-description => $matrix
						);
			}
			++$tmp;

			# To get each TFBS, do a loop
			my $site_iterator = $siteset->Iterator(-sort_by => "start");
			my $start = 0;
			my $end = 0;
			my $score = 0;
			while (my $site = $site_iterator->next())  {
				# If two TFBS overlap then draw just one of them (I do not take into account the score because that does not affect the graphical display, but it should be taken into acocunt when interpreting results);
				next unless ($end < $site->start());
				$start = $site->start();
				$end = $site->end();
				$score = $site->score();

				# In this case I add features to the track but it's also possible to declare the track based on already set features without even the need to create a track object
				my $feature = Bio::SeqFeature::Generic->new(-seq_id=>($start-$seq_length),
									-start=>($start-$seq_length),-end=>($end-$seq_length));
				$track->add_feature($feature);
				# Array for later building the general view of all TFBS
				push @tfhits, $start;
				# Hash for later building the clusters
				if (defined($tfstats{$start})) {
					my @x = split / /, $tfstats{$start};
					$end = $x[0] unless $end > $x[0];
					$score += $x[1];
					$tfstats{$start} = "$end $score";
				} else {
					$tfstats{$start} = "$end $score";
				}
			}
		}

		# Drawing of a general view of all TFBS. Useful mostly for large sequences.
		unless ($seq_length < 200) {
			$full_length->seq_id("$query TFBS distribution");
			$panel->add_track($full_length,
			-glyph  => 'generic',
			-bgcolor => 'blue',
			-label  => 1,
			);
			my $track = $panel->add_track(-glyph  => 'triangle',
					-bgcolor  => 'black',
					-bump     => +1,
					-height   => 10,
					-orient => N,
					-label    => 0,
					);

			while (@tfhits) {
				my $start = shift @tfhits;
				my $feature = Bio::SeqFeature::Generic->new( -seq_id=>"$score",
								-start=>($start-$seq_length),-end=>(($start+int($seq_length/100))-$seq_length));
				$track->add_feature($feature);
			}

			# Drawing a multi-colored gradient of TFBS intensity.
			# Then we add the labels with the scale, which is a new track
			my $track = $panel->add_track(-glyph  => 'arrow',
					-bgcolor => 'white',
					-fgcolor => 'white',
					-label  => 1,
					-height  => 0,
					);
			my $feature = Bio::SeqFeature::Generic->new(-seq_id=>"TFBS clusters",
								-start=>-$seq_length,-end=>-$seq_length);
			$track->add_feature($feature);

			my $limit = int($seq_length/100);
			# What follows is a very simple cluster identification loop, though it can be modified.
			foreach (1 .. $limit) {
				clusters($_*10, $seq_length);
			}
		}

		open OUT, ">$out_dir"."/"."$query";
		print OUT $panel->png;
	}
}

# Subroutine to process the clusters found at different levels of threshold and draw them on the $panel
sub clusters {
	my ($threshold, $seq_length) = @_;
	undef %cluster;
	undef @cluster;
	my $past = 0;
	my $past_score = 0;
	foreach $key (sort by_score keys %tfstats) {
		my @tmp =  split / /, $tfstats{$key};
		$end = shift @tmp;
		$score = shift @tmp;
		my $y = $score;
		if (($past+$threshold) > $key && $past > 0) {
			if (defined($cluster[$#cluster]) && (($cluster[$#cluster-1]+$threshold) > $key)) {
				$score += $cluster[$#cluster];
				$cluster[$#cluster-1] = $end;
				$cluster[$#cluster] = $score;
			} else {
				push @cluster, $past, $end, $score;
			}
		}
		$past = $key;
		$past_score = $y;
	}

	# Getting the highest score for reference
	undef @numbers;
	my $maximus = 0;

	foreach (1 .. (($#cluster+1)/3)) {
		my $cv = int($cluster[$_*3-1]/($cluster[$_*3-2]-$cluster[$_*3-3])*$threshold);
		unshift @numbers, $cv;
	}
	$maximus = $numbers[0];

	my $track = $panel->add_track(-glyph => 'graded_segments',
		-bgcolor => 'blue',
		-label  => 0,
		-font2color => 'black',
		-bump => 0,
		-height  => 8,
		-min_score => 0,
		-max_score => int($maximus+$threshold-10)
		);

while (@cluster) {
	my $start = shift @cluster;
	my $end = shift @cluster;
	my $temp = shift @cluster;
	my $score = int($temp/($end-$start)*$threshold);
	my $feature = Bio::SeqFeature::Generic->new( -score => $score, -seq_id=>"$score",
					-start=>($start-$seq_length),-end=>($end-$seq_length));
	$track->add_feature($feature);
}
}

# To get a list of matrices in the archive use the listall subroutine
sub TF_list {
	my $db_name = $_[0];

	# Connect to the database
	my $db = TFBS::DB::FlatFileDir->connect("$db_name");

	# Get the whole set of matrices, then create an iterator sorted by name (ID also useful)
	my $matrixset = $db->get_MatrixSet(-matrixtype=>"PFM");
	my $mx_iterator = $matrixset->Iterator(-sort_by=>'name');

	# to print heading if normal output
	printf("\n %-10s%15s%25s\n",
		'MatrixID', 'Name','Length');

	# Then it just loops through the matrices to print each of them
	while (my $pfm = $mx_iterator->next())  { # For each matrix in the set
		printf(" %-10s%15s%25s\n",
		$pfm->ID, $pfm->name, $pfm->length);
	}
}

# Subroutine that returns an array with all the names of TF stored in a local database
sub TF_names {
	# Connect to the database
	my $db = TFBS::DB::FlatFileDir->connect("$_[0]");

	# Get the whole set of matrices, then create an iterator sorted by name (ID also useful)
	my $matrixset = $db->get_MatrixSet(-matrixtype=>"PFM");
	my $mx_iterator = $matrixset->Iterator(-sort_by=>'name');

	# Then it just loops through the matrices to print each of them
	while (my $pfm = $mx_iterator->next())  { # For each matrix in the set
		my $name = $pfm->name;
		push @names, $name  unless ($name =~ /.+i+/);
	}
	return @names;
}

# To know how many matrices for a given TF are present
sub tf_max {
	my($db_name, $query) = @_;

	# Connect to the database
	my $db = TFBS::DB::FlatFileDir->connect("$db_name");

	# Get the whole set of matrices, then create an iterator sorted by name (ID also useful)
	my $matrixset = $db->get_MatrixSet(-matrixtype=>"PFM");
	my $mx_iterator = $matrixset->Iterator(-sort_by=>'name');

	# to print heading if normal output
	printf("\n %-10s%15s%25s\n",
		'MatrixID', 'Name','Length');

	# Then it just loops through the matrices to print each of them
	while (my $pfm = $mx_iterator->next())  { # For each matrix in the set
		my $name = $pfm->name;
		if ($name=~/$query.*/) {
			printf(" %-10s%15s%25s\n",
			$pfm->ID, $name, $pfm->length);
		}
	}
}

# To create a database of TFs
sub TF_db {
	my $db_name = shift @_;
	my(@names, $id, $hdb);
	my $count = 0;
	my @query_id = undef;

	print $_[0];
	if ($_[0] =~ /\s*(V\$.+)/ ) {
		@query_id = @_;
	} else {
		my $input_file = $_[0];
		open INPUT, "<"."$input_file";
	}

	# We start by connecting to TRANSFAC and creating or opening the FlatFileDir
	my $db = TFBS::DB::TRANSFAC->connect(-accept_conditions => 1);
	if (-e $db_name) {
		$hdb = TFBS::DB::FlatFileDir->connect("$db_name");
	} else {
		$hdb = TFBS::DB::FlatFileDir->create("$db_name");
	}

	# If we have a file, we open it and process all IDs we can find. Otherwise we assume we already have IDs
	if (defined($input_file)) {
		foreach $id (<INPUT>) {
			if ( $id =~ /\s*(V\$.+)/ ) {
				$id =~ s/\s//g;
				process_matrix ($id, $db, $hdb);
			}
		}
	} else {
		foreach (@query_id) {
			process_matrix ($_, $db, $hdb);
		}
	}
}

# The subroutine that fetches the matrix and then stores it in the local database. To be used from within TF_db
sub process_matrix {
	my($id, $db, $hdb) = @_;
	$pwm = $db->get_Matrix_by_ID("$id",'PFM');

	return unless $pwm; # Just a precaution in case something goes wrong and I can't get the matrix
	my $name = $pwm->name();

	# Checking if the TF is already in the archive. If it is, then we add a "i" to it and archive it with a new name
	$x = 0;
	foreach $tmp (@names) {
		if ($tmp eq $name) {
			$name .= "i";
			$pwm->name($name);
		}
	}
	print "Stored matrix: ".$name."\n";     # To make sure we are doing something worthwhile

	# If not, let's archive it
	if ( $x == 0) {
		push @names, $hdb;
		$hdb->store_Matrix($pwm);
	}
}

# Subroutines for sorting hashes and lists
sub by_score { $tfstats{$a} <=> $tfstats{$b} }
sub by_number { $a <=> $b }

1;
