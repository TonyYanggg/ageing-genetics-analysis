#!/usr/bin/perl
#
# Module to generate phylogenetic profiles and perform statistical analysis. This module will probably not be updated again
# since I'm working on more efficient implementations of the same functions.

use Bio::SeqIO;
use Bio::Seq;
use Bio::SeqFeatureI;
use Bio::AnnotationI;
use Bio::Tools::Run::RemoteBlast;

sub blast_here {
	my($in_file, $start, $end, @query_species) = @_;
	my %species = init_prof($in_file, $start, $end, @query_species);

	my $count = $start-1;
	my %total;
	for $x ($start .. $end) {
		print "\nProgress: $count of $end\n\n";
		my $blast_report = new Bio::SearchIO ('-format' => 'blast',
							'-file'   => $in_file.$x.".txt");

		my $result = $blast_report->next_result;
		process_result($result,  \%species, \%total, @query_species);
		$count++;
	}
	final_stats($count - $start + 1, \%total, @query_species);
}

sub blast_away {
	my($in_file, $in_format, $out_file, $start, $end, @query_species) = @_;
	my %species = init_prof($in_file, $start, $end, @query_species);

	# We load the sequences and create our sequence objects
	my $seqio  = Bio::SeqIO->new( '-format' => "$in_format" , -file => "$in_file");
	if ($start > 1) {
		for (2 .. $start) {
			$seq = $seqio->next_seq();
			$x++;
		}
	}

	while (my $seq = $seqio->next_seq()) {
		last if $x == $end;
		$ref_name{$seq->display_id()} = $seq->accession_number();	# Get the reference sequence's accession number, if available
		push @query, $seq;
		$x++;
	}

	$end = ($#query+1) unless defined($end);

	# Creates the object for the method
	my $factory = Bio::Tools::Run::RemoteBlast->new(
					'-prog'  => 'blastp',	# We will BLAST proteins, so BLASTP is advizable
					'-data' => 'nr',
					_READMETHOD => "Blast",
					'-expect' => '1e-10'	# This is the E cut-off value; if you are using evolutionary distant species, you should use a higher value
				);

	# Multiple submission to the NCBI server
	my $blast_report = $factory->submit_blast(\@query);

	my $max_number = 100;
	my $trial = 0;
	my %total;
	my $count = $start-1;
	# Retrieve the results from the queue at NCBI
	while ( my @rids = $factory->each_rid ) {

		print STDERR "\nSorry, maximum number of retries $max_number exceeded\n" if $trial >= $max_number;
		last if $trial >= $max_number;
		$trial++;

		print STDERR "waiting... ".(5*$trial)." units of time\n" ;
		# RID = Remote Blast ID (e.g: 1017772174-16400-6638)
		foreach my $rid ( @rids ) {
			my $rc = $factory->retrieve_blast($rid);
			if( !ref($rc) ) {
				if( $rc < 0 ) {
						# retrieve_blast returns -1 on error
						$factory->remove_rid($rid);
				}
				# retrieve_blast returns 0 on 'job not finished'
				sleep 5*$trial;
			} else {
				#---- Blast done ----
				# Save the resulting BLAST report
				print "\nProgress: $count of $end\n\n";
				$factory->save_output(">".$out_file.($count+1).".txt") if defined($out_file);

				$factory->remove_rid($rid);
				my $result = $rc->next_result;

				process_result($result,  \%species, \%total, @query_species);
				$count++;
			}
		}
	}
	final_stats($count - $start + 1, \%total, @query_species);
}

sub init_prof {
	my($in_file, $start, $end, @query_species) = @_;

	# In order to add another species, all you must do is add the corresponding pair as below
	my %species = ("man", "Homo sapiens",
			"mouse", "Mus musculus",
			"rat", "Rattus norvegicus",
			"fugu", "Takifugu rubripes",
			"worm", "Caenorhabditis elegans",
			"fly", "Drosophila melanogaster",
			"yeast", "Saccharomyces cerevisiae",
			"ecoli", "Escherichia coli");

	# If we don't start the query at the first sequence, it's assumed there was a connection problem and so we add the results we get to the already existing file instead of creating a new file
	my $tmp = ">";
	if ($start > 1) { $tmp = ">>"; }

	# The simple 1/0 phylogenetic profile will be saved as profile.txt
	open RESULT, "$tmp"."profile.txt";
	# A more detailed comparison between proteins is saved as comparison.txt
	open REPORT, "$tmp"."comparison.txt";

	# We create the headers in both output files
	unless ($start > 1) {
		my $xy = "Reference gene       ";
		foreach (@query_species) {
			my $int = length($_);
			foreach (1 .. (6 - $int)) {
				$xy .= " ";
			}
			$xy .= $_;
			#foreach (1 .. (6 - $int)) {
				$xy .= "   ";
			#}
		}
		print RESULT "$xy\n";
		print REPORT "Multiple species comparison\n";
	}
	return %species;
}

sub process_result {
	my $result = shift @_;
	my $species = shift @_;
	my $total = shift @_;
	my @query_species = @_;

	# We start by print the reference gene name and length
	my $query_name = $result->query_name();
	my $query_length = $result->query_length();
	if ($ref_name{$query_name} =~ /\w+/) {
		$query_acc = $ref_name{$query_name};
	} else {
		$query_acc = "Not available";
	}
	my $xy = "$query_name ($query_length aa):";
	foreach (1 .. (24 - length($xy))) {
		$xy .= " ";
	}
	print RESULT $xy;
	print REPORT "\nReference gene: $query_name, $query_length aa & $query_acc\n";

	undef %retrieve;
	# Here we use an object from the Bio::Search::Hit::GenericHit class
	while( $hit = $result->next_hit ) {
		foreach my $organism (@query_species) {
			if ($hit->description() =~ /$species->{$organism}/i ) {
				my @out = split / /, $retrieve{$organism};
				my $score = $hit->raw_score();
				if ( $score > $out[1] ) {
					undef @out;
					my $name = $hit->name();
					my $length = $hit->length();
					my $ident = $hit->frac_identical();
					my $conserv = $hit->frac_conserved();
					my $hitac = $hit->accession();
					print "$query_name hit found in $organism\nScore: ", $score, " frac_identical: ", $ident, "\n";	# So the user knows what is going on and keeps track of the progress
					@out = ($name, $score, $length, $ident, $conserv, $hitac);
					$retrieve{$organism} = "@out";
				}
			}
		}
	}

	# Now we treat the results statistically and print them into file
	foreach $org (@query_species) {
		undef @find;
		@find = split / /, $retrieve{$org};

		my $found = 0;
		# This is the formula that determines whether a hit is a functionally similar protein or not. At present, the hit must have at least 60% of the length of the query sequence and have a conservation of at least 40%
		if ( $find[1] != 0 && ($find[2]/$query_length) > 0.6 && $find[4] > 0.4) {
			$found = 1;
			my @z = undef;
			@z = split / /, $total->{$org};
			++$z[3];
			$total->{$org} = "@z";
		}
		print RESULT "$found        ";

		# The following piece of code prints a more detailed report
		print REPORT "$org: ";
		if  ( $find[1] != 0) {
			my @z = undef;
			@z = split / /, $total->{$org};
			++$z[0];		# Number of hits
			$z[1] += $find[3];		# Similarity
			$z[2] += $find[4];		# Positives
			$total->{$org} = "@z";
			print REPORT "$find[5], $find[2] aa, ";
			printf REPORT "Similarity = %.1f%%, Positives = %.1f%%\n", $find[3]*100, $find[4]*100;
		} else {
			print REPORT "Not found\n";
		}

	}

	print RESULT "\n";
}

sub final_stats {
	my $count = shift @_;
	my $total = shift @_;
	my @query_species = @_;

	# Printing the final stats
	print RESULT "Final results         ";
	foreach $org (@query_species) {
		my @z = undef;
		@z = split / /, $total->{$org};
		$z[3] = 0 unless defined($z[3]);
		my $xy = $z[3]."/".$count;
		foreach (1 .. (9 - length($xy))) {
			$xy .= " ";
		}
		print RESULT $xy;
	}

	print RESULT "\nPercentage             ";
	print REPORT "\nFinal results:\n";
	foreach $org (@query_species) {
		my @z = undef;
		@z = split / /, $total->{$org};
		my $temp = "";
		my $xy = ($z[3]/$count)*100;
		printf RESULT "%.0f%%", $xy;
		foreach (1 .. (8 - length(int("$xy")))) {
			$temp .= " ";
		}
		print RESULT $temp;
		# For the detailed results
		print REPORT "$org: ";

		my $xy = ($z[0]/$count)*100;
		if ($z[0] == 0) {
			$aa_sim = 0;
			$aa_con = 0;
		} else {
			$aa_sim = $z[1]/$z[0]*100;
			$aa_con = $z[2]/$z[0]*100;
		}
		printf REPORT "Homologs found: %.1f%%, Overall aa similarity = %.1f%%, Overall aa positives = %.1f%%\n", $xy, $aa_sim, $aa_con;
	}
	print RESULT "\n";
	print REPORT "\n";
}

# Subroutine to find 1:1 orthologs, assuming man to be the reference organism
sub find_orthologs {
	my ($organism, @user) = @_;
	my %query;	# hash with the reference => query accession number pairs
	my %names;	# hash with the names of the genes we are loooking for

	foreach (@user) {
		my @x = split /\&/, $_;
		my $query_acc = "";
		my $find_acc = "";
		foreach (@x) {
			my @z = split /,/, $_;
			if ($z[0] eq "man") {	# At present, we get the query accession number from the man hit because I normally keep the old accession number (often from a DNA sequence) in the EMBL sequence file. This can be changed in the future quite easily.
				$query_acc = $z[1];
				$query_acc =~ s/ //;
				my @y = split /,/, $x[0];
				$names{"$query_acc"} = $y[0];
			}
			if ($z[0] eq $organism) {
				$find_acc = $z[1];
				$find_acc =~ s/ //;

			}
		}
		$query{"$query_acc"} = $find_acc;	# %query holds the pair of accession numbers we want to compare
	}

	foreach (keys %query) {
		my $query_id = $query{$_};
		if ( $query_id =~ /NP(\w)+/) {
			$database = new Bio::DB::RefSeq;
		} else {
			$database = new Bio::DB::GenPept;
		}
		$seq = $database->get_Seq_by_acc("$query_id");
		push @query_seq, $seq;
	}

	my $factory = Bio::Tools::Run::RemoteBlast->new(
					'-prog'  => 'blastp',
					'-data' => 'nr',
					_READMETHOD => "Blast",
					'-expect' => '1e-20'	# Since we are looking for human proteins, the score can be quite high
				);

	my $blast_report = $factory->submit_blast(\@query_seq);

	my $max_number = 100;
	my $trial = 0;
	my $count = 0;
	# Retrieve the results from the queue at NCBI
	while ( my @rids = $factory->each_rid ) {

		print STDERR "\nSorry, maximum number of retries $max_number exceeded\n" if $trial >= $max_number;
		last if $trial >= $max_number;
		$trial++;

		print STDERR "waiting... ".(5*$trial)." units of time\n" ;
		# RID = Remote Blast ID (e.g: 1017772174-16400-6638)
		foreach my $rid ( @rids ) {
			my $rc = $factory->retrieve_blast($rid);
			if( !ref($rc) ) {
				if( $rc < 0 ) {
						# retrieve_blast returns -1 on error
						$factory->remove_rid($rid);
				}
				# retrieve_blast returns 0 on 'job not finished'
				sleep 5*$trial;
			} else {
				#---- Blast done ----
				$factory->save_output(">orthol".($count+1).".txt");	# Save blast reports

				$count++;
				undef @hits;
				print "Progress: $count of $end\n";
				$factory->remove_rid($rid);
				my $result = $rc->next_result;
				while( $hit = $result->next_hit ) {
					if ($hit->description() =~ /Homo sapiens/i ) {
						my $score = $hit->raw_score();
						my $hitac = $hit->accession();
						@hits = ($hitac, $score) unless ($hits[1] > $score);
					}
				}
				print "The best human homologue is: $hits[0]\n";
				foreach (keys %query) {
					if ($_ eq $hits[0]) {
						print "Found match of $names{$_}\n";
						push @orthologs, $_;
					}
				}

			}
		}
	}

	open OUT, ">orthologs.txt";
	print OUT "List of 1:1 orthologs in $organism:\n\n";
	foreach $find (keys %query) {
		print OUT "$names{$find}: ";
		my $tmp = 0;
		foreach (@orthologs) {
			$tmp = -1 unless $_ ne $find;
		}
		if ($tmp = -1) {
			print OUT "found 1:1 ortholog\n";
		} else {
			print OUT "did not found 1:1 ortholog\n";
		}
	}
}

1;
