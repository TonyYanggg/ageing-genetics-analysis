#!/usr/bin/perl
#
# Module for data-mining. This is a restricted version based on our published results. As soon as more of my
# results using this module are published, I will include more subroutines.
#
# Note: All the files from a given type and organisms must be in the name of the organism as will be specified in the query

# algorithm for data-mining
sub mine_data {
	my $dir_file = shift @_;
	my @g_names = split / /, shift @_;
	my($out_file, @groups) = @_;
	my @in_files = dir_or_file($dir_file);
	my $nuc_fix = 0;
	$nuc_fix = 30 unless $dir_file =~ /.*\/CDS\/$/;	# more stringent for nucleotides

 	foreach my $in_file (@in_files) {
		my $count = 0;
		foreach my $group (@groups) {
			++$count;
			my $out_path = ">".$g_names[$count-1]."_".$out_file;
			$out_path = ">".$out_path if -e $g_names[$count-1]."_".$out_file;
			open OUT, $out_path;
       			print OUT "Group correction algorithm\n\n";

			# first we get the alignment of all the sequences in the reference alignment
			my $in = Bio::AlignIO->new(-file => "<$in_file",
						'-format' => 'clustalw');
			my $ref_aln = $in->next_aln();

			my $num_seq = $ref_aln->no_sequences();	# number of sequences in the overall alignment

			# and we get its statistics
			undef my @ref_matches;
			undef my @ref_strong;
			undef my @ref_weak;
			undef my @ref_mismatches;
			my $ref_sim = align_stats($ref_aln, \@ref_matches, \@ref_strong, \@ref_weak, \@ref_mismatches);
			my $ref_consensus = $ref_aln->consensus_string(30+$nuc_fix);

			# then we keep only the sequences under study (@$group)
			re_align(1, \$ref_aln, @$group);

			# and get their stats
			undef my @query_matches;
			undef my @query_strong;
			undef my @query_weak;
			undef my @query_mismatches;
			my $query_sim = align_stats($ref_aln, \@query_matches, \@query_strong, \@query_weak, \@query_mismatches);
			my $query_consensus = $ref_aln->consensus_string(30+$nuc_fix);
			my $control_consensus = $ref_aln->consensus_string(60+$nuc_fix);

			# if the query sequences are very divergent, we use a less stringent threshold
			if ($query_sim < 60) {
				@query_matches = (@query_matches, @query_strong, @query_weak);
			} elsif ($query_sim < 70) {
				@query_matches = (@query_matches, @query_strong);
			} elsif ($query_sim > 90) {
				@query_mismatches = (@query_mismatches, @query_weak);
			}

			# the data-mining is divided into finding disrupted funcitonal motifs and finding novel motifs

			# to find novel functional motifs, we check whether the conserved sequences in
			# the query sequence are different in the other sequences

			# we check the sequences without those in the query
			my $new_in = Bio::AlignIO->new(-file => "<$in_file",
						'-format' => 'clustalw');
			my $comp_aln = $new_in->next_aln();
			re_align(0, \$comp_aln, @$group);	# alignment with all sequences except those in the query

			# probably a bug in Bioperl obliges me to use this system
			save_align("temp.txt", $comp_aln);
			$new_in = Bio::AlignIO->new(-file => "<temp.txt",
					'-format' => 'clustalw');
			$comp_aln = $new_in->next_aln();

			undef my @comp_matches;
			undef my @comp_strong;
			undef my @comp_weak;
			undef my @comp_mismatches;
			my $comp_sim = align_stats($comp_aln, \@comp_matches, \@comp_strong, \@comp_weak, \@comp_mismatches);
			my $comp_seq = $comp_aln->consensus_string(40+$nuc_fix);

			# we exclude the positions that are similar between both query and reference alignments
			undef my @false_positives;
			for (0..length($query_consensus)) {
				push @false_positives, $_ if (substr($query_consensus, $_, 1) eq substr($comp_seq, $_, 1));
				push @false_mismatches, $_ if (substr($query_consensus, $_, 1) ne "?");
				push @false_positives, $_ if (substr($control_consensus, $_, 1) eq "?");
				push @false_positives, $_ if (substr($ref_consensus, $_, 1) eq "?" && $query_sim > 80);
			}
			my @results_novel = unique_positions(\@query_matches, \@false_positives, 0, 0);

			# we also exclude the positions that are repeated in the sequences not in the query
			undef @false_positives;
			for (1..($comp_aln->no_sequences())) {
				my $seq = $comp_aln->get_seq_by_pos($_);
				my $seq_string = $seq->seq();
				for (0..length($ref_consensus)) {
					push @false_positives, $_ if substr($query_consensus, $_, 1) eq substr($seq_string, $_, 1);
				}
			}
			@results_novel = unique_positions(\@results_novel, \@false_positives, $num_seq, 10);

			# and we can exclude the positions that are not conserved in the reference alignment
			# this is advizable if using very similar query sequences
			if ($query_sim > 80) {
				@results_novel = unique_positions(\@results_novel, \@comp_mismatches, 0, 0);
			} elsif ($query_sim > 90) {
				@results_novel = unique_positions(\@results_novel, \@comp_weak, 0, 0);
			}

	      		# findinding disrupted functional motifs takes the mismatches of our query and
			# checks whether these are conserved in the other sequences

			my @results_disrupted = unique_positions(\@query_mismatches, \@comp_mismatches, 0, 0);
			@results_disrupted = unique_positions(\@results_disrupted, \@false_mismatches, 0, 0);

			# to save the resulting alignment to file (useful as a control)
			#save_align("m".$out_file, $ref_aln);

			# print out results
			$in_file =~ /^.*\/(\S+)al.fa$/;
			my $name = $1;

			print OUT "  ### ".$g_names[$count-1]." group ###\n\n";
			print OUT "$name (length = ".length($query_consensus)."): $ref_sim% identity in reference alignment and $query_sim% in query\n\n";
			evolutionary_selection(\@query_matches, \@query_strong, \@query_weak, \@query_mismatches, \@ref_matches, \@ref_strong, \@ref_weak, \@ref_mismatches) unless $nuc_fix == 30 || scalar(@$group) <= 1;
			output_unique($name, "+", $query_consensus, \@ref_matches, \@results_novel);
			output_unique($name, "-", $query_consensus, \@ref_matches, \@results_disrupted);
			print OUT "\n";
			close OUT;
		}
     	}
	unlink "temp.txt";
}

# simple subroutine to save alignments, used mostly as a control
sub save_align {
	my($out_file, $aln) = @_;
	my $out  = Bio::AlignIO->new(-file => ">$out_file",
			'-format' => 'clustalw');
	$out->write_aln($aln);
}

# subroutine to print out the results into OUT, which must be previously opened
sub output_unique {
	my($name, $type, $consensus, $temp_m, $temp_r) = @_;
	my @results = @$temp_r;
	my @matches = @$temp_m;
	if (@results && scalar(@results) != 0) {
		if ($type eq "-") {
			print OUT scalar(@results)." unique mismatches found in $name (".int(length($consensus)+1)." aa): ";
			$type = "|";
		} else {
			print OUT scalar(@results)." unique aa found in $name (".int(length($consensus)+1)." aa): ";
			$type = "^";
		}
		foreach (1..scalar(@results)) {
			$results[$_-1] += 1;
		}
		foreach (1..scalar(@matches)) {
			$matches[$_-1] += 1;
		}
		@results = sort s_by_number @results;
		print OUT "@results"."\n";
		print OUT "Resulting sequence:\n";
		if (length($consensus) > 50) {
			for my $count (0..int(length($consensus)/60)) {
			  	print OUT substr($consensus, $count*60, 60)."\n";
				my $cursor = 0;
				# we print the results
				foreach (@results) {
					unless ($_ <= $count*60 || $_ > $count*60+60) {
						printf OUT "%".($_-$count*60-$cursor)."s", $type;						$cursor = $_-$count*60;
					}
				}
				print OUT "\n" unless $cursor == 0;
				# we print the conserved regions in the reference alignment
				#my $cons_cursor = 0;
				$cursor = 0;
				foreach (@matches) {
					unless ($_ <= $count*60 || $_ > $count*60+60) {
						printf OUT "%".($_-$count*60-$cursor)."s", "-";						$cursor = $_-$count*60;
					}
				}
				print OUT "\n" unless $cursor == 0;
			}
		} else {
			print OUT $consensus."\n";
		}

		print OUT "\n";
	} else {
		if ($type eq "-") {
			print OUT "No unique mismatches found in $name\n";
		} else {
			print OUT "No unique matches found in $name\n";
		}
	}
}

sub s_by_number { $a <=> $b }

# subroutine to check matches between two aligments
sub compare_matches {
	my($ref, $query) = @_;
	undef my @found;
	foreach my $ref_pos (@$ref) {
		foreach my $pos (@$query) {
			push @found, $pos if ($pos == $ref_pos);
		}
	}
	return @found;
}

#subroutine that removes ($remove_or_keep == 0) or keeps ($remove_or_keep == 1) all sequences except those names in @names from an aligment
sub re_align {
	my($remove_or_keep, $aln, @names) = @_;
	my $num_seq = $$aln->no_sequences();
	for (0..($$aln->no_sequences()-1)) {
		my $seq = $$aln->get_seq_by_pos($num_seq-$_);
		my $flag = 0;
		foreach (@names) {
			$flag = 1 if $seq->id() eq $_;
		}
		$$aln->remove_seq($seq) unless $flag == $remove_or_keep;
	}
}

# subroutine to find unique mismatches
# it returns the numbers in $query not found in $ref
sub unique_positions {
	my($query, $ref, $num_seq, $threshold) = @_;
	undef my @results;
	foreach my $f_pos (@$query){
		my $flag = -($num_seq-1)*($threshold/100);
		foreach my $ref_pos (@$ref) {
			++$flag if $f_pos == $ref_pos;
		}
		push @results, $f_pos unless $flag > 0;
	}
	return @results;
}

# to get a sequence's string from an aligment just use the next subroutine
sub get_sequence {
	my($name, $aln) = @_;
	foreach my $seq ( $$aln->each_seq_with_id($name) ) {
		return $seq->seq() if $seq->id() eq $name;
	}
}

# checks whether the query is a file or a directory and returns the file(s) to query in an array
# it also checks which files in a directory are alignments in ARCT general format
sub dir_or_file {
	my($dir_file) = @_;

	if (-d $dir_file) {
		$dir_file .= "/" unless ($dir_file =~ /\/$/);
	}

	undef my @in_files;
	if ($dir_file =~ /\/$/) {
		my @files = <$dir_file*>;
	      	foreach (@files) {
			push @in_files, $_ if (/al\.fa$/ || /\.clw$/);	# finds alignments
		}
	} else {
		push @in_files, $dir_file;
	}
	return @in_files;
}

# subroutine to analyse an aligment for each type of matches and mismatches
# it returns the identity of the alignment
sub align_stats {
	my($aln, $mt, $st, $wk, $ms) = @_;
	my $line = $aln->match_line();
	for (0..(length($line)-1)) {
		my $char = substr($line, $_, 1);
		if ($char eq ":") {
			push @$st, $_;
		} elsif ($char eq ".") {
			push @$wk, $_;
		} elsif ($char eq "*") {
			push @$mt, $_;
		} else {
			push @$ms, $_;
		}
	}
	return int($aln->percentage_identity());
}

# algorithm to search each group of organisms or all organisms grouped into one array
sub compare_sequences {
	my ($ref, $in_path, $out_path, $temp_t, $temp_o) = @_;

	# first we prepare the directory tree
	mkdir $out_path, 0755 unless -e  $out_path;

	my @types = @$temp_t;
	my @organisms = @$temp_o;

	# then we prepare the data
	prepare_data($in_path, $out_path, $ref, \@organisms, \@types);

	# finally, we make the calculations
	process_data($out_path, 80, "10+12",  1, @types);
}

# subroutine that initializes and prepares the files to be aligned from a dataset
sub prepare_data {
	my ($path, $out_path, $ref, $species, $types) = @_;
	$path = "" if $path eq "/";
	foreach (@$types) {
		mkdir $out_path."/".$_, 0755 unless -e  $out_path.$_;
		my $in  = Bio::SeqIO->new(-file => "<".$path.$_.$ref.".fa" , '-format' => 'Fasta');
		while (my $seq = $in->next_seq) {
			my $name = $seq->id();
			$name .= "i" if (-e $out_path.$_."/".$name.".fa");
			my $out  = Bio::SeqIO->newFh(-file =>  ">".$out_path."/".$_."/".$name.".fa" , '-format' => "Fasta");
			$seq->id("$ref");
			print $out $seq;
		}
	}
	foreach my $animal (@$species) {
		undef my %control;
		foreach (@$types) {
			my $in  = Bio::SeqIO->new(-file => "<".$path.$_.$animal.".fa" , '-format' => 'Fasta');
			while (my $seq = $in->next_seq) {
				my $name = $seq->id();
				$name .= "i" if ($control{$name} && $control{$name} == 1);
				my $out  = Bio::SeqIO->newFh(-file =>  ">>".$out_path."/".$_."/".$name.".fa" , '-format' => "Fasta");
				$seq->id("$animal");
				print $out $seq;
				$control{$name} = 1;
			}
		}
	}
}

# subroutine that makes the calculations using clustalw and gibbs
sub process_data {
	my($dir, $cs_threshold, $g_lengths, $g_cutoff, @types) = @_;

	foreach (@types) {
		my $nuc_aa = "aa";
		$nuc_aa = "n" unless $_ eq "CDS";

		my $directory .= $dir."/".$_."/";

		my @in_files = <$directory*>;

		foreach my $file (@in_files) {
			my @name = split /\./, $file;
			my $out = $name[0]."al.".$name[1];

			# ClustalW
			clustalw($file, "Fasta", $out);
			my $out_seq = $name[0]."cs.".$name[1];
			seek_align($out,$out_seq,$nuc_aa,$cs_threshold);	# gets consensus sequence

			# Gibbs
			my $g_path = '/usr/local/bin/gibbs';	# Your location of Gibbs
			my $expected = 4 if $nuc_aa eq "n";	# Omit if loooking for aa motifs
			gibbs($g_path,$file,$g_lengths,$nuc_aa,$expected,$g_cutoff);
		}
	}
}

1;
