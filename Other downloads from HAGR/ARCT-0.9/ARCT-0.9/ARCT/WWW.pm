#!/usr/bin/perl -w
#
# WWW subroutines used at the Human Ageing Genomics Resources (http://genomics.senescence.info)

# I-CpGI or Identifier of CpG Islands

sub i_cpgi {
	my($cg, $length, $ratio, $dna) = @_;
	my $dna_length = length($dna);
	die "Invalid sequence. Please verify your query sequence.\n" unless $dna_length > $length;
	undef my %islands;
	my $cpg_end = 0;
	my $cpg_start = 0;

	# following the method by Takai and Jones, what follows are steps A to D
	while ($cpg_start != -1) {
		my $cursor = 0;
		$cursor = $cpg_end + 1 unless ($cpg_end == 0);
		$cpg_start = -1;

		my $temp = 0;
		# we start checking in the 3' direction
		while (($dna_length > $cursor + $length) && ($temp != -1)) {
			my $seq = substr($dna, $cursor, $length);
			my $c_content = ($seq =~ s/c/c/gi)/$length;
			my $g_content = ($seq =~ s/g/g/gi)/$length;
			if ((($c_content + $g_content) > $cg) && (find_cpg($c_content, $g_content, $ratio, $seq) == 1)) {
				$cpg_start = $cursor unless ($cpg_start != -1);
				$cpg_end = $cursor + $length;
				$cursor += $length;
			} else {
				++$cursor;
				$temp = -1 unless ($cpg_start == -1);
			}
		}

		# and then step E, seeking in the 5' direction to find the end of the island
		$cursor = $cursor + $length;
		$cursor = $dna_length unless ($cursor < $dna_length);

		while ($cpg_end <= $cursor && $cpg_start != -1) {
			my $seq = substr($dna, $cursor-$length, $length);
			my $c_content = ($seq =~ s/c/c/gi)/$length;
			my $g_content = ($seq =~ s/g/g/gi)/$length;
			if ((($c_content + $g_content) > $cg) && find_cpg($c_content, $g_content, $ratio, $seq)) {
				$cpg_end = $cursor;
				$cursor -= $length;
			}
			$cursor--;
		}

		# now for steps G and H
		my $flag = 0;
		while ($flag == 0 && $cpg_start != -1) {
			my $seq = substr($dna, $cpg_start, $cpg_end - $cpg_start);
			my $c_content = ($seq =~ s/c/c/gi)/length($seq);
			my $g_content = ($seq =~ s/g/g/gi)/length($seq);
			if ((($c_content + $g_content) > $cg) && find_cpg($c_content, $g_content, $ratio, $seq)) {
				$flag = 1;
			} else {
				$cpg_start++;
				$cpg_end--;
			}
		}

		$islands{"$cpg_start"} = "$cpg_end" unless $cpg_start == -1;
	}

	# now for step I
	my $past_end = 0;
	my $past_start = 0;
	foreach (sort by_number keys %islands) {
		if ($past_end + $length/2 > $_ && $past_end != 0) {
			$islands{$past_start} = $islands{$_};
			$past_end = $islands{$_};
			delete $islands{$_};
		} else {
			$past_end = $islands{$_};
			$past_start = $_;
		}
	}

	# graphical display
	print "I-CpGI\n\n";
	printf "%g%72g\n", 0, $dna_length;
	print "|----------------------------------------------------------------------|\n";
	foreach my $cp_start (sort by_number keys %islands) {
		my $start = int($cp_start/$dna_length*70);
		my $end = int($islands{$cp_start}/$dna_length*70);
		for ($start..$end) {
			printf "%s", "^";
		}
		print "\n";
		printf "%".$start."s%".$end."s\n", $cp_start, $islands{$cp_start};
	}
	print "\n".scalar(keys %islands)." CpG island(s) found\n";
}

sub find_cpg {
	my($c_ratio, $g_ratio, $ratio, $query_seq) = @_;
	my $observed = ($query_seq =~ s/cg/cg/gi);
	my $expected = $c_ratio*$g_ratio*length($query_seq);
	return 0 unless ($c_ratio > 0 && $g_ratio > 0);
	if ($observed/$expected >= $ratio) {
		return 1;
	} else {
		return 0;
	}
}

# Interaction Graphic Display
#
# subroutine to draw protein-protein interactions

sub draw_interactions {
	my($core, $inter, $names, $type) = @_;
	my $fix = 0;
	$fix = 1 if ($core != 0);
	undef my %coordx;
	undef my %coordy;
	my $radius = 50+(scalar(keys %$names)-$fix)*10;
	my $im = new GD::Image(3*$radius,3*$radius);
	my $brush = new GD::Image(2,2);
	my $ct_x = 3*$radius/2;
	my $ct_y = 3*$radius/2;
	my $white = $brush->colorAllocate(255,255,255);
	my $blue = $brush->colorAllocate(0,0,255);
	my $red = $brush->colorAllocate(255,0,0);
	my $green = $brush->colorAllocate(0,255,0);
	my $black = $brush->colorAllocate(0,0,0);
	$brush->transparent($white);
	$brush->fill(1,1,$blue);
	$im->setBrush($brush);
	$im->transparent($white);

	if ($fix == 1) {
		my $colour = $black;
		$colour = $red if ($$type{"$core"} eq "human");
		$colour = $green if ($$type{"$core"} eq "model");
		$colour = $blue if ($$type{"$core"} eq "mammal");
		if ($$type{"$core"} eq "downstream") {
			$im->string(gdSmallFont,$ct_x-length($$names{$core})*2.75,$ct_y-8,$$names{$core},$colour);
		} else {
			$im->string(gdMediumBoldFont,$ct_x-length($$names{$core})*3.5,$ct_y-8,$$names{$core},$colour);
		}
		$cordx{$core} = $ct_x;
		$cordy{$core} = $ct_y;
		$im->arc($ct_x,$ct_y,50,50,0,360,$black);
	}

	my $divider = 6.2832/(scalar(keys %$names)-$fix);
	my $counter = 0;
	foreach (keys %$names) {
		next if ($fix == 1 && $_ == $core);
		++$counter;
		my $x = $ct_x+int(cos($divider*$counter)*$radius)+length($$names{$_})*3.5;
		my $y = $ct_y+int(sin($divider*$counter)*$radius)+8;
		my $colour = $black;
		$colour = $red if ($$type{"$_"} eq "human");
		$colour = $green if ($$type{"$_"} eq "model");
		$colour = $blue if ($$type{"$_"} eq "mammal");
		if ($$type{"$_"} eq "downstream") {
			$im->string(gdSmallFont,$x-length($$names{$_})*2.75,$y-8,$$names{$_},$colour);
		} else {
			$im->string(gdMediumBoldFont,$x-length($$names{$_})*3.5,$y-8,$$names{$_},$colour);
		}
		$cordx{$_} = $x;
		$cordy{$_} = $y;
		$im->arc($x,$y,50,50,0,360,$black);
	}
	foreach (keys %$inter) {
		my @temp = split/\+/, $$inter{$_};
		my $xa = $cordx{$temp[0]};
		my $ya = $cordy{$temp[0]};
		my $xb = $cordx{$temp[1]};
		my $yb = $cordy{$temp[1]};
		my $distance = sqrt(($xa-$xb)*($xa-$xb)+($ya-$yb)*($ya-$yb));
		if ($fix == 1) {
			$im->line($xa+int(($xb-$xa)/$distance*25),$ya+int(($yb-$ya)/$distance*25),$xb+int(($xa-$xb)/$distance*25),$yb+int(($ya-$yb)/$distance*25),gdBrushed);
		} else {
			$im->line($xa+int(($xb-$xa)/$distance*25),$ya+int(($yb-$ya)/$distance*25),$xb+int(($xa-$xb)/$distance*25),$yb+int(($ya-$yb)/$distance*25),$black);
		}
	}
	print "Content-type:image/png\n\n";
	binmode(STDOUT);
	print $im->png;
}

# Phylogenetic Tree Plotter

sub draw_phylum {
	my ($divergence, $draw_common, $draw_names, $tree_type, @phylogeny) = @_;
	my $count = 0;
	my (%done, %past, %present);
	my $previous_div = 0;
	foreach ((keys %{$phylogeny[5]})) {
		$done{$_} = $phylogeny[5]->{$_}." ".$phylogeny[6]->{$_};
	}
	my $im = new GD::Image(800+$draw_common*100,100+60*$divergence);
	my $white = $im->colorAllocate(255,255,255);
	my $black = $im->colorAllocate(0,0,0);
	$im->transparent($white);
	for $present (@phylogeny) {
		++$count;
		my %hagrid = reverse %$present;
		my @names = sort sort_alpha keys %hagrid;

		if (scalar(@names) > $previous_div && $previous_div > 0) {
			for (1..scalar(@names)) {
				my $num_nam = $_-1;
				undef my @old;
				undef my @new;
				foreach my $key ( keys %$present) {
					if ($$present{$key} eq $names[$num_nam]) {
						push @new, $key;
					}
				}
				my $backref = $past{"$new[0]"};
				foreach my $key ( keys %$present) {
					if ($past{$key} eq $backref) {
						push @old, $key;
   					}
				}
				if ((@old) != (@new) && defined($new[0])) {
					undef my %branch;
					foreach (@old) {
						$branch{"$$present{$_}"} += 1;
					}
					my $counter = 0;
					my $div_conv = 0;
					my $fix_count = 1;
					foreach (keys %branch) {
						++$counter;
	       					my $branching = $branch{$_}-1;
 						$div_conv = $div_conv+$counter*$fix_count;
						$fix_count = -$fix_count;
						foreach my $key (@old) {
							if ($$present{$key} eq $_) {
								$past{$key} += $div_conv-$branching*$fix_count;
								$im->line((($count-1)*100),(40*$divergence+($past{$key}-$div_conv+$branching*$fix_count)*20),(($count-1)*100),(40*$divergence+$past{$key}*20),$black);
								if ($branching == 0 && $tree_type == 1) {
									$im->line((($count-1)*100),(40*$divergence+$past{$key}*20),(($count)*100),(40*$divergence+$past{$key}*20),$black);
					       			       	$im->string(gdSmallFont,($count*100+10),(40*$divergence+$past{$key}*20-5),$done{$key},$black);
									$done{$key} = -1;
								}
								delete $$present{$key};
   							}
						}
				   	}

					my $hagrid_key = $hagrid{$names[$num_nam]};
					if ($done{$hagrid_key} != -1 && $names[$num_nam] =~ /^[A-Z]/) {
						$im->line((($count-1)*100),(40*$divergence+$past{$hagrid_key}*20),(($count)*100),(40*$divergence+$past{$hagrid_key}*20),$black);
						$im->string(gdSmallFont,(($count-1)*100+3),(40*$divergence+$past{$hagrid_key}*20),$names[$num_nam],$black) unless ($draw_names == 0);
					}
				} else {
					my $hagrid_key = $hagrid{$names[$num_nam]};
					if ($done{$hagrid_key} != -1 && $names[$num_nam] =~ /^[A-Z]/) {
						$im->line((($count-1)*100),(40*$divergence+$past{$hagrid_key}*20),(($count)*100),(40*$divergence+$past{$hagrid_key}*20),$black);
						$im->string(gdSmallFont,(($count-1)*100+3),(40*$divergence+$past{$hagrid_key}*20),$names[$num_nam],$black) unless ($draw_names == 0);
					}
				}
			}
		} else {
			if ($draw_common == 1) {
				for (1..scalar(@names)) {
					my $hagrid_key = $hagrid{$names[$_-1]};
					$im->line((($count-1)*100),(40*$divergence+$past{$hagrid_key}*20),(($count)*100),(40*$divergence+$past{$hagrid_key}*20),$black) unless ($done{$hagrid_key} == -1);
					$im->string(gdSmallFont,(($count-1)*100+3),(40*$divergence+$past{$hagrid_key}*20),$names[$_-1],$black) unless ($draw_names == 0 || $done{$hagrid_key} == -1);
				}
			} else {
				--$count;
			}
		}
		$previous_div = scalar(@names);
	}
	if ($tree_type != 1) {
		foreach (keys %past) {
			$im->line((($count-1)*100),(40*$divergence+$past{$_}*20),(($count)*100),(40*$divergence+$past{$_}*20),$black);
			$im->string(gdSmallFont,($count*100+10),(40*$divergence+$past{$_}*20-5),$done{$_},$black);
		}
	}
	print "Content-type:image/png\n\n";
	binmode(STDOUT);
	print $im->png;
}

sub sort_alpha { "\L$a" cmp "\L$b" }

sub by_number { $a <=> $b }

# my death subroutine, in case something goes wrong
sub death {
	print "Content-type:text/plain\n\n";
	print "Error in script!\nPlease make sure you have entered your query correctly.\nIf the error persists, please contact the system administrator at feedback\@genomics.senescence.info\n";
	exit;
}

1;
