#!/usr/bin/perl
#
# Series of text parsers, warappers, and other utilities for use with ARCT. Basically, all the auxilliary modules are here.

use Bio::DB::GenPept;
use Bio::DB::GenBank;
use Bio::DB::RefSeq;
use Bio::SeqIO;
use Bio::Seq;
use Bio::SeqFeatureI;
use Bio::AnnotationI;

#get_gene($in_file, $in_format, $query, $out_file, $out_format);

# subroutine to get a gene or protein from a file with many sequences
sub get_gene  {
	my($in_file, $in_format, $query, $out_file, $out_format) = @_;

	my $seqio  = Bio::SeqIO->new( '-format' => "$in_format" , -file => "$in_file");
	my $out  = Bio::SeqIO->newFh(-file => ">$out_file" , '-format' => "$out_format");
	while ($seq = $seqio->next_seq()) {
		# For finding the query sequence according to accession number
		#if ($query =~ /.*[0-9].*/ && $seq->display_id eq "$query") {
		#	print "$query found: ".$seq->desc."\n";
		#	print $out $seq;
		#}

		# For finding the query sequence according to the description
		if ($seq->desc =~ /$query/) {
			# To change the display_id
			if ($query =~ /.*\[(\w+) (\w+).*\].*/) {
				$name = $1."_".$2;
			}
			print "$name found: ".$seq->desc."\n";
			$seq->display_id($name);
			print $out $seq;
		}
	}
}

# the next subroutine just puts all files in a directory together in the same file
sub join_dir {
	my ($path, $out_file) = @_;
	$path = "" if $path eq "/";
	my @files = <$path*>;
	print ">".$path.$out_file."\n";
	open TEMP, ">".$path.$out_file or die "Could not open output: $!\n";
	close TEMP;
	my $out  = Bio::SeqIO->newFh(-file =>  ">>$path".$out_file , '-format' => "Fasta");
	foreach (@files) {
		my $in  = Bio::SeqIO->new(-file => "<".$_ , '-format' => 'Fasta');
		while (my $seq = $in->next_seq) {
			my $name = $seq->id();
			$seq->id("$name");
			print $out $seq;
		}
		unlink $_;
	}
}

# This is the main parser of the report from profiles.pl. It takes a report and returns an hash with species as key and accession numbers, similarity, etc. as values.

# Subroutine to parse the comparison.txt report and return an array.
sub report_parser {
	open IN, "<".$_[0];
	my @output;

	while (<IN>) {
		if (/Reference gene:(.*)/) {
			my $tmp = $1;
			if ($tmp =~ /\s*(\w+), ([0-9]+) aa & (\w+)/) {
				my $ref_name = $1;
				my $ref_length = $2;
				my $ref_acc = $3;
				push @output, "$ref_name, $ref_length, $ref_acc";
			}
		} elsif (/Final.*/) {
			return @output;
		} elsif (/(\w+):(.*)/) {
			my $org = $1;
			my $tmp = $2;
			if ($tmp =~ /\s*(\w+), ([0-9]+) aa, Similarity = ([0-9.]+)%, Positives = ([0-9.]+).*/) {
				my $hit_acc = $1;
				my $hit_length = $2;
				my $hit_sim = $3;
				my $hit_pos = $4;
				$output[$#output] .= "&$org, $1, $2, $3, $4";
			}
		}
	}
	return @output;	# Just in case the Final Results are not in the file
}

# Small subroutine to get the sequence upstream of the translation start site. Since it only gets the sequence in the annotated databases, it will most likely get <1000 bp. The way it works is by getting the sequence upstream of the translation start site, which can be large, small, or non-existent--depending on the annotation. Please note that errors may emerge from using the subroutine and so use carefully.
sub get_promoter {
	my $filename = shift @_;
	my @query = @_;

	my $out  = Bio::SeqIO->newFh(-file => ">$filename" , '-format' => 'Fasta');
	foreach $query (@query) {

		if ( $query =~ /NM(\w)+/ ) {
			$database = new Bio::DB::RefSeq;
		} else {
			$database = new Bio::DB::GenBank;
		}
		my $seq = $database->get_Seq_by_acc("$query");

		foreach $feat ( $seq->top_SeqFeatures() ) {
			foreach $tag ( $feat->primary_tag() ) {
				if ( $tag eq 'CDS' ) {
					$start = $feat->start;
					$promo = 1;
					print "CDS found from ", $start, " to ", $feat->end, " \n";
					$promo = $start - 999 if ($start > 1000);
				}
			}
		}

		$promoter=$seq->subseq($promo,$start-1);

		$seqobj = Bio::Seq->new( -display_id => $seq->display_id,
						-desc => $seq->desc,
						-seq => $promoter);
		print $out $seqobj;
		print "Saved upstream sequence (".length($promoter)." bp) of ".$seq->display_id." (".$seq->accession_number.")\n";
		print "Description: ".$seq->desc()."\n";
	}
}

# Subroutine to download protein or gene sequences from databases
sub get_single {
	my ($start, $format,  $out_file, $db, $query_in) = @_;
	my @query = split / /, $query_in;

	if ($start > 1) {
		$out_file = ">"."$out_file";
		foreach (1 .. ($start-1)) {
			shift @query;
		}
	}

	my $out  = Bio::SeqIO->newFh(-file => ">$out_file" , '-format' => "$format");

	foreach my $query_id ( @query ) {
		# Access RefSeq, GenBank, or GenPept (other databases are possible, of course, just check the Bioperl documentation for a list of available databases)
		if ($db == 1) {	# User wants genes
			if ( $query_id =~ /N(M|P)(\w)+/) {
				$database = new Bio::DB::RefSeq;
			} else {
				$database = new Bio::DB::GenBank;
			}
		} else {	# If you're dealing with proteins, you'll have to use GenPept instead of GenBank
			if ( $query_id =~ /N(M|P)(\w)+/) {
				$database = new Bio::DB::RefSeq;
			} else {
				$database = new Bio::DB::GenPept;
			}
		}

		# Retrieve the sequence by accession number
		my $seq = undef;
		my $retries = 5;
		my $tmp = 0;
		until ( defined($seq) || $tmp >= $retries ) {
			$seq = $database->get_Seq_by_acc("$query_id");
			++$tmp;
		}
		die "Could not access database for $query_id\n" unless (defined($seq));

		# Extract the protein from the entire record
		my @prot=0;
		foreach $feat ( $seq->all_SeqFeatures() ) {
			foreach $tag ( $feat->all_tags() ) {
				if ( $tag eq 'translation' ) {
					@prot = $feat->each_tag_value($tag);
				}
				# If you are downloading genes and want to use the protein's accession number, you can use
				#if ( $tag eq 'protein_id' ) {
				#	@id = $feat->each_tag_value($tag);
				#}
			}
		}

		# We only keep the protein sequence, accession number and ID
		if ($db == 1) {	# Even if we are downloading genes we just keep the protein
			$newseq = $prot[0];
		} else {
			$newseq = $seq->seq();
		}
		my $seqobj = Bio::Seq->new( -display_id => $seq->display_id,
								-accession_number => $seq->accession_number,
								-seq => $newseq);

		# Save results to file
		print $out $seqobj;
		++$progress;
		print "Progress: $progress of ".($#query+1)."\n";
	}
}

# Batch mode (recommended unless we are using a small number of different types of sequences)
sub get_batch {
	my ($start, $format,  $out_file, $db, $query) = @_;

	if ($db == 1) {
		$database = new Bio::DB::RefSeq(-retrievaltype => 'tempfile');
	} elsif ($db == -1) {
		$database = new Bio::DB::GenPept(-retrievaltype => 'tempfile');
	} else {
		$database = new Bio::DB::GenBank(-retrievaltype => 'tempfile');
	}

	my $seqio = $database->get_Stream_by_id($query);

	if ($start > 1) {
		$out_file = ">"."$out_file";
		foreach (1 .. ($start-1)) {
			$seqio->next_seq;
		}
	}

	my $out  = Bio::SeqIO->newFh(-file => ">$out_file" , '-format' => "$format");

	while( my $seq  =  $seqio->next_seq ) {
		# Extract the protein
		my @prot=0;
		foreach $feat ( $seq->all_SeqFeatures() ) {
			foreach $tag ( $feat->all_tags() ) {
				if ( $tag eq 'translation' ) {
					@prot = $feat->each_tag_value($tag);
				}
			}
		}

		if (defined($prot[0])) {
			$newseq = $prot[0];
		} else {
			$newseq = $seq->seq();
		}

		# Save results to file as FASTA
		my $seqobj = Bio::Seq->new( -display_id => $seq->display_id(),
					-accession_number => $seq->accession_number(),
					-seq => $newseq);
		print $out $seqobj;
	}
}

# Subroutine to get promoters from EPD (http://www.epd.isb-sib.ch/). The $query can either be an accession number
# or a HUGO symbol. If a promoter is found at the EPD, the subroutine returns the EPD accession plus the -499 to
# 100 bp of sequence. See bellow on how to adjust these parameters.
# Examples:
# print get_prom("NM_021081");
# print get_prom("GHR");
#
sub get_prom {

	my ($query) = @_;


	my $URL ='http://www.epd.isb-sib.ch/cgi-bin/epd_quick_search.pl?query_str='.$query.'&database=epd';


	my $retrieve = get($URL);

	my $tree = HTML::TreeBuilder->new_from_content($retrieve);

	$formatter = HTML::FormatText->new(leftmargin => 0, rightmargin => 50);

 	$text = $formatter->format($tree);



	my @lines = split/\n/, $text;

	my $acc = "";

	foreach (@lines) {

		if (/(HS_\w+)\s*/) {

			$acc = $1;

		}

	}



	if ($acc ne "") {
		# You can adjust the sequence to retrieve in here

		$URL = 'http://www.epd.isb-sib.ch/cgi-bin/query_result.pl?out_format=FASTA&from=-499&to=100&Entry_0='.$acc;

		my $retrieve = get($URL);

		my @cut = split/\n/, $retrieve;

		my $seq = "";

		foreach (@cut) {

			if (/^[A|C|G|T|N]/g) {

				s/\s//g;

				$seq .= $_;

			}

		}

		return $acc, $seq;

	}
}

# Subroutine that creates an hash with the Gene Ontologies from the OBO format flat file obtained
# from http://geneontology.org/. It returns an hash with the GO and, if you pass a reference to
# an hash, that hash will feature the type for each GO.
# Example:
# my %go_type = parse_go("gene_ontology.obo", \%type);
sub parse_go {
	my ($file, $hash) = @_;
	my ($go, $name, $type) = "";
	my %results = ();

	open IN, "<$file" or die("Could not open $file\n");
	while ( <IN> ) {
		if ( /\[Term\]/ ) {
			if ( $go > 0) {	# Define the type of GO to be saved
				$results{"$go"} = $name;
				if ( $type eq "biological_process" ) {
					$$hash{$go*1} = "P";
				} elsif ( $type eq "molecular_function" ) {
					$$hash{$go*1} = "F";
				} elsif ( $type eq "cellular_component" ) {
					$$hash{$go*1} = "C";
				}
			}
			$go = $1;
		} elsif ( /^name: (.*)/ ) {
			$name = $1;
		} elsif ( /^namespace: (\w+)/ ) {
			$type = $1;
		} elsif ( /^id: GO:(\d+)/ ) {
			$go = $1;
		}

	}
	return %results;
}

# Very simple subroutine to find a gene's (or protein's) position in a file with many sequences.
sub parse_gene {
	my($in_file, $in_format, $query) = @_;
	my $tmp = 0;

	my $seqio  = Bio::SeqIO->new( '-format' => "$in_format" , -file => "$in_file");
	while ($seq = $seqio->next_seq()) {
		++$tmp;
		if ($seq->display_id eq "$query") {
			print "$query found at position $tmp\n";
		}
	}
}

# subroutine to parse the mtDNA files
sub parse_mtdna {
	my ($out_file, $accession) = @_;
	my $database = new Bio::DB::RefSeq;
	my $seq = $database->get_Seq_by_id("$accession");
	open TRNA, ">tRNA".$out_file.".fa";
	open RRNA, ">rRNA".$out_file.".fa";
	open CDS, ">aa_". $out_file.".fa";
	my $promo_end = 16000;
	my $promo_last = 0;
	foreach $feat ( $seq->top_SeqFeatures() ) {
		foreach $tag ( $feat->primary_tag() ) {
			if ( $tag eq 'tRNA' ) {
				$start = $feat->start;
				$end = $feat->end;
				print "tRNA found from ", $start, " to ", $end, " \n";
				my @x = $feat->each_tag_value('product');
				print TRNA ">".$x[0]."\n";
				print TRNA $seq->subseq($start,$end)."\n";
				$promo_end = $start unless ($start > $promo_end);
				$promo_last = $end unless ($end < $promo_end);
			}
			if ( $tag eq 'rRNA' ) {
				$start = $feat->start;
				$end = $feat->end;
				print "rRNA found from ", $start, " to ", $end, " \n";
				my @y = $feat->each_tag_value('product');
				print RRNA ">".$y[0]."\n";
				print RRNA $seq->subseq($start,$end)."\n";
			}
			if ( $tag eq 'CDS' ) {
				$start = $feat->start;
				$end = $feat->end;
				print "CDS found from ", $start, " to ", $end, " \n";
				my @z = $feat->each_tag_value('gene');
				print CDS ">".$z[0]."\n";
			}
			undef $tag;
			foreach $tag ( $feat->all_tags() ) {
				if ( $tag eq 'translation' ) {
					my @w = $feat->each_tag_value('translation');
					print CDS $w[0]."\n";
				}
			}
		}
	}

	open PROMOTER, ">promo".$out_file.".fa";
	print PROMOTER ">promoter\n";
	$temp = $seq->length();
	print PROMOTER $seq->subseq($promo_last, $temp);
	print PROMOTER $seq->subseq(1,$promo_end) unless ($promo_end < 2);

	open OUT, ">mtDNA".$out_file.".txt";
	print OUT ">mtDNA\n";
	print OUT $seq->seq();
}

1;
