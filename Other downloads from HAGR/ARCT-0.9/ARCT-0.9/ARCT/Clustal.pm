#!/usr/bin/perl
#
# Module to use ClustalW locally and analyse alignments

# This is the path where you MUST have clustalw (tested with ClustalW 1.8.3)
$ENV{CLUSTALDIR} = '/usr/local/bin';
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::SeqIO;

# Subroutine to align sequences using a locally installed version of ClustalW
sub clustalw {
	my($in_file, $format, $out_file) = @_;
	undef my @seq_array;

	# We get the sequences from file
	my $in  = Bio::SeqIO->new(-file => "<$in_file" , '-format' => "$format");
	while (my $seq = $in->next_seq) {
		push @seq_array, $seq;
	}

	#  Build a clustalw alignment factory
	my @params = ('ktuple' => 2, 'matrix' => 'BLOSUM');	# You can specify'OUTFILE' => 'test.txt' but its safer to work with AlignIO
	my $factory = Bio::Tools::Run::Alignment::Clustalw->new(@params);

	my $seq_array_ref = \@seq_array;
	# where @seq_array is an array of Bio::Seq objects
	my $aln = $factory->align($seq_array_ref);

	# At time of writing, we must save the alignment using AlignIO (to save as ClustalW 1.81). Otherwise the format (from ClustalW 1.83) will be incompatible with the other applications in Bioperl
	my $out  = Bio::AlignIO->new(-file => ">$out_file",
				'-format' => 'clustalw');
	$out->write_aln($aln);
}

# Subroutine to align alignments. Beware that this subroutine eliminates the first sequence from each alignment but the first one to make sure only one homo sequence is present
# at present it only aligns alignments to either alignments or sequences, not a combination of both. Also, one sequence/alignment per file, please
sub combine_align {
	my($out_file, $ref_file, $input_format, @files) = @_;
	undef my @queries;
	my $input;
	my $update;
	my $new_align;

	my $ref_in  = Bio::AlignIO->new(-file => "<$ref_file",
				'-format' => "clustalw");
	$new_aln = $ref_in->next_aln();

	my @params = ('ktuple' => 2, 'matrix' => 'BLOSUM');	# You can specify'OUTFILE' => 'test.txt' but its safer to work with AlignIO
	my $factory = Bio::Tools::Run::Alignment::Clustalw->new(@params);

	# We get the sequences/alignments from file
	foreach my $count (1..scalar(@files)) {
		if ($input_format =~ /^fasta$/i) {
			$input  = Bio::SeqIO->new(-file => "<$files[$count-1]" , '-format' => "Fasta");
			$update = $input->next_seq();
		} elsif ($input_format =~ /^clustalw$/i) {

			$input  = Bio::AlignIO->new(-file => "<$files[$count-1]",
						'-format' => "clustalw");
			my $aln = $input->next_aln();
			my $seq = $aln->get_seq_by_pos(1);	# rmove first sequence since this should be a repeated human sequence
			$aln->remove_seq($seq);
			$update = $aln;
		} else {
			die "Error. Please check your options and paths\n";
		}
		$new_aln = $factory->profile_align($new_aln,$update);
	}

	my $out  = Bio::AlignIO->new(-file => ">$out_file",
				'-format' => 'clustalw');
	$out->write_aln($new_aln);
}

# Subroutine to analyze ClustalW alignments. You can easily modify it to accept other formats.
sub seek_align {
	my($in_file, $out_file, $type, $consensus, $iupac) = @_;
	$type = "bp" unless $type eq "aa";

	# We get the alignment from file
	my $in  = Bio::AlignIO->new(-file => "<$in_file",
				'-format' => 'clustalw');
	my $aln = $in->next_aln();

	open OUT, ">$out_file";

	# We get some stats from it
	print OUT "Alignment length: ".$aln->length." $type\n";
	print OUT "Average identity: ".int($aln->percentage_identity)."%\n";
	#print OUT "Match line:\n".$aln->match_line()."\n";

	# the $iupac flag is set to get IUPAC ambiguity codes while comparing nucleotides
	if ($iupac && $iupac == 1) {
		print OUT $aln->consensus_iupac();
	} else {
		print OUT $aln->consensus_string("$consensus")."\n";	# Gives the consensus sequence found in at least 80% of aligned sequences. This can be very useful for comparing different alignments and you can create alignments of consensus sequences or aligments of aligments.
	}
}

1;
