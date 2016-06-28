# collect.pl summarizes the output from the GANN core software.
# Copyright (C) 2000-2002 Robert G. Beiko

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.


# CollectScores.pl summarizes the results of a GANN run.
#	Required arguments:
#	$ARGV[0] -> Name of index file
#	$ARGV[1] -> Name of score file
#	Output files:
#	indexList.csv -> the index names, in order
#	scoresum.csv -> a summary of all scores for each chromosome that contains a given index
#	paramsum.csv -> a summary of neural network parameters over all generations

$| = 1;

if (scalar(@ARGV) != 3) { die "Usage: perl collect.pl [indexFile] [scoreFile] [num_rounds]"; }

$indexFile = $ARGV[0];
$scoreFile = $ARGV[1];
$num_rounds = $ARGV[2];

$outList = "indexList.csv";

$nonIndCols = 3; # Number of leading columns that do not have indices (i.e. PosNeg)

# ------------------------------ Create indexList.csv ---------------------------------------

print STDOUT "Creating indexList.csv...";

open (INDEX, "$indexFile") or die "Could not open $indexFile\n";

$index_count = <INDEX>; chomp $index_count;
$set_list = <INDEX>; chomp $set_list;
@set_sizes = split (/,/,$set_list);
$header_list = <INDEX>; chomp $header_list;
@headers = split (/,/,$header_list);

print STDOUT "Number of indices: $index_count\n";

open (OUT, ">$outList");
for $i ($nonIndCols...$index_count + $nonIndCols) {
	print OUT "$headers[$i]\n";
}
close (OUT);

print STDOUT "Done.\n";

# ----------------------- Extract the score from each OGA round -----------------------------

print STDOUT "Extracting OGA chromosome scores...";

open (INSCORE, "$scoreFile") or die "Could not open $scoreFile\n";
@allscores = <INSCORE>;
close (INSCORE);

$add = 1; # Add = 1: add line to array

# Find the score lines sandwiched between lines of text, and record them
foreach $line (@allscores) {
	if ($line =~ m/C/) {
		$add = 1 - $add; # Switch!
	} else {
		if ($add == 1) {
			chomp $line;
			push (@score_array, $line);
		}
	}
}
print STDOUT "Done.\n";

# ----------------- Identify which indices were used in each OGA chromosome ----------------

print STDOUT "Identifying indices...";

open (INPARAM, "nnparams.csv") or die "Can't find nnparams.csv\n";
@allparams = <INPARAM>;
close (INPARAM);

$found = 0; $i = 0;
while (!$found) {
	if ($allparams[$i] =~ m/C/) { # Find the start point from the end of the headers
		@col_throw = split ',', $allparams[$i];
		$col_start = $#col_throw + 2;
		$found = 1;
	} else { ++$i; } 
}
$found = 0; $i = 0;
while (!$found) {
	if (!($allparams[$i] =~ m/C/)) { # Find the end point
		@col_throw = split ',', $allparams[$i];
		$col_stop = $#col_throw - 1;
		$found = 1;
	} else { ++$i; }
}
$total_input = $col_stop - $col_start + 1;
print STDOUT "Number of inputs: ", $total_input, "\n";

$max = 0; # Max input value

foreach $line (@allparams) { # Count the occurrences of each index in each OGA round
	if (!($line =~ m/C/)) {
		$info = $line;
		chomp $info;
		@bits = split /,/, $info;
		$OGAr = $bits[0];
		for $j ($col_start..$col_stop) {
			$inp = $bits[$j];
			if ($inp > $max) {
				$max = $inp;
			}
			++$counts[$OGAr][$inp];
		}
	}
}

open (COUNTOUT, ">input-sum.csv");

print COUNTOUT "Inp#,";
for $printRound (0..$#counts) {
	print COUNTOUT $printRound, ",";
}
print COUNTOUT "\n";

for $row_out (0..$max) {
	print COUNTOUT "$headers[$row_out + $nonIndCols] ($row_out),";
	for $round_out (0..$#counts) {
		print COUNTOUT "$counts[$round_out][$row_out],";
	}
	print COUNTOUT "\n";
}
close (COUNTOUT);

print STDOUT "Done.\n";

# ------------------------ BUILD THE SCORE SUMMARY FILE --------------------------------

print STDOUT "Building score summary...";

open (IN, "nnparams.csv") or die "Can't find nnparams.csv\n";
$chr_num = 0;
while (<IN>) { # Get the input identifiers
	chomp $_;
	if (!($_ =~ m/C/)) {
		@thisline = split /,/,$_;
		if ($thisline[0] == $num_rounds - 1) { # Count goes from 0
			for $i ($col_start..$col_stop) {
				$scores[$chr_num][$i  - $col_start][0] = $thisline[$i];
			}
			++$chr_num;
		}
	}
}
close (IN);

$chr_num  = 0;
foreach $line (@score_array) { # Collect the set of scores associated with each replicate...
	@scoreline = split /,/,$line;
	if ($scoreline[0] == $num_rounds - 1) {
		push(@{$score_set[$scoreline[1]]}, $scoreline[-1]);
	}
}
for $eachSet (0..$#score_set) { # ...then
	for $i (0..$total_input) {
		$scores[$eachSet][$i][1] = MAX (@{$score_set[$eachSet]});
	}
}

# Convert the inputs and scores to a score summary file

$outscore = "scoresum.csv";

$max_index = 0;
# Store scores for each input
for $i (0..$#scores) {
	for $j (0..$total_input - 1) {
		push (@{$scoresum[$scores[$i][$j][0]]}, $scores[$i][$j][1]);
		if ($scores[$i][$j][0] > $max_index) { $max_index = $scores[$i][$j][0]; }
	}
}


open (OUT, ">$outscore");

# Average the scores for each input
for $i (0..$max_index) {
	$avg[$i] = AVERAGE(@{$scoresum[$i]});
	$smax[$i] = MAX(@{$scoresum[$i]});
	$smin[$i] = MIN(@{$scoresum[$i]});
	print OUT "$i,$headers[$i + $nonIndCols],$smax[$i],$avg[$i],$smin[$i],,";
	@sortedscores = sort {$b <=> $a} @{$scoresum[$i]};
	foreach $bit (@sortedscores) { printf OUT ("%f",$bit); print OUT ","; }
	print OUT "\n";
}

close (OUT);

print STDOUT "Done.\n";

# -------------------------- Summarize the parameters in nnparams.csv ------------------------

print STDOUT "Summarizing nnparams.csv...";

$outfile = "paramSum.csv";

# Return the average value of an array
sub AVERAGE {
	if (scalar (@_) == 0) { return -1; }
	$sum = 0;
	foreach $elem (@_) {
		$sum += $elem;
	}
	return $sum / scalar (@_);
}

# Find the max value of an array (from Christiansen, p.113)
sub MAX {
	my $max = $_[0];
	foreach $elem (@_) {
		if ($elem > $max) {
			$max = $elem;
		}
	}
	return $max;
}

# Find the min value of an array (see above)
sub MIN {
	my $min = $_[0];
	foreach $elem (@_) {
		if ($elem < $min) {
			$min = $elem;
		}
	}
	return $min;
}

sub STDEV {
	my $sum = 0;
	my $sumsq = 0;
	my $n = scalar (@_);

	foreach $elem (@_) {
		$sum += $elem;
		$sumsq += $elem * $elem;
	}
	return sqrt(($sumsq - $sum * $sum / $n) / ($n - 1));
}

open (INPARAMSUM, "nnparams.csv") or die "Can't open $infile\n";

@allsum = <INPARAMSUM>;

close(INPARAMSUM);

$found = 0; $i = 0;
while (!($found)) {
	if ($allsum[$i] =~ m/C/) {
		$header = $allsum[$i];
		$found = 1;
	} else { ++$i; }
}

foreach $line (@allsum) {
	if (!($line =~ m/C/)) {
		$lbuf = $line;
		chomp $lbuf;
		@lsplit = split /,/, $lbuf;
		push (@{$vals[$lsplit[0]]}, [ @lsplit ]);
	}
}

$cols = $#{$vals[0][0]}; # Number of columns to consider
$ocs = $#{$vals[0]}; # Number of generations

for $i (0..$#vals) {
	for $j (0..$cols) {
		@val_list = ();
		for $getval (0..$ocs) {
			push (@val_list, $vals[$i][$getval][$j]);
		}
		push (@{$AVG[$i]}, AVERAGE @val_list);
		push (@{$MAX[$i]}, MAX @val_list);
		push (@{$MIN[$i]}, MIN @val_list);
		push (@{$STD[$i]}, STDEV @val_list);
	}
}

open (OUTPARAM, ">$outfile");

print OUTPARAM "$header\n\n";

print OUTPARAM "Averages:\n";
for $q (0..$#vals) {
#	print OUTPARAM "$q,";
	foreach $column (@{$AVG[$q]}) {
		print OUTPARAM "$column,";
	}
	print OUTPARAM "\n";
}
print OUTPARAM "\n";

print OUTPARAM "Max:\n";
for $q (0..$#vals) {
#	print OUTPARAM "$q,";
	foreach $column (@{$MAX[$q]}) {
		print OUTPARAM "$column,";
	}
	print OUTPARAM "\n";
}
print OUTPARAM "\n";

print OUTPARAM "Min:\n";
for $q (0..$#vals) {
#	print OUTPARAM "$q,";
	foreach $column (@{$MIN[$q]}) {
		print OUTPARAM "$column,";
	}
	print OUTPARAM "\n";
}
print OUTPARAM "\n";

print OUTPARAM "Stdev:\n";
for $q (0..$#vals) {
#	print OUTPARAM "$q,";
	foreach $column (@{$STD[$q]}) {
		print OUTPARAM "$column,";
	}
	print OUTPARAM "\n";
}
print OUTPARAM "\n";

close(OUTPARAM);

print STDOUT "Done.\n";

