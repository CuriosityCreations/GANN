# whittle.pl extracts a subset from a given index file.
# Copyright (C) 2000-2002 Robert G. Beiko

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

$| = 1;

$USAGE = "Usage: perl whittle.pl [ infile ] [ -r range | -f file ] [ outfile ]";

$numNonInd = 3;

sub AVERAGE (@) {
	my $size = scalar(@_);
	my $sum = 0;
	if ($size == 0) { return 0; }
	foreach $val (@_) {
		$sum += $val;
	}
	return $sum / $size;
}

if (scalar(@ARGV) != 4) { die "$USAGE\n"; }

$infile = $ARGV[0];
$inflag = $ARGV[1];
$inlist = $ARGV[2];
$outfile = $ARGV[3];

open (IN, "$infile") or die "Can't open $infile\n";
$header = <IN>;
$sizes = <IN>;
@everything = <IN>;
$last_seq = $#everything;
$i = 0;
foreach $line (@everything) {
	chomp $line;
	@allind = split (/,/, $line);
	@{$allsplit[$i++]} = @allind;
}
close(IN);

if ($inflag eq '-r') { # Build list from the command line
	if (!($inlist =~ m/\-/)) { die "Need a hyphen in the number list.\n"; }

	@startStop = split '-',$inlist;

	if ($inlist =~ m/^\-/) { $start = 0; } else { $start = $startStop[0]; }
	if ($inlist =~ m/\-$/) { $finish = $header - 1; } else {$finish = $startStop[1];}

	print STDOUT "Collecting indices between ", $start, " and ", $finish, "\n";

	foreach $i ($start..$finish) { # Build the list of indices to capture
		push (@inputs, $i);
	}
} elsif ($inflag eq '-f') { # Build the list from a file containing index numbers

	open (INLIST, "$inlist") or die "Can't open listfile $inlist\n";
	@inputs = <INLIST>;
	close (INLIST);

} else {
	die "Invalid command line flag. Legal options are -l and -f.\n";
}

foreach $inp (@inputs) {
	$pos = $inp + $numNonInd; # To correct for presence of PosNeg and TrainTest
	for $j (0..$last_seq) {
		push (@{$new_ind[$j]}, $allsplit[$j][$pos]);
	}
}

$header = scalar (@inputs);

open (OUT, ">$outfile");

print OUT "$header\n";
print OUT $sizes;

# Average the scores for each input
for $k (0..$last_seq) {
	for $kk (0..$numNonInd - 1) { # PosNeg and TrainTest and whatever else
		print OUT "$allsplit[$k][$kk],";
	}
	foreach $bit (@{$new_ind[$k]}) {
		print OUT "$bit,";
	}
	print OUT "\n";
}

close (OUT);
