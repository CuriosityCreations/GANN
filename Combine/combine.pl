use strict;

$| = 1;

sub MAX (@) {
	my $max = $_[0];
	for my $i (1..$#_) {
		if ($_[$i] > $max) { $max = $_[$i]; }
	}
	return $max;
}

sub MIN (@) {
	my $min = $_[0];
	for my $i (1..$#_) {
		if ($_[$i] < $min) { $min = $_[$i]; }
	}
	return $min;
}

sub MEAN (@) {
	my $mean = $_[0];
	for my $i (1..$#_) {
		$mean += $_[$i];
	}
	return ($mean / (scalar @_));
}

# Generate a random column with m zeroes and n ones
sub RANDCOL ($$) {
	my $colSize = $_[0] + $_[1];
	my @generateCol = ();
	for my $i (0..$colSize - 1) { $generateCol[$i] = 0; }
	my $numOnes = 0;
	while ($numOnes < $_[1]) {
		my $checkMe = int rand $colSize;
		if ($generateCol[$checkMe] != 1) {
			$generateCol[$checkMe] = 1;
			++$numOnes;
		}
	}
	return @generateCol;
}

# Process the command line flags
# Command line options:
#	-n = no Z scores
#	-g = Global values
#	-h xxx = restrict to files containing xxx
#	-w n1 n2 = only consider windows from n1 to n2
#	-s = standardize columns
#	-e n = number of experimental sets (with random train / test partitions)
#	-f n = number of negative control sets

my $usage = "Usage: perl combine2004.pl outfile [-g] [-noZ] [-w lowWin hiWin] [-e numExp] [-f numNeg] [-h handle]\n";
my $Zscores = 1;
my $global = 0;
my $handle = "";
my $restrictWin = 0;
my $maxWin;
my $minWin;
my $standardize = 0;
my $numExp = 1;
my $numNeg = 0;
my $dataSetSize;

if (scalar(@ARGV) < 1) { die "$usage\n"; }
if ($ARGV[0] =~ /^\-/) { die "$usage\n"; }
my $outfile = $ARGV[0];
for my $i (1..$#ARGV) {
	if ($ARGV[$i] =~ /^\-/) {
		my $thisChar = substr($ARGV[$i],1,1);
		if ($thisChar eq 'n') {	$Zscores = 0; }
		elsif ($thisChar eq 'g') { $global = 1; }
		elsif ($thisChar eq 'h') { $handle = $ARGV[++$i]; }
		elsif ($thisChar eq 'w') { $minWin = $ARGV[++$i]; $maxWin = $ARGV[++$i]; $restrictWin = 1; }
		elsif ($thisChar eq 'e') { $numExp = $ARGV[++$i]; }
		elsif ($thisChar eq 'f') { $numNeg = $ARGV[++$i]; }
	}
}
if ($numExp + $numNeg == 0) { die "Nothing to do!\n"; }

# Which files need to be read?
my @fileList = ();
opendir(DIR,'.');
if ($handle ne "") { # Get all of a specific file type
	while (my $thing = readdir(DIR)) {
		if ($thing =~ m/$handle/) {
			push @fileList, $thing;
		}
	}
} else {
	print "Which files should be combined?\n\n";
	while (my $thing = readdir(DIR)) {
		my $answered = 0;
		while ($answered == 0) {
			print "$thing? ";
			my $answer = <STDIN>;
			my $firstBit = substr(uc $answer, 0, 1);
			if ($firstBit eq 'Y') {	push @fileList, $thing; $answered = 1; }
			elsif ($firstBit eq 'N') { $answered = 1; }
			else { print "Please type Y or N.\n"; }
		}
	}
}
closedir(DIR);

# Create a 'master list' of entities
my $standard = $fileList[0]; # The first file sets the standard
my %masterList = ();
my $posSize = 0;
my $negSize = 0;
open (IN, "$standard") or die "Aaargh!\n";
my $toss = <IN>; # Not checking headings at the moment
while (<IN>) {
	chomp;
	if (/\,/) {
		s/\,$//g;
		my @currentLine = split /\,/,$_;
		if (!(defined($masterList{$currentLine[1]}))) {
			$masterList{$currentLine[1]} = $currentLine[0];
			if ($currentLine[0] == 0) { ++$negSize; } else { ++$posSize; }
		} else { die "Error creating master list. $currentLine[1] found twice in $standard.\n"; } 
	}
}

$dataSetSize = scalar (keys %masterList);
# Read each file, check for consistency
my %indices = ();
my @successFiles = ();
foreach my $thisFile (@fileList) {
	my %assoc = ();
	my %foundNames = ();
	my %colsToUse = (); # Only used if window IDs are restricted
	open (IN, "$thisFile") or die "Aaargh!\n";
	print "\nOpening $thisFile...\t";
	my $header = <IN>;
	chomp $header;
	$header =~ s/\,$//g;
	my @headerBits = split /\,/,$header;
	my $success = 0;
	READFILE: if (($headerBits[0] eq "PosNeg") && ($headerBits[1] eq "SeqID")) {
		for my $thisHeading (2..$#headerBits) { 
			if ($restrictWin == 1) {
				if ($headerBits[$thisHeading] =~ /n(\d+)$/) {
					if (($1 >= $minWin) && ($1 <= $maxWin)) {
						$assoc{$thisHeading} = $headerBits[$thisHeading];
					}
				}
			} else { $assoc{$thisHeading} = $headerBits[$thisHeading]; }
		}
		while (<IN>) {
			chomp;
			s/\,$//g;
			my @thisLine = split /\,/,$_;
			if ($#thisLine != $#headerBits) {
				print "$thisLine[1] does not have the right number of indices.\n";
				last READFILE;
			} else {
				for my $addToList (2..$#thisLine) {
					if (defined($assoc{$addToList})) {
						$indices{$assoc{$addToList}}{$thisLine[1]} = $thisLine[$addToList];
# print "$addToList\t$assoc{$addToList}\t$thisLine[1]\n";
					}
				}
			}
			if (!(defined($foundNames{$thisLine[1]}))) {
				$foundNames{$thisLine[1]} = 1;
			} else { 
				print "$thisLine[1] found multiple times in file $thisFile.\n";
				last READFILE;
			}
		}
		# Check foundNames for consistency with the master list
		foreach my $inMaster (keys %masterList) {
			if (!(defined($foundNames{$inMaster}))) {
				print "$inMaster entity missing from file $thisFile.\n";
				last READFILE;
			}
		}
		$success = 1;
	} 
	if ($success == 0) { print "$thisFile does not conform to the required format, skipping...\n"; }
	else { push @successFiles, $thisFile; }
	close (IN);
}
# Create the output set of indices, calculating Z-scores, standardizing and globals 
my %masterOut = ();
my @varList = sort keys %indices;
if ($global == 1) { # Globals?
	# my @globalize = (); # List of variable ranges to 'globalize'
	my %globalHandles = ();
	for my $i (0..$#varList) { 
		if ($varList[$i] =~ m/(.*)n(\d+)$/) { # Then it has a window number
			my $thisHandle = $1;
			push (@{$globalHandles{$thisHandle}}, $varList[$i]);
		}
	}
	# Generate the global vars for each set of windows
	foreach my $globalSet (keys %globalHandles) {
		foreach my $line (keys %{$indices{$varList[0]}}) {
			my @currentSet = ();
			foreach my $getVar (@{$globalHandles{$globalSet}}) {
				push (@currentSet,$indices{$getVar}{$line});
			}
			$indices{$globalSet . 'MAX'}{$line} = MAX(@currentSet);
			$indices{$globalSet . 'MIN'}{$line} = MIN(@currentSet);
			$indices{$globalSet . 'MEAN'}{$line} = MEAN(@currentSet);
		}
	}
} # End globals
if ($Zscores == 1) { # Then standardize
	foreach my $var (keys %indices) {
		my ($sum, $sumSq) = 0;
		foreach my $getName (keys %{$indices{$var}}) {
			my $getVal = $indices{$var}{$getName};
			$sum += $getVal;
			$sumSq += $getVal * $getVal;
		}
		my $mean = $sum / $dataSetSize;
		my $stdev;
		if ($dataSetSize <= 1) { $stdev = 0; }
		else { $stdev = sqrt (($sumSq - ($sum * $sum / $dataSetSize)) / ($dataSetSize - 1)); }
		foreach my $replaceVal (keys %{$indices{$var}}) {
print $indices{$var}{$replaceVal};
			$indices{$var}{$replaceVal} -= $mean;
			if ($stdev != 0) { $indices{$var}{$replaceVal} /= $stdev; }
print "\t$indices{$var}{$replaceVal}\n";
		}
	}
}

my $posTrainSize = int ($posSize * 0.66);
my $posTestSize = $posSize - $posTrainSize;
my $negTrainSize = int ($negSize * 0.66);
my $negTestSize = $negSize - $negTrainSize;

# Generate the requested number of experimental output files
for my $i (1..$numExp) {
	open (OUT,">${outfile}E$i.txt");
	print OUT scalar keys %indices,"\n";
	print OUT "$negTrainSize,$posTrainSize,$negTestSize,$posTestSize\n";
	my @trainTest;
	@{$trainTest[1]} = RANDCOL ($posTrainSize,$posTestSize);
	@{$trainTest[0]} = RANDCOL ($negTrainSize,$negTestSize);
	my @outSet = (); my $outRow = 0;
	print OUT "PosNeg,TrainTest,SeqID";
	foreach my $varName (sort { $a <=> $b } keys %indices) {
		print OUT ",$varName";
	}
	print OUT "\n";
	foreach my $line (keys %masterList) {
		my $PosNeg = $masterList{$line};
		print OUT "$PosNeg,$trainTest[$PosNeg][$outSet[$PosNeg]++],$line";
		foreach my $varName (sort { $a <=> $b } keys %indices) {
			print OUT ",$indices{$varName}{$line}";
		}
		print OUT "\n";
	}
	close (OUT);
}
# Generate the requested number of negative control output files
for my $i (1..$numNeg) {
	open (OUT,">${outfile}N$i.txt");
	print OUT scalar keys %indices,"\n";
	print OUT "$negTrainSize,$posTrainSize,$negTestSize,$posTestSize\n";
	my @newPosNeg = RANDCOL ($negSize,$posSize);
	my @trainTest;
	@{$trainTest[1]} = RANDCOL ($posTrainSize,$posTestSize);
	@{$trainTest[0]} = RANDCOL ($negTrainSize,$negTestSize);
	my @outSet = (); my $outRow = 0;
	print OUT "PosNeg,TrainTest,SeqID";
	foreach my $varName (sort { $a <=> $b } keys %indices) {
		print OUT ",$varName";
	}
	print OUT "\n";
	foreach my $line (keys %masterList) {
		my $PosNeg = $newPosNeg[$outRow++];
		print OUT "$PosNeg,$trainTest[$PosNeg][$outSet[$PosNeg]++],$line";
		foreach my $varName (sort { $a <=> $b } keys %indices) {
			print OUT ",$indices{$varName}{$line}";
		}
		print OUT "\n";
	}

	close (OUT);
}