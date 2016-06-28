// indices.cpp converts DNA sequence to a series of numbers using user-specified mapping rules.
// Copyright (C) 2001-2002 Robert G. Beiko

// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// indices.cpp takes a sequence file and a file containing a series of maps as input,
// then converts the set of sequences to a series of indices, representing
// the conversions specified in the map files, nucleotide interval values, and
// nucleotide frequencies.

#include"indices.h"

int main(int args, char * argv[]) {

	string seqfile, outfile, mapfile, oligofile, MGAfile, win_size_txt = "";
	int win_size, overlap, nmer_size, nucint_size, num_reps = 0;
	bool seqOpen, mapOpen, oligoOpen, MGAopen = false; // open status for seqfile, outfile, mapfile, oligofile
	bool useType[5] = { false }; // Which of the four types of analyses to use

	map<char,int> nucDegen; // How many different nucleotides are represented by each degen letter?
	nucDegen['A'] = 1; nucDegen['C'] = 1; nucDegen['G'] = 1; nucDegen['T'] = 1;
	nucDegen['R'] = 2; nucDegen['Y'] = 2; nucDegen['M'] = 2; nucDegen['K'] = 2; nucDegen['S'] = 2; nucDegen['W'] = 2;
	nucDegen['B'] = 3; nucDegen['D'] = 3; nucDegen['H'] = 3; nucDegen['V'] = 3;
	nucDegen['N'] = 4; nucDegen['X'] = 4;

	// Processing inputs
	if (args < 6) {
		cout << USAGE;
		exit(-1);
	} else {
		// Record the essential inputs
		seqfile = argv[1];
		outfile = argv[2];
		win_size = atoi(argv[3]);
		win_size_txt = argv[3];
		overlap = atoi(argv[4]);
		// Record the optional inputs
		int inpcount = 5;
		while (inpcount < args) {
			string flag = argv[inpcount++];
			if (flag[0] == '-') {
				switch (flag[1]) {
					case 'm':
						mapfile = argv[inpcount++];
						useType[0] = true;
						break;
					case 'o':
						oligofile = argv[inpcount++];
						useType[1] = true;
						break;
					case 'n':
						nmer_size = atoi(argv[inpcount++]);
						useType[2] = true;
						break;
					case 'i':
						nucint_size = atoi(argv[inpcount++]);
						useType[3] = true;
						break;
					case 'g': // MGA representation
						MGAfile = argv[inpcount++];
						useType[4] = true;
						break;
					case 'r':
						num_reps = atoi(argv[inpcount++]);
						break;
					default:
						cout << "Invalid option " << flag << "\n\n" << "Usage: " << USAGE << "\n";
						exit(-1);
				}
			} else { // Flag does not start with '-'
				cout << "Invalid option " << flag << "\n\n" << "Usage: " << USAGE << "\n";
				exit(-1);
			}	
		}
	}		

	// Open the sequence file
	cout << "Reading sequence...";
	vector<string> seqVec;
	vector<string> IDVec;
	int numSeq = 0;
	seqOpen = readSeqFile(seqfile, seqVec, IDVec, numSeq);
	if (!seqOpen) {
		cout << "Couldn't open sequence file " << seqfile << endl;
		exit(-1);
	}
	cout << "Done.\n";

	// Chop the sequences into windows, create Z-score replicates if appropriate
	cout << "Chopping sequence into windows..." << flush;
	vector<vector<vector<string> > > seqFrags(numSeq);
	chopSequence(seqVec, win_size, overlap, num_reps, seqFrags);
	cout << "Done.\n";

	// Print the chopped sequences to a file
	seqOutput(seqFrags, IDVec, outfile);

	// Open and interpret the map file
	vector<map<string, double> > ruleVec;
	vector<string> ruleNames;
	if (useType[0]) { // Mapfile specified
		cout << "Reading map file..." << flush;
		mapOpen = readMapFile(mapfile, ruleVec, ruleNames, nucDegen);
		if (!mapOpen) {	
			cout << "Couldn't open map file " << mapfile << endl;
			exit(-1);
		}
		cout << "Done.\n";
	}

	// Open and interpret the oligo file
	vector<string> oligoVec;
	if (useType[1]) { // Oligofile specified
		oligoOpen = readOligoFile(oligofile, oligoVec);
		if (!oligoOpen) {
			cout << "Couldn't open oligo file " << oligofile << endl;
			exit(-1);
		}
	}
	
	// Open and interpret the MGA file
	vector<Chromosome*> MGAvec;
	if (useType[4]) { // MGAfile specified
		ifstream MGAin(MGAfile.c_str());
		if (!MGAin) { cout << "Couldn't open MGA file " << MGAfile << endl; exit (-1); }
		readMGAfile(MGAin, MGAvec, win_size, ruleVec, ruleNames);
	}

	if (useType[0]) { mapConvert(seqFrags, ruleVec, ruleNames, IDVec, outfile + "win" + win_size_txt); }
	if (useType[1]) { 
		cout << "Counting oligonucleotides..." << flush;
		oligoConvert(seqFrags, oligoVec, IDVec, outfile + "win" + win_size_txt, "oligos"); 
		cout << "Done.\n";
	}
		
	if (useType[2]) { nmerConvert(seqFrags, nmer_size, IDVec, outfile + "win" + win_size_txt); }
	if (useType[3]) { nucintConvert(seqFrags, nucint_size, IDVec, outfile + "win" + win_size_txt); }
	
	if (useType[4]) { MGAconvert(seqFrags, MGAvec, ruleVec, ruleNames, IDVec, outfile + "win" + win_size_txt + "MGA"); }

	for (int delChr = 0; delChr < MGAvec.size(); ++delChr) {
		delete MGAvec[delChr];
		MGAvec[delChr] = 0;
	}
	cout << "Goodbye!\n";

	return  0;
}

// Read the sequence file, split each entry into (category + name), (sequence)
bool readSeqFile(const string& seqfile, vector<string>& seqVec, vector<string>& IDVec, int& numSeq) {

	string currLine;
	int delim;

	ifstream inSeq(seqfile.c_str());
	if (!inSeq) return false;

	while (!inSeq.eof()) {
		getline(inSeq, currLine);
		delim = currLine.rfind(","); // Find the last comma, which should separate ID from seq
		if (delim != string::npos) {
			IDVec.push_back(currLine.substr(0, delim)); // Store ID
			seqVec.push_back(currLine.substr(delim + 1)); // Store sequence
			++numSeq;
		}
	}

	inSeq.close();

	return true;
}

// Read the map file, build each map and store the name of each set of rules
// The SWITCH statement in this function is ugly but tolerant of variation in the structure of the input file
bool readMapFile(const string& mapfile, vector<map<string, double> >& ruleVec, vector<string>& ruleNames, map<char,int>& nucDegen) {

	string currLine;
	// int thePos;	

	int readPos = 0; // Define the data expected
	int pos = 0, start = 0, finish = 0; // Position being searched
	int len = 0; // Length of the current string
	int power, num_cases = 0;
	int num_counter = 0, mapCount = 0;
	int casesCovered = 0;
	bool searchForName; // Search for name or conversion value?
	pair<string, double> toAdd; 
	string nmer;

	map<string, double> inProgress; // The current map being read, which must be (DEEP!) copied to the vector when complete

	ifstream inMap(mapfile.c_str());
	if (!inMap) return false;
	
	bool finished = false;
	bool getNextLine = true;

	while (!finished) {
		if (getNextLine) { 
			getline(inMap, currLine); 
			if (readPos > 0) getNextLine = false;
		}

		if (inMap.eof() && readPos == 0) { finished = true; } else { finished = false; } // Only stop if at EOF
		switch(readPos) {
			case 0: // Outside the map specification
				if (currLine.find("{") != string::npos) { // If { is found
					currLine.erase(0, currLine.find("{") + 1); // Toss everything before the {
					getNextLine = false; // Look for stuff on this line
					finished = false;
					++readPos;
				}
				break;
			case 1: // Inside specification, looking for NAME
				pos = 0;				
				len = currLine.length();
				while (!getNextLine && (readPos == 1)) {
					if (pos == len) { 
						getNextLine = true; 
					} else {
						if (isgraph(currLine[pos])) { // Found a legal name character (letter, number, punctuation)
							start = pos;
							while (isgraph(currLine[pos])) { ++pos; }
							finish = pos;
							ruleNames.push_back(currLine.substr(start, finish - 1)); // Add the name to the list
							currLine.erase(0, finish);
							readPos = 2;
						} else { 
							++pos; 
						}
					}
				}
				break;
			case 2: // Inside specification, looking for exponent 4^n (n is a SINGLE DIGIT)
				pos = 0;
				len = currLine.length();
				while (!getNextLine && (readPos == 2)) {
					if (pos == len) {
						getNextLine = true;
					} else {
						if (isdigit(currLine[pos])) { // Found a digit
							power = atoi((currLine.substr(pos, 1)).c_str()); // Record the exponent
							num_cases = (int) pow(4.0,power); // Hopefully 4 or 16 or 64 and not 4^9
							currLine.erase(0,pos + 1);
							++readPos;
						}
						++pos;
					}
				}
				break;
			case 3: // Inside specification, looking for { to open map definition
				if (currLine.find("{") != string::npos) { // If { is found
					currLine.erase(0, currLine.find("{") + 1); // Toss everything before the {
					searchForName = true;
					++readPos;
				} else {
					getNextLine = true;
				}
				break;
			case 4: // Inside the definition, looking for members to add to the map. Members
					// must be specified in the form (degen_seq rule) with any (ignored) separators
				pos = 0;
				len = currLine.length();
				while (!getNextLine && (readPos == 4)) {
					if (pos == len) {
						getNextLine = true;
					} else {
						char toCheck = currLine[pos];
						if (toCheck == '}') { // Close the map definition, fail if incorrect # of examples were counted
							if (num_counter != num_cases) {
								cout << "Incorrect number of cases (" << num_counter << ") specified in map definition " << mapCount + 1 << endl;
								exit(-1);
							}
							currLine.erase(0, currLine.find("}") + 1);
							++readPos; 
						} else {
							if (searchForName) { // Looking for degen nucleotides
								if (isalpha(toCheck)) { // Found the start of the degen word
									casesCovered = 1; // Number of specific n-mers described by degen example
									nmer = currLine.substr(pos, power);
									toAdd.first = nmer;
									for (int countD = 0; countD < nmer.length(); ++countD) { // How many cases are represented by this nmer?
										casesCovered *= nucDegen[nmer[countD]];
									}
									num_counter += casesCovered;
									pos += power - 1;
									searchForName = false; // Start looking for a floating-point value
								} else { ++pos; }					
							} else { // Looking for floating-point value
								if (isdigit(toCheck) || toCheck == '-') { 
									start = pos;
									char here = currLine[++pos];
									while (isdigit(here) || here == '.') { 
										here = currLine[++pos]; 
									}
									toAdd.second = atof((currLine.substr(start, pos - start)).c_str());
									string buildNmer = "";
	 								buildMap(toAdd, 0, power, buildNmer, inProgress); // Recursive function to add all non-degenerate words to map
	 								searchForName = true; // Switch to looking for the next word 
								} else { ++pos; }
							} 
						}
					}
				}
				break;
			case 5: // Looking for a } to terminate the map specification
				getNextLine = true;
				if (currLine.find("}") != string::npos) { // If } is found
					currLine.erase(0, currLine.find("}") + 1); // Toss everything before the }
					ruleVec.push_back(inProgress);
					inProgress.erase(inProgress.begin(),inProgress.end()); // Reset inProgress
					++mapCount;
					readPos = 0; // Next!
					num_counter = 0;
				}
				break;
			} // Switch statement 
		}	
	inMap.close();
	return true;

}

// Take a degenerate n-mer, and add all of its non-degenerate members to a map using recursion
void buildMap(const pair<string, double>& toAdd, const int currentPos, const int finalPos, string buildNmer, map<string, double>& inProgress) {
	
	if (currentPos == finalPos) {
		map<string, double>::iterator IP = inProgress.find(buildNmer);
		if (IP == inProgress.end()) { // If not found in map
			inProgress[buildNmer] = toAdd.second;
		} else {
			cout << "Duplicate word " << buildNmer << " in map." << endl;
 		    exit(-1);
 		}
		// pair<string, double> addThis(buildNmer, toAdd.second);
		// pair<map<string, double>::iterator, bool> success = inProgress.insert(addThis);
		// if (!success.second) { // If failed, then a duplicate key was specified
 		// 	cout << "Duplicate word " << addThis.first << " in map." << endl;
 		// 	exit(-1);
 		// }
 	} else {
		int bases[4] = { 0, 0, 0, 0 }; // A, C, G, T
		char theChar = (char) toupper((toAdd.first)[currentPos]);
		bool trueChar = false;
		if (theChar == 'N' || theChar == 'X') { // Any character. Don't add 4 entries to the map
			buildMap(toAdd, currentPos + 1, finalPos, buildNmer + "N", inProgress);
		} else {
			switch (theChar) {
				case 'A': case 'R': case 'W': case 'M': case 'D': case 'H': case 'V': // Adenine
					bases[0] = 1;
					trueChar = true;
			}
			switch (theChar) {
				case 'C': case 'Y': case 'S': case 'M': case 'B': case 'H': case 'V': // Cytosine
					bases[1] = 1;
					trueChar = true;
			}
			switch (theChar) {
				case 'G': case 'R': case 'S': case 'K': case 'B': case 'D': case 'V': // Guanine
					bases[2] = 1;
					trueChar = true;
			}
			switch (theChar) {
				case 'T': case 'Y': case 'W': case 'K': case 'B': case 'D': case 'H': // Thymine
					bases[3] = 1;
					trueChar = true;
			}
			if (!trueChar) {
					cout << "Invalid character " << theChar << " in word " << toAdd.first << endl;
					exit (-1);
			}
			if (bases[0]) { buildMap(toAdd, currentPos + 1, finalPos, buildNmer + "A", inProgress); }
			if (bases[1]) { buildMap(toAdd, currentPos + 1, finalPos, buildNmer + "C", inProgress); }
			if (bases[2]) { buildMap(toAdd, currentPos + 1, finalPos, buildNmer + "G", inProgress); }
			if (bases[3]) { buildMap(toAdd, currentPos + 1, finalPos, buildNmer + "T", inProgress); }
		}
	}
}

// Read the oligo file
bool readOligoFile(const string& oligofile, vector<string>& oligoVec) {

	ifstream inOligo(oligofile.c_str());
	if (!inOligo) { return false; }
	
	string currLine;
	int len;
	int pos;
	int start, finish = 0;

	while (!inOligo.eof()) {
		getline(inOligo, currLine);
		len = currLine.length();
		pos = 0;
		while (pos < len) {
			if (isalpha(currLine[pos])) { // The beginning of an n-mer
				start = pos;
				while (isalpha(currLine[pos])) { 
					++pos; 
				} // Increment pos until no more alpha
				finish = pos; 
				oligoVec.push_back(currLine.substr(start, finish));
			} else { ++pos; } // Nothing found: keep looking
		}
	}
	inOligo.close();
	return true;
}

// Function to cut sequence into fragments and randomly shuffled replicates
void chopSequence(const vector<string>& seqVec, const int win_size, const int overlap, const int numReps,
	vector<vector<vector<string> > >& seqFrags) {
	
	int seqLen = (seqVec[0]).length();
	int firstPos = seqLen - win_size - 2 - 1; // Explicit: 2 for flanking bases, 1 for array subscripting
	int seqNum = seqVec.size();
	int decrement = win_size - overlap;
	int currentSeq, startPos, totalWin = 0;

	int winCount = firstPos / decrement + 1;

	for (currentSeq = 0; currentSeq < seqNum; ++currentSeq) {
		(seqFrags[currentSeq]).resize(winCount);
		int winNum = 0;
		for (startPos = firstPos + 1; startPos >= 0; startPos -= decrement) {
			(seqFrags[currentSeq][winNum]).resize(1 + numReps);
			seqFrags[currentSeq][winNum][0] = (seqVec[currentSeq]).substr(startPos, win_size + 1);
			++winNum;
			// FOR A TOTAL OF WIN_SIZE + 2 NUCLEOTIDES
		}
	}
	if (numReps > 0) {
		totalWin = (seqFrags[0]).size(); // How many windows?
		for (currentSeq = 0; currentSeq < seqNum; ++currentSeq) {
			for (startPos = 0; startPos < totalWin; ++startPos) {
				for (int rep = 1; rep <= numReps; ++rep) { // Randomize the sequence numReps times
					seqFrags[currentSeq][startPos][rep] = seqShuffle(seqFrags[currentSeq][startPos][0], win_size);
				}
			}
		}
	}
}

// Print the sequence fragments to a file
void seqOutput (const vector<vector<vector<string> > >& theFragments, const vector<string>& IDVec, const string& outfile) {

	string theOut = outfile + ".seqs";

	ofstream seqOut(theOut.c_str());

	int numSeq = theFragments.size();
	int numFrags = (theFragments[0]).size();
	
	for (int titleOut = 0; titleOut < numFrags; ++titleOut) {
		seqOut << titleOut << ",";
	}
	seqOut << endl;

	for (int i = 0; i < numSeq; ++i) {
		seqOut << IDVec[i] << ",";
		for (int j = 0; j < numFrags; ++j) {
			seqOut << theFragments[i][j][0] << ",";
		}
		seqOut << endl;
	}
	seqOut.close();
}

// Read the MGA file and generate a bunch of Chromosomes
void readMGAfile(ifstream& inFile, vector<Chromosome*>& MGAvec, const int winSize, 
	vector<map<string, double> >& ruleVec, vector<string>& ruleNames) {

	string chrLine = "";
	getline(inFile, chrLine);
	while (!(inFile.eof())) {
		Chromosome* nextChr = new Chromosome(chrLine, winSize, ruleVec, ruleNames);
		if (nextChr->maxLen() < winSize) { 
			MGAvec.push_back(nextChr);
		} else {
			delete nextChr;
			nextChr = 0;
		}
		getline(inFile, chrLine);
	}
}

// Shuffle a sequence
string seqShuffle (const string& toShuf, const int win_size) {

	string product = toShuf;

	for (int i = 0; i < win_size; ++i) {
		int shufWith = (int) rand() % win_size;
		char temp = product[i];
		product[i] = product[shufWith];
		product[shufWith] = temp;
	}
	
	return product;
}

// Function to convert sequence to numbers using a complete map representation
void mapConvert(const vector<vector<vector<string> > >& seqFrags, const vector<map<string, double> >& ruleVec, 
	const vector<string>& ruleNames, const vector<string>& IDVec, const string outFile) {
	
	int numRules = ruleVec.size();
	int numSeq = seqFrags.size();
	int numWin = (seqFrags[0]).size();
	int numRep = (seqFrags[0][0]).size() - 1; // -1 since the first element is the real sequence fragment
	int numVar = (numRep > 0) ? numWin * 2: numWin; // * 2 if Z-scores are to be generated
	int mult = (numRep > 0) ? 2 : 1;
	int winSize = (seqFrags[0][0][0]).length();

	ostringstream outStr;

	cout << "Converting map: " << flush;

	for (int i = 0; i < numRules; ++i) {

		int numOligos = (ruleVec[i]).size();

		vector<vector<double> > genIndices(numSeq);
		for (int sCount = 0; sCount < numSeq; ++sCount) {
			(genIndices[sCount]).resize(numVar);
		}
		map<string, double> thisMap = ruleVec[i]; // To reduce confusion later on
		int ruleLength = ((thisMap.begin())->first).length();
		int lenMinusOne = ruleLength - 1;
		int winFinish = winSize - lenMinusOne;

		string thisName = ruleNames[i];
		vector<string> varNames;
		varNames.push_back("PosNeg"); varNames.push_back("SeqID");
		for (int winNames = 0; winNames < numWin; ++winNames) { // Create the list of window names
			outStr << thisName << "ws" << winSize - 1 << "n" << winNames << '\0';
			varNames.push_back(outStr.str()); // Variable value
			if (numRep > 0) {
				outStr.seekp((int) outStr.tellp() - 1);
				outStr << "Z" << '\0';
				varNames.push_back(outStr.str()); // Z-score
			}
			outStr.seekp(0);
		}

		cout << thisName << "..." << flush;

		for (int j = 0; j < numSeq; ++j) { // For each sequence
			for (int k = 0; k < numWin; ++k) { // For each window
				double runningTotal = 0.0;
				string toCalc = seqFrags[j][k][0]; // Real sequence
				for (int pos = 0; pos < winFinish; ++pos) {
					runningTotal += thisMap[toCalc.substr(pos, ruleLength)];
				}
				genIndices[j][k * mult] = runningTotal;
				if (numRep > 0) {
					double repSum = 0.0;
					double repSumSq = 0.0;
					for (int repCount = 1; repCount <= numRep; ++repCount) {
						string repToCalc = seqFrags[j][k][repCount];
						double thisRep = 0;
						for (int repPos = 0; repPos < winFinish; ++repPos) {
							double addVal = thisMap[repToCalc.substr(repPos, ruleLength)];
							thisRep += addVal;
						}
						repSum += thisRep;
						repSumSq += thisRep * thisRep;
					}
					double Zsc = Zscore(runningTotal, repSum, repSumSq, numRep);
					genIndices[j][k * mult + 1] = Zsc;
				} // Done with reps
			} // Done with this window of sequence
		} // Done with this sequence
		indOutput(varNames, genIndices, IDVec, outFile + thisName);
		varNames.resize(0);
		genIndices.resize(0);
		
	} // Done with this map representation
	cout << "Done.\n";

}

// Calculate oligo counts within sequences. Every oligo count will be output to the same file
void oligoConvert(const vector<vector<vector<string> > >& seqFrags, const vector<string>& oligoVec, 
	const vector<string>& IDVec, const string outFile, const string outID) {

	int numOligos = oligoVec.size();
	int numSeq = seqFrags.size();
	int numWin = (seqFrags[0]).size();
	int numRep = (seqFrags[0][0]).size() - 1; // -1 since the first element is the real sequence fragment
	int numVar = (numRep > 0) ? numWin * 2: numWin; // * 2 if Z-scores are to be generated
	int mult = (numRep > 0) ? 2 : 1;
	int winSize = (seqFrags[0][0][0]).length(); // = the requested win_size + 1

	string varName = "";
	ostringstream outStr;

	// Build the list of variables to be output
	vector<string> OligoNames;
	OligoNames.push_back("PosNeg"); OligoNames.push_back("SeqID");
	for (int oName = 0; oName < numOligos; ++oName) {
		string thisName = oligoVec[oName];
		for (int winNames = 0; winNames < numWin; ++winNames) { // Create the list of window names
			outStr << thisName << "ws" << winSize - 1 << "n" << winNames << '\0';
			OligoNames.push_back(outStr.str()); // Variable value
			if (numRep > 0) {
				outStr.seekp((int) outStr.tellp() - 1);
				outStr << "Z" << '\0';
				OligoNames.push_back(outStr.str()); // Z-score
			}
			outStr.seekp(0);
		}
	}
	
	vector<vector<double> > genIndices(numSeq); // Will hold all of the oligo representations
	for (int sCount = 0; sCount < numSeq; ++sCount) {
		(genIndices[sCount]).resize(numVar * numOligos);
	}

	for (int thisOligo = 0; thisOligo < numOligos; ++thisOligo) {

		string oligo = oligoVec[thisOligo];
		// Build a list of non-degenerate oligos
		map<string, double> currentOligoRep;
		pair<string, double> thePair; thePair.first = oligo; thePair.second = 1.0;

		int oLength = oligo.length();
		string inProgress = "";
		buildMap(thePair, 0, oLength, inProgress, currentOligoRep);
		int oligoBase = numVar * thisOligo;
		for (int j = 0; j < numSeq; ++j) { // For each sequence
			for (int k = 0; k < numWin; ++k) { // For each window
				double total = 0.0;
				string toCalc = seqFrags[j][k][0]; // Real sequence
				// Iterate over every non-degen word, count them all (INEFFICIENT)
				for (CI whichSearch = currentOligoRep.begin(); whichSearch != currentOligoRep.end(); ++whichSearch) {
					total += enumerateHits(whichSearch, toCalc);
				} 
				genIndices[j][oligoBase + k * mult] = total;

				// Deal with the replicates and generate a Z-score
				if (numRep > 0) {
					double repSum = 0.0;
					double repSumSq = 0.0;
					double thisRep = 0.0;
					for (int repCount = 1; repCount <= numRep; ++repCount) {
						thisRep = 0.0;
						string repToCalc = seqFrags[j][k][repCount];
						for (CI repSearch = currentOligoRep.begin(); repSearch != currentOligoRep.end(); ++repSearch) {
							thisRep += enumerateHits(repSearch, repToCalc);
						}
						repSum += thisRep;
						repSumSq += thisRep * thisRep;
					} // Done with reps
					double Zsc = Zscore(total, repSum, repSumSq, numRep);
					genIndices[j][oligoBase + k * mult + 1] = Zsc;
				} // Again, repeat for all degenerate nmers
			} // Done with this window of sequence
		} // Done with this sequence
	}
	indOutput(OligoNames, genIndices, IDVec, outFile + outID);
} 

// Count the number of hits for a given n-mer against a target sequence. 'N' is treated as a 'don't care' character
double enumerateHits (const CI& whichSearch, const string& toCalc) {

	// Convert the search string into a representation of specific oligos and 'N's
	string theString = whichSearch->first;
	double theVal = whichSearch->second;
	
	double runningTotal = 0.0;

	int startPt = 0;
	int findPos = degen_find(theString, toCalc.substr(startPt));
	while (findPos != string::npos) {
		findPos += startPt;
		runningTotal += theVal;
		startPt = findPos + 1;
		findPos = degen_find(theString, toCalc.substr(startPt)); // Must check for npos, THEN add start_pt
	}
	return runningTotal;
}

// Search for degenerate patterns in a string. Eventually it would be good to implement a real text-searching
// algorithm (i.e. BNDM), but this would probably be best implemented waaaay back in BuildMap() to store bit masks
// in the map
int degen_find(const string& theString, const string& theTarget) {

	bool found = false;
	bool goOn = true;

	int targLen = theTarget.length();
	int searchLen = theString.length();
	int stop = targLen - searchLen + 1;

	int p = 0;
	int i = 0;

	while (!found) {
		if (i == stop) return -1; // String not found
		while (goOn) {
			if (p == searchLen) return i; // String found
			if (theString[p] == theTarget[i + p] || theString[p] == 'N') { ++p; } // OK: Increment search index
			else { goOn = false; } // No good: Go to next i
		}
		p = 0;
		++i;
		goOn = true;
	}
}

// Count the frequency of each nmer up to size (nmer_size) in the sequence. Relies heavily on oligoConvert
// to do the dirty work
void nmerConvert(const vector<vector<vector<string> > >& seqFrags, const int nmer_size, 
	const vector<string>& IDVec, const string outFile) {

	ostringstream oNmerName;

	cout << "Counting n-mers, n = " << flush;

	// For each size of nmers from 1 to nmer_size
	for (int i = 1; i <= nmer_size; ++i) {

		cout << i << "..." << flush;

		int num_vars = (int) pow((double) i, 4);

		// Create a vector of all non-degenerate nmers
		vector<string> nmers;
		fillVector(i, 0, "", nmers);
	
		oNmerName << "nmer" << i << '\0';

		// Call the countOligos function to count all of the nmers
		oligoConvert(seqFrags, nmers, IDVec, outFile, oNmerName.str());

		oNmerName.seekp(0);

		nmers.resize(0);
	}	

	cout << "Done.\n";

}

// Recursively fill a vector of strings with all possible non-degenerate nmers	
void fillVector (const int num_vars, const int currenPos, const string& growString, vector<string>& nmers) {

	if (currenPos == num_vars) {
		nmers.push_back(growString);
	} else {
		char nuc[4] = { 'A', 'C', 'G', 'T' };
		for (int i = 0; i < 4; ++i) {
			fillVector(num_vars, currenPos + 1, growString + nuc[i], nmers);
		}
	}
}


// Calculate nucleotide interval counts from sequence
void nucintConvert(const vector<vector<vector<string> > >& seqFrags, const int nucint_size, 
	const vector<string>& IDVec, const string outFile) {

	int num_int = 16; // every possible N..N
	char nuc[4] = { 'A', 'C', 'G', 'T' };

	ostringstream oNIname;	

	cout << "Counting nucleotide intervals of size " << flush;

	// for each size of nuc_int, up to the maximum
	for (int interval = 1; interval <= nucint_size; ++interval) {
		
		cout << interval << "..." << flush;

		// Build the 16 vectors that represent the degenerate nuc_ints
		vector<string> nuc_ints(num_int);
		for (int buildVec = 0; buildVec < num_int; ++buildVec) {
			nuc_ints[buildVec] += nuc[(int) (buildVec / 4)];
			for (int addN = 0; addN < interval; ++addN) {
				(nuc_ints[buildVec]).append("N");
			}
			nuc_ints[buildVec] += nuc[buildVec % 4];
		}

		oNIname << "nuc-int" << interval << '\0';

		// Call the countOligos function to count all of the nuc-int values
		oligoConvert(seqFrags, nuc_ints, IDVec, outFile, oNIname.str());

		oNIname.seekp(0);

		nuc_ints.resize(0);
	}
	cout << "Done.\n";

}

// Convert the MGA reps to counts
void MGAconvert(const vector<vector<vector<string> > >&seqFrags, const vector<Chromosome*>& MGAvec, 
				vector<map<string, double> >& ruleVec, const vector<string>& ruleNames,
				 const vector<string>& IDVec, const string& outname) {

	// Create a map that relates rule names to rule numbers
	map<string, int> nameToNum;
	for (int getRuleNum = 0; getRuleNum < ruleNames.size(); ++getRuleNum) {
		nameToNum[ruleNames[getRuleNum]] = getRuleNum;
	}

	cout << "Finding MGA Chromosomes..." << endl;
	
	int numSeq = seqFrags.size();
	int numWin = (seqFrags[0]).size();
	int numRep = (seqFrags[0][0]).size() - 1; // -1 since the first element is the real sequence fragment
	int numVar = (numRep > 0) ? numWin * 2: numWin; // * 2 if Z-scores are to be generated
	int mult = (numRep > 0) ? 2 : 1;
	int numMaps = MGAvec.size();
	int winSize = seqFrags[0][0][0].length();

	// Build the list of variables to be output
	ostringstream outStr;
	vector<string> MGAheaders;
	MGAheaders.push_back("PosNeg"); MGAheaders.push_back("SeqID"); 
	for (int MName = 0; MName < numMaps; ++MName) {
		string thisName = MGAvec[MName]->itsName();
		for (int winNames = 0; winNames < numWin; ++winNames) { // Create the list of window names
			outStr << thisName << "ws" << winSize - 1 << "n" << winNames << '\0';
			MGAheaders.push_back(outStr.str()); // Variable value
			if (numRep > 0) {
				outStr.seekp((int) outStr.tellp() - 1);
				outStr << "Z" << '\0';
				MGAheaders.push_back(outStr.str()); // Z-score
			}
			outStr.seekp(0);
		}
	}

	vector<vector<double> > MGAinds(numSeq);
	for (int setInds = 0; setInds < numSeq; ++setInds) {
		MGAinds[setInds].resize(numVar * numMaps);
	}

	vector<vector<int> > geneBounds(2);
	for (int currMap = 0; currMap < numMaps; ++currMap) { 
		if (currMap % 25 == 0) { cout << currMap << "..." << flush; }
		cout << currMap << "..." << flush;
		geneBounds[0].resize(MGAvec[currMap]->itsNumGenes());
		geneBounds[1].resize(MGAvec[currMap]->itsNumGenes());
		for (int getGeneBounds = 1; getGeneBounds < MGAvec[currMap]->itsNumGenes(); ++getGeneBounds) {
			geneBounds[0][getGeneBounds] = MGAvec[currMap]->minLen(getGeneBounds);
			geneBounds[1][getGeneBounds] = MGAvec[currMap]->maxLen(getGeneBounds);
		}
		int currBase = currMap * numVar;
		for (int curSeq = 0; curSeq < numSeq; ++curSeq) {
			for (int currWindow = 0; currWindow < numWin; ++currWindow) {
				double realVal = findMatches(MGAvec[currMap], seqFrags[curSeq][currWindow][0], geneBounds, ruleVec, nameToNum);
				MGAinds[curSeq][currBase + currWindow * mult] = realVal;
				if (numRep > 0) {
					double curCount = 0.0;
					double sum = 0.0;
					double sumSq = 0.0;
					for (int curRep = 1; curRep <= numRep; ++curRep) {
						curCount = findMatches(MGAvec[currMap], seqFrags[curSeq][currWindow][curRep], geneBounds, ruleVec, nameToNum);
						sum += curCount;
						sumSq += curCount * curCount;
					}
					MGAinds[curSeq][currBase + currWindow * mult + 1] = Zscore(realVal, sum, sumSq, numRep);
				}
			}
		}
	}
	cout << "Done.\n";
	indOutput(MGAheaders, MGAinds, IDVec, outname);
}
	
// Count the number of matches between el chromosome and el text string
double findMatches(const Chromosome* searchKey, const string& target, 
	const vector<vector<int> >& geneBounds, vector<map<string, double> >& ruleVec, 
	map<string, int>& nameToNum) {

	double numMatches = 0.0;

	int stopPos = searchKey->itsSize() - searchKey->maxLen();
	for (int startPos = 0; startPos < stopPos; ++startPos) {
		if (isFound(target, searchKey, 1, startPos, geneBounds, ruleVec, nameToNum)) {
			++numMatches;
			startPos += searchKey->minLen();
		}
	}
	return (double) numMatches;
}

// Is the pattern found from the current position?
bool isFound(const string& searchString, const Chromosome* searchKey, const int motifCheck, 
	const int seqPos, const vector<vector<int> >& geneBounds, 
	vector<map<string, double> >& ruleVec, map<string, int>& nameToNum) {

//if (motifCheck > 6) cout << motifCheck << flush;

	if (motifCheck == searchKey->itsNumGenes()) {  /* cout << "!" << flush; */ return true; } // Successfully found each component, trues will cascade backward
	
	// (else)
	int startPos = geneBounds[0][motifCheck];
	int stopPos = geneBounds[1][motifCheck];

	for (int i = startPos; i <= stopPos; ++i) {
		if (searchKey->isCompatible(searchString, seqPos, motifCheck, i, ruleVec, nameToNum)) {
			if (isFound(searchString, searchKey, motifCheck + 1, seqPos + i, geneBounds, ruleVec, nameToNum)) {
				return true;
			}
		}
	} 
	return false; 
}

// Calculate a Z-score, given the real value, sum, sum of squares and population size
double Zscore (const double total, const double sum, const double sumSq, const int popSize) {

	double mean = sum / popSize;

	double stdev = sqrt(sumSq / popSize - mean * mean);

	if (stdev == 0) { return 0; }

	double score = (total - mean) / stdev;

	return score;
}

// Function to write indices to an output file. Takes a header including PosNeg & TrainTest
void indOutput(const vector<string>& headers, const vector<vector<double> >& outIndices, 
	const vector<string>& IDVec, const string fileToOpen) {
	
	ofstream fOut(fileToOpen.c_str());

	// Print out the variable names
	int lastCol = headers.size();
	for (int i = 0; i < lastCol; ++i) {
		fOut << (headers[i]).c_str() << ',';
	}
	fOut << endl;
	
	// For each line of data, print the sequence ID followed by the list of variables
	int numSeq = IDVec.size();
	int lastVar = (outIndices[0]).size();
	for (int j = 0; j < numSeq; ++j) {
		fOut << IDVec[j] << ',';
		for (int k = 0; k < lastVar; ++k) {
			fOut << (outIndices[j])[k] << ',';
		}
		fOut << endl;
	}
	
	fOut.close();
}
