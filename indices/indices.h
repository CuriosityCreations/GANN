// indices.h: Header file for indices.cpp
// Copyright (C) 2001-2002 Robert G. Beiko

// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

#include"classChromosome.h"

using namespace std;

typedef map<string, double>::const_iterator CI;

const string USAGE = "Usage: indices seqfile outfile win-size overlap [-m mapfile]  [-o oligofile]  [-n num-nmers]  [-i num-nucints]  [-g MGAfile] [-r num-reps]\n";

// Functions that read and store input files
bool readSeqFile(const string& seqfile, vector<string>& seqVec, vector<string>& IDVec, int& numSeq);
bool readMapFile(const string& mapfile, vector<map<string, double> >& ruleVec, vector<string>& ruleNames, map<char,int>& nucDegen);
bool readOligoFile(const string& oligofile, vector<string>& oligoVec);

// Function to cut sequence into fragments and randomly shuffled replicates
void chopSequence(const vector<string>& seqVec, const int win_size, const int overlap, const int num_reps,
	vector<vector<vector<string> > >& seqFrags);

// Print the sequence fragments to a file
void seqOutput (const vector<vector<vector<string> > >& theFragments, const vector<string>& IDVec, const string& outfile);

// Take a degenerate n-mer, and add all of its non-degenerate members to a map using recursion
void buildMap(const pair<string, double>& toAdd, const int currentPos, const int finalPos, string buildNmer, map<string, double>& inProgress);
	
// Recursively fill a vector of strings with all possible non-degenerate nmers	
void fillVector (const int num_vars, const int currenPos, const string& growString, vector<string>& nmers);

// Read MGA Chromosomes from a file
void readMGAfile(ifstream& inFile, vector<Chromosome*>& MGAvec, const int winSize,
	vector<map<string, double> >& ruleVec, vector<string>& ruleNames);

// Find matches to a given MGA Chr
double findMatches(const Chromosome* searchKey, const string& target, 
	const vector<vector<int> >& geneBounds, vector<map<string, double> >& ruleVec, 
	map<string, int>& nameToNum);

// Is the pattern found from the current position?
bool isFound(const string& searchString, const Chromosome* searchKey, const int motifCheck, 
	const int seqPos, const vector<vector<int> >& geneBounds, 
	vector<map<string, double> >& ruleVec, map<string, int>& nameToNum);

// Function to shuffle a string
string seqShuffle (const string& toShuf, const int win_size);

// Functions that calculate and output indices
void mapConvert(const vector<vector<vector<string> > >& seqFrags, const vector<map<string, double> >& ruleVec, 
	const vector<string>& ruleNames, const vector<string>& IDVec, const string outFile); 
void oligoConvert(const vector<vector<vector<string> > >& seqFrags, const vector<string>& oligoVec, 
	const vector<string>& IDVec, const string outFile, const string outID); 
void nmerConvert(const vector<vector<vector<string> > >& seqFrags, const int nmer_size, 
	const vector<string>& IDVec, const string outFile); 
void nucintConvert(const vector<vector<vector<string> > >& seqFrags, const int nucint_size, 
	const vector<string>& IDVec, const string outFile); 
void MGAconvert(const vector<vector<vector<string> > >&seqFrags, const vector<Chromosome*>& MGAvec, 
	vector<map<string, double> >& ruleVec, const vector<string>& ruleNames,
	const vector<string>& IDVec, const string& outname);
	
// Count the number of hits within a given sequence
double enumerateHits (const CI& whichSearch, const string& toCalc);

// Search for a string of text containing exact nucleotides and 'N's
int degen_find(const string& theString, const string& theTarget);

// Calculate a Z-score, given the real value, sum, sum of squares and population size
double Zscore (const double total, const double sum, const double sumSq, const int popSize);

// Function to write indices to an output file
void indOutput(const vector<string>& headers, const vector<vector<double> >& outIndices, 
	const vector<string>& IDVec, const string fileToOpen);

