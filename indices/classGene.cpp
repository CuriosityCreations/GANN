// classGene.cpp defines various types of Gene objects used by Chromosomes in the genetic algorithm.
// GANN CORE
// Copyright (C) 2000-2002 Robert G. Beiko

// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

#include"classGene.h"

// ------------- POSGENE ------------------
posGene :: posGene(const int seqLen) : GeneBase() {

	firstPos = 0;
	lastPos = seqLen - 1;
	
}

void posGene :: fixIfIllegal(const int maxMotif, const int seqLen) {

	if (lastPos < firstPos) { // If the two variables are in the wrong order, then swap them
		int temp = lastPos;
		lastPos = firstPos;
		firstPos = temp;
	}
	
	int excess = lastPos + maxMotif - seqLen + 1;
	if (excess > 0) { // Prevent the motif search from passing the end of the sequence
		lastPos -= excess;
//		firstPos -= excess;
if (lastPos < firstPos) firstPos = lastPos;
	}
	if (firstPos < 0) {
		firstPos = 0;
	}
}

void posGene :: printComponent(ofstream& outfile) {

	outfile << "StartPoint" << firstPos << "-" << lastPos << "|";
	
}

// --------------------------- MOTIFGENE: ABSTRACT CLASS --------------------------

// -------------------------------- DISCRETEGENE ----------------------------------

// Create from a text string
discreteGene :: discreteGene(const string& createFrom) {

	vector<string> bases = split(createFrom, "/");
	minLength = bases.size() - 1; // Subtract because of the threshold value
	maxLength = minLength;
	threshold = atof(bases[minLength].c_str());
	for (int i = 0; i < minLength; ++i) {
		if (bases[i].find(" ") != string::npos) { // Then it's a weightedBaseRep
			baseSeq.push_back(new weightedBaseRep(bases[i]));
		} else { 
			if (bases[i].length() == 1) { // Then it's a discreteBaseRep
				baseSeq.push_back(new degenBaseRep(bases[i]));
			} else { // Huh?
				die("Failed on baseRep creation at " + bases[i]);
			}
		}
	}
}

/*discreteGene :: discreteGene(const discreteGene& rhs) { 

	vector<baseRep*> copyRep = rhs.repAtPos();
	setLength(rhs.minLength, rhs.maxLength); 
	for (int i = 0; i < rhs.itsMinLength(); ++i) {
		baseSeq.push_back((copyRep[i])->clone());
	}
	threshold = rhs.threshold;
} */

discreteGene :: ~discreteGene() {

	for (int i = 0; i < minLength; ++i) {
		delete baseSeq[i];
		baseSeq[i] = 0;
	}
}

// Determine whether a subsequence is compatible with this motif component
bool discreteGene :: isCompatible(const string& seq, const int startPt, 
		const int lengthToCheck, vector<map<string, double> >& ruleVec, 
		map<string, int>& nameToNum) const {

	if (lengthToCheck == 0) { return true; }

// cout << "D" << flush;
	double runningTotal = 0.0;

	for (int i = 0; i < lengthToCheck; ++i) { 
		runningTotal += baseSeq[i]->isCompatibleWith(seq[startPt + i]);
	}

	return (runningTotal >= threshold) ? true : false;
};

// Print the base sequence
void discreteGene :: printComponent(ofstream& outfile) {

	if (minLength == 0) { 
		outfile << "(X)";
	} else {
		for (int i = 0; i < minLength; ++i) {
			baseSeq[i]->printComponent(outfile);
			outfile << "/";
		}
	}
	outfile << threshold;
	outfile << "|";
}

void discreteGene :: fixIfIllegal() {

//	if (threshold < MIN_THRESHOLD * (double) minLength) {
//		threshold = minLength;
//	}
}

// -------------------------------- CONTINUOUSGENE ---------------------------

// Create from a string of the form 'N1-1' (negative lengths are not allowed)
continuousGene :: continuousGene(const string& createFrom) {

	string toChange = createFrom;
	toChange.erase(0,1);
	vector<string> minMax = split(toChange, "-");
	minLength = atoi(minMax[0].c_str());
	maxLength = atoi(minMax[1].c_str());

}

// Check the boundary conditions
void continuousGene :: fixIfIllegal() {
	
/*	if (minLength > MAX_GAP_LEN) minLength = MAX_GAP_LEN; 
	else if (minLength < 0) minLength = 0; 
	if (maxLength > MAX_GAP_LEN) maxLength = MAX_GAP_LEN; 
	else if (maxLength < 0) maxLength = 0;
	
	if (minLength > maxLength) { 
		int temp = minLength;
		minLength = maxLength;
		maxLength = temp;
	}
	int diff = maxLength - minLength;
	if (diff > MAX_CONT_RANGE) {
		if (maybe) {
			maxLength -= (diff - MAX_CONT_RANGE);
		} else {
			minLength += (diff - MAX_CONT_RANGE);
		}
	}	*/
}

// ------------------------------ RANGEDCONTINUOUSGENE ------------------------

rangedContinuousGene :: rangedContinuousGene(const string& createFrom, 
	vector<map<string, double> >& ruleVec, vector<string>& ruleNames) {

	vector<string> components = split(createFrom, " ");
	
cout << "RCG." << endl;

	// Get the name
	ruleName = components[0];

	bool nameFound = false;
	for (int i = 0; i < ruleNames.size(); ++i) {
		if (ruleName == ruleNames[i]) {
			ruleType = i;
			ruleSize = (int) pow(ruleVec[i].size(), 0.25); // The fourth root
			nameFound = true;
		}
	}
	if (!nameFound) { cout << "Rule name " << ruleName << " not in input file.\n"; exit(-1); }

	// Get the length
	components[1].erase(0,1);
	vector<string> minMax = split(components[1], "-"); // OK since length >= 0
	minLength = atoi(minMax[0].c_str());
	maxLength = atoi(minMax[1].c_str());
	
	// Get the range (careful! Can be negative)
	components[3].erase(0,1);
	vector<string> range = split(components[3], "-");
	if (components[3].find("--") != string::npos) { // Both values are negative
		minVal = atof(range[0].c_str()) * -1;
		maxVal = atof(range[1].c_str()) * -1;
	} else { // Second value is positive
		minVal = atof(range[0].c_str()) * ((components[3][0] == '-') ? -1 : 1);
		maxVal = atof(range[1].c_str());
	}
cout << "RuleName " << ruleName << " MinLength: " << minLength << " MaxLength: " << maxLength << " MinVal: " << minVal << " MaxVal: " << maxVal << endl;
}

// Determine whether the specified subsequence is compatible with the motif component
bool rangedContinuousGene :: isCompatible(const string& seq, const int startPt, 
		const int lengthToCheck, vector<map<string, double> >& ruleVec, 
		map<string, int>& nameToNum) const {

//cout << ruleName << " " << ruleType << " ";
// cout << "R" << flush;
// if (ruleName == "BRUNAKPROP") cout << "(" << lengthToCheck << ")" << flush;
	if (lengthToCheck < ruleSize) { return true; } // We can't check it in this case, so return 'ok'

	double runningTotal = 0.0;
	int totalLen = lengthToCheck - ruleSize + 1;
	int stopPoint = startPt + totalLen;	

	for (int i = startPt; i < stopPoint; ++i) {
//		runningTotal += ruleVec[nameToNum[ruleName]][seq.substr(i, ruleSize)];
		runningTotal += ruleVec[ruleType][seq.substr(i, ruleSize)];
//if (ruleName == "BRUNAKPROP") cout << "{" << seq.substr(i, ruleSize) << "}" << flush;
	}
//if (ruleName == "BRUNAKPROP") cout << endl;

	double stupidVar = runningTotal / (double) totalLen;

	cout << flush;
//	runningTotal /= (double) totalLen;
//cout << runningTotal << " " << totalLen << " " << minVal << " " << maxVal;
//cout << "[" << runningTotal << "]" << endl;
//cout << flush;
// if (ruleName == "BRUNAKPROP") if (lengthToCheck == 4) cout << "{" << seq.substr(startPt, 4) << "}" << endl; //<< "[" << runningTotal << "]" << endl;
/*	if (stupidVar >= minVal && stupidVar <= maxVal) { 
if (ruleName == "BRUNAKPROP") cout << endl;
		return true;
	} 
	return false; */
	return (stupidVar >= minVal && stupidVar <= maxVal);
}

// Ensure that the definition of rangedContinuousGene is legal
void rangedContinuousGene :: fixIfIllegal() {

	continuousGene::fixIfIllegal(); // Check the length

	if (minVal > maxVal) { // Check the boundary values
		double temp = minVal;
		minVal = maxVal;
		maxVal = temp;
	}
	if (fabs(minVal) < 0.0001) { minVal = 0.0; }
	if (fabs(maxVal) < 0.0001) { maxVal = 0.0; }
}

// Print the motif component
void rangedContinuousGene :: printComponent(ofstream& outfile) {

	outfile << ruleName << " L" << minLength << "-" << maxLength << " * V" << minVal << "-" << maxVal << "|";
	
}
