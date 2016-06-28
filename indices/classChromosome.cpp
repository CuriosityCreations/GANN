// classChromosome.cpp defines the Chromosome objects used by the genetic algorithm
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

// Modified from file 'classChromosome.cpp' created by Dr. Robert L. Charlebois, NeuroGadgets, Inc.

#include "classChromosome.h"

// Take a string representing an MGA Chromosome, and convert it to a real MGA Chromosome
Chromosome :: Chromosome(const string& chrRep, const int seqLen, vector<map<string, double> >& ruleVec, vector<string>& ruleNames) {

	seqLen_ = seqLen;
	clearScores();

	vector<string> chrBits = split(chrRep, ",");
	
	theName = chrBits[0];

	// Do the PosGene!
	int start, stop;
	if (chrBits[1][0] == '-') { start = 0; stop = seqLen - 1; } // '-' indicates complete coverage
	else {
		vector<string> startStop = split(chrBits[1], "-");
		start = atoi(startStop[0].c_str());
		stop = max(atoi(startStop[1].c_str()), seqLen_ - 1);
	}
	theGenes_.push_back(new posGene(start, stop));
	// Build the motif components
	vector<string> mComp = split(chrBits[3], "|");
	numComp = mComp.size();
cout << "Start." << endl;
	for (int i = 0; i < numComp; ++i) {
		if (mComp[i].find("*") != string::npos) { // Then it's a rangedContinuousGene
			theGenes_.push_back(new rangedContinuousGene(mComp[i], ruleVec, ruleNames));
		} else {
			if (mComp[i].find("/") != string::npos) { // Then it's a discreteGene
				theGenes_.push_back(new discreteGene(mComp[i]));
			} else {
				if (mComp[i].find("N") != string::npos) { // Then it's a gap!
					theGenes_.push_back(new continuousGene(mComp[i]));
				} else { 
					if (mComp[i].find("(X)") != string::npos) {
						; // Shouldn't have to do anything, this is a placeholder
					} else { cout << "Failed on " << mComp[i]; exit(-1);}
				}
			}
		}
	}
	numGenes_ = theGenes_.size();
	startReadHere_ = 0;
//	setConvertArray(startReadHere_);
	checkConstraints();
};

/*void Chromosome :: setConvertArray(const int readStart) {

	int motifLook;
	for (int i = 1; i < numComp + 1; ++i) {
		motifLook = i + readStart - 1;
		if (motifLook >= numGenes_) {
			motifLook -= (numGenes_ - 1);
		}
		searchPos[i] = motifLook;
	}
} */

Chromosome :: Chromosome(const Chromosome& rhs)
{

	sumScores_ = rhs.sumScores_;
	sumSquares_ = rhs.sumSquares_;
	numTrials_ = rhs.numTrials_;
	numGenes_ = rhs.numGenes_;
	theGenes_.resize(numGenes_);
	seqLen_ = rhs.seqLen_;
	maxMotifLength_ = rhs.maxMotifLength_;
	minMotifLength_ = rhs.minMotifLength_;

	// Create a new GeneBase* pointer for each Gene
	for (int i = 0; i < numGenes_; ++i) {
		theGenes_[i] = rhs.theGenes_[i]->clone();
	}
	if (!(isConsistent())) {
		cout << "Aaargh!" << endl;
	}
	startReadHere_ = rhs.startReadHere_;
//	setConvertArray(startReadHere_);
};

Chromosome :: ~Chromosome() {

	for (int i = 0; i < numGenes_; ++i) {
		delete theGenes_[i];
		theGenes_[i] = 0;
	}
}

Chromosome& Chromosome :: operator=(const Chromosome& rhs) {

	sumScores_ = rhs.sumScores_;
	sumSquares_ = rhs.sumSquares_;
	numTrials_ = rhs.numTrials_;
	numGenes_ = rhs.numGenes_;
	theGenes_.resize(numGenes_);
	startReadHere_ = rhs.startReadHere_;
//	setConvertArray(startReadHere_);

	// Create a new GeneBase* pointer for each Gene
	for (int i = 0; i < numGenes_; ++i) {
		theGenes_[i] = rhs.theGenes_[i]->clone();
	}
	
	return *this;
	
};

void Chromosome :: clearScores()
{
	sumScores_ = 0.0;
	sumSquares_ = 0.0;
	numTrials_ = 0;
};

// Check the validity of each motif component and the boundary
void Chromosome :: checkConstraints() {

	maxMotifLength_ = 0;
	minMotifLength_ = 0;

	// Calculate the minimum length of the motif
	for (int i = 1; i < numGenes_; ++i) {
		theGenes_[i]->fixIfIllegal();
		maxMotifLength_ += theGenes_[i]->itsMaxLength();
		minMotifLength_ += theGenes_[i]->itsMinLength();
	}
	theGenes_[0]->fixIfIllegal(maxMotifLength_, seqLen_);
}

bool Chromosome :: isConsistent() const {

	int testMin = 0, testMax = 0;
	// Calculate the minimum length of the motif
	for (int i = 1; i < numGenes_; ++i) {
		testMax += theGenes_[i]->itsMaxLength();
		testMin += theGenes_[i]->itsMinLength();
	}
	return (testMax == maxMotifLength_ && testMin == minMotifLength_);
}

// Return a text string representing all of the genes
void Chromosome :: toText (ofstream& outFile) const {

	theGenes_[0]->printComponent(outFile);
	outFile << "(" << startReadHere_ << ")|";
	int theIndex;
	for (int i = 1; i < numGenes_; ++i) {
		theIndex = i + startReadHere_ - 1;
		theGenes_[theIndex % numGenes_ + (int) (theIndex / numGenes_)]->printComponent(outFile);
	}
//	outFile << " F:" << itsFitness(); 
}

// Shuffle a vector
void vecShuf (vector<Chromosome*>& toShuf) {

	int swapWith = 0;	

	int len = toShuf.size();
	// Swap each element with another random element
	for (int i = 0; i < len; ++i) {
		swapWith = (int) rnd((double) len);
		if (swapWith != i) {
			Chromosome* temp = toShuf[i];
			toShuf[i] = toShuf[swapWith];
			toShuf[swapWith] = temp;
		}
	}
}
