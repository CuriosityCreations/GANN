// classChromosome.h: header file for clasChromosome.cpp
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

// Modified from file 'classChromosome.h' created by Dr. Robert L. Charlebois, NeuroGadgets, Inc.

#include "classGene.h"

class Chromosome {
private:

	double fitness_;
	double sumScores_;	// used in evaluating fitness
	double sumSquares_;	// used in evaluating fitness
	int numTrials_;		// used in evaluating fitness
	int numGenes_;
	int seqLen_;
	int maxMotifLength_;
	int minMotifLength_;
	int startReadHere_; // In a circular Chromosome representation, start reading the first gene at this position
	
	int numComp;

	string theName;

//	int searchPos[MOTIF_SIZE + 1]; // Converts the requested position into the correct motif component (required by startReadHere_;)

	vector<GeneBase*> theGenes_;

public:

	Chromosome(const string& chrRep, const int seqLen, vector<map<string, double> >& ruleVec, vector<string>& ruleNames);

//	void setConvertArray(const int readStart);

	Chromosome(const Chromosome& rhs);

	~Chromosome();
	Chromosome& operator=(const Chromosome& rhs);
	
	int itsNumGenes() const { return numGenes_; };
	int itsNumTrials() const { return numTrials_; };

	int itsOrigin() const { return startReadHere_; };
	int itsStart() const { return theGenes_[0]->itsMinLength(); };
	int itsStop() const { return theGenes_[0]->itsMaxLength(); };
	int maxLen() const { return maxMotifLength_; };
	int minLen() const { return minMotifLength_; };
	int minLen(const int index) const { return theGenes_[index]->itsMinLength(); };
	int maxLen(const int index) const { return theGenes_[index]->itsMaxLength(); };
	int itsSize() const { return seqLen_; };
	string itsName() const { return theName; };
	
	int convertedIndex(const int position) const { return (position + startReadHere_ <= numGenes_) ? (startReadHere_ + position - 1) : (startReadHere_ + position - numGenes_); };

	// Is the motif found at a given point in the sequence?
	bool isCompatible(const string& searchString, const int seqPos, const int motifToCheck, const int lengthToCheck,
		vector<map<string, double> >& ruleVec, map<string, int>& nameToNum) const
		{ return theGenes_[motifToCheck]->isCompatible(searchString, seqPos, lengthToCheck, ruleVec, nameToNum); };
	
	GeneBase* theGene(const int index) { return theGenes_[index]; };

	void clearScores();
		// Zeroes the scoring information

	void clearGenes() { for (int i = 0; i < numGenes_; ++i) { delete theGenes_[i]; theGenes_[i] = 0; } };

	void computeFitness();

	// Check the validity of each motif component and the boundary
	void checkConstraints();
	
	bool isConsistent() const;

	// Print the Chromosome's genes to the ofstream
	void toText (ofstream& outFile) const;

	GeneBase* cloneGene(const int index) { return theGenes_[index]->clone(); };
	
	friend class Evolvable;

};

void vecShuf (vector<Chromosome*>& toShuf);

