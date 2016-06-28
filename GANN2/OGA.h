// OGA.h: header file for OGA.cpp
// GANN CORE
// Copyright (C) 2000-2004 Robert G. Beiko

// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

#include "MLP.h"

////////////////////////////////////////////////////////////////////
// CLASS OGA -------------------------------------------------------
//
//
//

class OGA : public SelectableBase {

public:

	// Constructors
	OGA(const int numSets, const int numSuperSets, const int numChr, const int sizeD, 
		const int sizeB, const int sizeI, const int sizeC, const std::vector<double>& range,
		const bool isGA, const bool IVO, const string inFileName, long& setSeed, const pStr& thePars);
		
	// Destructor
	~OGA(); 
	
	// Create input arrays
	void setInputs(const string inFileName);

	// Output methods
	void startRun(ofstream& outFile);
	void genOutput(ofstream& outFile, const int generation);

	// Check for best or worst times
	void checkTime (const long theTime) { if (theTime > worstTime) worstTime = theTime;
										  if (theTime < bestTime) bestTime - theTime;
										}

	// Run a generation
	void execute(const short OGAgen, const std::vector<double> GARange, 
		const int numParams, const int cNum, const pStr& thePars);

	// Run a single member in a generation
	void OCrun (const short OGAgen, std::vector<double>& currentGAPar, const int numParams, const pStr& thePars);
	// Evolve
	virtual void optimizeParameters (const std::vector<double> GARange);

	void updateBestNets(const vector<bool>& whichReplace);

	void bestNetOutput(ofstream& outFile);

private:

	void createNN(std::vector<double>& currentGAPar, const int numParams, const int OGAgen, vector<int>& inputsToUse, const pStr& thePars); // Determine which inputs to use

	ofstream saveNetFile; // Save NNs to file

	bool isGANN; // True if GANN, false if BPNN
	bool vecOp; // True if input vector optimization is used
	double eTable[ETABLESIZE * 2]; // Table of e^x values
	MLP* neuralNet; // Pointer to the neural network
	Chromosome* toSave; // The best Chromosome
	oCscore* chrScores;  // Scores of outer Chromosomes
	
	int totalOC;
	int OGArounds;
	
	double recRate, recRepl;
	double mutRate, mutAmt, mutProp;
	double killFrac, killDiff;

	double bestTime;
	double worstTime;
	
	int BPrepl, GAnumRec;

	// Store the best Chromosome
	vector<int> bestNetConfig[3]; // The configuration of the best nets
	vector<vector<double> > bestNet; // The best nets
	vector<double> best; // The best scores

	// Training and test arrays
	int howManyCategories; // = Number of categories to classify into
	int howManySets; // = Number of categories * 2
	int numInputs;
	vector<int> setSizes; // Sizes of all sets
	vector<string> titles; // Names of all inputs
	vector<vector<vector<double> > > sets; // The training and test sets
						  // Index: 0 = neg. train
						  //        1 = pos. train
						  //        2 = neg. test
						  // 	    3 = pos. test
};
