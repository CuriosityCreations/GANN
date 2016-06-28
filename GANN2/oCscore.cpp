// oCscore.cpp defines the objects that store the training and testing scores achieved by neural networks.
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

#include"oCscore.h"

// CLASS OCSCORE
oCscore :: oCscore (const int numInnerRounds) {

	totalRounds = numInnerRounds;	

	for (int j = 0; j < 4; ++j) bestScores[j] = 0.0;
	for (int i = 0; i < 5; ++i) oCoutputs[i].open((OCOUTFILES[i]).c_str());
}

oCscore :: ~oCscore () {

	for (int i = 0; i < 5; ++i) oCoutputs[i].close();
}

void oCscore :: startOfRound(const int numIn, const int numHid, const int numOut, 
	const int theChr, const int gen, const int howBig) {
	
	for (int i = 0; i < 4; ++i) { (bestRunsOfOC[i]).resize(howBig); bestScoresOfOC[i] = 0.0; }

	chrConfig[0] = numIn;
	chrConfig[1] = numHid;
	chrConfig[2] = numOut;
	size = numIn * numHid + numHid * numOut;
		
	theTime = 1000000000; // This will (hopefully) change quickly	
	chrID = theChr;
	currentRound = gen;	

	// Initialize the scoring array
	for (int j = 0; j < 4; ++j) {
		scoreDev[j].resize(howBig);
		for (int k = 0; k < howBig; ++k) scoreDev[j][k] = OC_SENTRY;
	}
}

void oCscore :: endOfRound(const std::vector<double>& currentGAPar, const vector<int>& inputsToUse) {

	// Output scores
	chrScoresInit();
	chrScoresOutput(true);

	// Output NN parameters
	paramOutput(currentGAPar, inputsToUse);
}

// paramOutput sends the identifiers, parameters and best scores and time of each to outFile
void oCscore :: paramOutput(const std::vector<double>& paramToOutput, const vector<int>& whichInputs) {

	oCoutputs[4] << currentRound << "," << chrID << ","; // Identifiers
	for (int i = 0; i < nParams; ++i) oCoutputs[4] << paramToOutput[i] << ","; // GA parameters
	for (int k = 0; k < 4; ++k) {
		oCoutputs[4] << bestScores[k] << ",";
	}
	oCoutputs[4] << "," << theTime << ","; // Best generalization scores and time
	for (int j = 0; j < whichInputs.size(); ++j) {
		oCoutputs[4] << whichInputs[j] << ",";
	}
	oCoutputs[4] << endl;
}

bool oCscore :: ifNewBestTime(const long aTime) {

		if (aTime < theTime) {
		theTime = aTime;
		return true;
	}
	return false;
}

// Sort a copy of the scoreDev array, and return the nth best element
double oCscore :: getBestScores(const int type, const int nthBest) const {

	vector<double> copyOf = scoreDev[type];
	std::sort(copyOf.begin(), copyOf.end());
	return copyOf[copyOf.size() - nthBest - 1];

}



// Record the best generalization scores of the neural network, return true if max generalization is achieved
bool oCscore :: recordScores(const double minGen) {
	
	for (int checkBest = 0; checkBest < 4; ++checkBest) {
		if (roundScores[checkBest] > bestScores[checkBest]) {
			bestScores[checkBest] = roundScores[checkBest];
		}
	}
	if ((bestScores[2] == roundScores[2]) && (bestScores[2] > minGen)) return true;
	return false;
} // oCscore :: recordScores
///////////////////////////////////////////////////////////////////////
// CLASS BPSCORE
//
//
//
BPscore :: BPscore (const int numInnerRounds) : oCscore (numInnerRounds) {

	nParams = N_BP_PARAMS;

	for (int i = 0; i < 4; ++i) {
		scoreDev[i].resize(totalRounds);
		for (int j = 0; j < totalRounds; ++j) {
			scoreDev[i][j] = 0.0;
		}
	}
}

BPscore :: ~BPscore () {}

void BPscore :: nextConfig(const std::vector<double>& params, const int chrRound, const int Ogen, const int inNum) {	// Initialize the best scores array

	int howBig = totalRounds;
	startOfRound(inNum, (int) params[XHIDNODE], (int) params[XOUTNODE], chrRound, Ogen, howBig);
}

void BPscore :: reset() {
	innerRound = 0;
	for (int i = 0; i < 4; ++i) {
		bestScores[i] = 0.0;
		roundScores[i] = 0.0;
		for (int j = 0; j < totalRounds; ++j) {
			scoreDev[i][j] = OC_SENTRY;
		}
	}
}
	

// paramInit sends the output headers for paramOutput to outFile
void BPscore :: paramInit() {

	oCoutputs[4] << "Current oGA Round, which oC, ";
	for (int i = 0; i < nParams; ++i) oCoutputs[4] << BP_PARAM_NAMES[i] << ",";
	oCoutputs[4] << "trainFP, trainDist, testFP, testDist, bestTime";
	oCoutputs[4] << endl;
}

// chrScoresInit sends the headings for chrScoresOutput to outFile
void BPscore :: chrScoresInit() {

	for (int i = 0; i < 4; ++i) {
		oCoutputs[i] << "Current oGA Round, which oC, ";
		oCoutputs[i] << endl;
	}
}

// chrScoresOutput sends the identifiers of a given set of scores to outFile
void BPscore :: chrScoresOutput(bool isEndOfRound) {

	double scoreOut = OC_SENTRY;

	for (int i = 0; i < 4; ++i) {
		oCoutputs[i] << currentRound << "," << chrID << ","; // Identifiers
		double maxScore = OC_SENTRY;
		for (int j = 0; j < totalRounds; ++j) { 
			scoreOut = (isEndOfRound) ? bestRunsOfOC[i][j] : scoreDev[i][j];
			if (scoreOut > maxScore) { maxScore = scoreOut; }
			if (scoreOut != OC_SENTRY) { oCoutputs[i] << scoreOut << ","; } else break; // No more scores to output
		}
		
		// Print out the best gen score
		oCoutputs[i] << "," << maxScore;

		// Scores for all training rounds
		oCoutputs[i] << endl;

		if (!isEndOfRound) {
			// Is this the best score from this round? If so, save it in bestOfOC
			if (maxScore > bestScoresOfOC[i]) {
				bestScoresOfOC[i] = maxScore;
				bestRunsOfOC[i] = scoreDev[i];
			}
		}
	}
}

void BPscore :: findBestScore(const bool isTest) {

	int inda, indb;
	if (isTest) {
		inda = 2;
		indb = 3;
	} else {
		inda = 0;
		indb = 1;
	}
	roundScores[inda] = scoreDev[inda][innerRound];
	roundScores[indb] = scoreDev[indb][innerRound];
}
//////////////////////////////////////////////////////////////////////
// CLASS GASCORE
//
//
GAscore :: GAscore (const int numInnerRounds) : oCscore (numInnerRounds) {

	nParams = N_GA_PARAMS;
}

GAscore :: ~GAscore () {}

void GAscore :: nextConfig(const std::vector<double>& params, const int chrRound, const int Ogen, const int inNum) {	// Initialize the best scores array

	numEvo = (int) params[XNUMEVO];
	numChr = (int) params[XNUMCHR] / (int) params[XNUMEVO]; // Number of Chromosomes PER EVOLVABLE
	
	const int howBig = numChr * numEvo;

	startOfRound(inNum, (int) params[XHIDNODE], (int) params[XOUTNODE], chrRound, Ogen, howBig);
	
}

// Reset scores
void GAscore :: reset() {
	innerRound = 0;
	for (int i = 0; i < 4; ++i) {
		bestScores[i] = 0.0;
		roundScores[i] = 0.0;
	}
}

// paramInit sends the output headers for paramOutput to outFile
void GAscore :: paramInit() {

	oCoutputs[4] << "Current oGA Round, which oC, ";
	for (int i = 0; i < nParams; ++i) oCoutputs[4] << GA_PARAM_NAMES[i] << ",";
	oCoutputs[4] << "trainFP, trainDist, testFP, testDist, bestTime";
	oCoutputs[4] << endl; 
}

// chrScoresInit sends the headings for chrScoresOutput to outFile
void GAscore :: chrScoresInit() {

	int howMany = numEvo * numChr;

	for (int type = 0; type < 4; ++type) {
		oCoutputs[type] << "Current oGA Round, which oC, ";
		for (int i = 0; i < 10; ++i) oCoutputs[type] << i << ",";
		oCoutputs[type] << ",avg";
		oCoutputs[type] << endl;
	} 
}

// chrScoresOutput sends the identifiers of a given set of scores to outFile
void GAscore :: chrScoresOutput(bool isEndOfRound) {

	int howBig = numChr * numEvo;
	for (int i = 0; i < 4; ++i) {
		vector<double> outVec = (isEndOfRound) ? bestRunsOfOC[i] : scoreDev[i];
		std::sort(outVec.begin(), outVec.end(), std::greater<double>());
		double average = 0.0;
		for (int q = 0; q < howBig; ++q) {
			average += outVec[q];
		}
		average /= howBig;
		oCoutputs[i] << currentRound << "," << chrID << ","; // Identifiers
		for (int j = 0; j < 10; ++j) oCoutputs[i] << outVec[j] << ","; // Scores for all iCs
		oCoutputs[i] << "," << average;
		oCoutputs[i] << endl;

		if (!isEndOfRound) {
			if (average > bestScoresOfOC[i]) {
				bestScoresOfOC[i] = average;
				bestRunsOfOC[i] = outVec;
			}
		}
	}	
}

void GAscore :: findBestScore(const bool isTest) {

	int inda, indb;
	int howBig = numChr * numEvo;
	if (isTest) {
		inda = 2;
		indb = 3;
	} else {
		inda = 0;
		indb = 1;
	}
	for (int i = 0; i < howBig; ++i) {
		if (scoreDev[inda][i] > roundScores[inda]) { roundScores[inda] = scoreDev[inda][i]; }
		if (scoreDev[indb][i] > roundScores[indb]) { roundScores[indb] = scoreDev[indb][i]; }
	}
}
