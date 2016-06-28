// classChromosome.cpp defines the Chromosome objects used by the genetic algorithm
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

// Modified from file 'classChromosome.cpp' created by Dr. Robert L. Charlebois, NeuroGadgets, Inc.

#include "classChromosome.h"

#define DEBUG_ON false

const bool MY_RAN = true;

Chromosome :: Chromosome(const int sizeD, const int sizeB, const int sizeI, const int sizeC)
{
	commonInit(sizeD, sizeB, sizeI, sizeC);
};

Chromosome :: Chromosome(const int sizeD, const int sizeB, const int sizeI, const int sizeC, const pStr& thePars)
{
	Btrue = elem(thePars.intGlob,"NUM_INPUTS");
	commonInit(sizeD, sizeB, sizeI, sizeC);
};

Chromosome :: Chromosome(const int sizeD, const int sizeB, const int sizeI, const int sizeC, const double& loBound, const
	double& hiBound, long& setSeed, const pStr& thePars)
{

	Btrue = elem(thePars.intGlob,"NUM_INPUTS");

	commonInit(sizeD, sizeB, sizeI, sizeC);
	int j = 0;
	for (int newGenes = 0; newGenes < boolStart_; ++newGenes) {
		j = newGenes * 2;
		geneDoubleArray_[newGenes].rangeStart(loBound, hiBound);
	}
	int totalTrue = 0;
	vector<int> trueGenes(sizeB);
	for (int ii = 0; ii < sizeB; ++ii) {
		trueGenes[ii] = 0;
	}
	int current = 0;
	
	// Select NUM_INPUTS random genes to be true if sizeB is > 0
	if (sizeB >= Btrue) {
		while (totalTrue < Btrue) {
			current = (int) (ran1(&setSeed) * sizeB);
			if (trueGenes[current] != 1) {
				trueGenes[current] = 1;
				++totalTrue;
			}
		}
	}
		
	// Randomly select NUM_INPUTS worth of genes to make positive
	for (int newBool = 0; newBool < numBool_; ++newBool) {
		if (trueGenes[newBool] == 0) { // If not an input
			geneBoolArray_[newBool].forceVal(0.0);
		} else { // If an input
			geneBoolArray_[newBool].forceVal(1.0);
		}
	}
};



Chromosome :: Chromosome(const int sizeD, const int sizeB, const int sizeI, const int sizeC, const std::vector<double>&
	range, long& setSeed, const pStr& thePars)

{

	Btrue = elem(thePars.intGlob,"NUM_INPUTS");

	commonInit(sizeD, sizeB, sizeI, sizeC);

	int j = 0;

	for (int newGenes = 0; newGenes < boolStart_; ++newGenes) {
		j = newGenes * 2;
		geneDoubleArray_[newGenes].rangeStart(range[j], range[j+1]);
	}

	for (int newBool = 0; newBool < numBool_; ++newBool) {
			geneBoolArray_[newBool].ranStart();
	}

	int totalTrue = 0;
	vector<int> trueGenes(sizeB);
	for (int ii = 0; ii < sizeB; ++ii) {
		trueGenes[ii] = 0;
	}	
	
	int current = 0;
	
	// Select NUM_INPUTS random genes to be true
	if (sizeB >= Btrue) {
		while (totalTrue < Btrue) {
			current = (int) (ran1(&setSeed) * sizeB);
			if (trueGenes[current] != 1) {
				trueGenes[current] = 1;
				++totalTrue;
			}
		}
	}		
	
	// Randomly select NUM_INPUTS worth of genes to make positive
	for (int newBool = 0; newBool < numBool_; ++newBool) {
		if (trueGenes[newBool] == 0) { // If not an input
			geneBoolArray_[newBool].forceVal(0.0);
		} else { // If an input
			geneBoolArray_[newBool].forceVal(1.0);
		}
	}

};


Chromosome :: Chromosome(const int sizeD, const int sizeB, const int sizeI, const int sizeC, std::ifstream& inFile) throw(FileErr)
{
	commonInit(sizeD, sizeB, sizeI, sizeC); 
};

Chromosome :: Chromosome(const Chromosome& rhs)
{

	commonInit(rhs.numDouble_, rhs.numBool_, rhs.numInt_, rhs.numChar_);
	sumScores_ = rhs.sumScores_;
	sumSquares_ = rhs.sumSquares_;
	numTrials_ = rhs.numTrials_;

//	gene_ = rhs.gene_;
	int i;
	for (i = 0; i < numDouble_; ++i) {
		geneDoubleArray_[i].forceVal(rhs.geneDoubleArray_[i].getVal());
	}
	for (i = 0; i < numBool_; ++i) {
		geneBoolArray_[i].forceVal(rhs.geneBoolArray_[i].getVal());
	}
	for (i = 0; i < numInt_; ++i) {
		geneIntArray_[i].forceVal(rhs.geneIntArray_[i].getVal());
	}
	for (i = 0; i < numChar_; ++i) {
		geneCharArray_[i].forceVal(rhs.geneCharArray_[i].getVal());
	}
};

// This was written specifically to save values from a BPNet
Chromosome :: Chromosome(const int num_double, const double * const the_links) {

	commonInit(num_double, 0, 0, 0);
	
	for (int i = 0; i < num_double; ++i) {
		geneDoubleArray_[i].forceVal(the_links[i]); // Force the linkvals
	}
}

Chromosome :: ~Chromosome() {

	if (boolStart_ != 0) {
		delete [] geneDoubleArray_;
		geneDoubleArray_ = 0;
	}
	if (intStart_ != boolStart_) {
		delete [] geneBoolArray_;
		geneBoolArray_ = 0;
	}
	if (charStart_ != intStart_) {
		delete [] geneIntArray_;
		geneIntArray_ = 0;
	}
	if (numGenes_ != charStart_) {
		delete [] geneCharArray_;
		geneCharArray_ = 0;
	}
}

void Chromosome :: clearGenes() {

	int i;
	
	for (i = 0; i < numDouble_; ++i) {
		geneDoubleArray_[i].forceVal(0.0);
	}
	for (i = 0; i < numBool_; ++i) {
		geneBoolArray_[i].forceVal(0.0);
	}
	for (i = 0; i < numInt_; ++i) {
		geneIntArray_[i].forceVal(0.0);
	}
	for (i = 0; i < numChar_; ++i) {
		geneCharArray_[i].forceVal(0.0);
	}
}

Chromosome& Chromosome :: operator=(const Chromosome& rhs) throw(Err)
{

	#if DEBUG_ON

		if (numGenes_ != rhs.numGenes_) throw Err("Chromosome::operator=(), incompatible Chromosomes");
	#endif

	Btrue = rhs.Btrue;

	sumScores_ = rhs.sumScores_;
	sumSquares_ = rhs.sumSquares_;
	numTrials_ = rhs.numTrials_;

	int i;

	for (i = 0; i < numDouble_; ++i) {
		geneDoubleArray_[i].forceVal(rhs.geneDoubleArray_[i].getVal());
	}
	for (i = 0; i < numBool_; ++i) {
		geneBoolArray_[i].forceVal(rhs.geneBoolArray_[i].getVal());
	}
	for (i = 0; i < numInt_; ++i) {
		geneIntArray_[i].forceVal(rhs.geneIntArray_[i].getVal());
	}
	for (i = 0; i < numChar_; ++i) {
		geneCharArray_[i].forceVal(rhs.geneCharArray_[i].getVal());
	}

	return *this;
};

void Chromosome :: clearScores()
{
	sumScores_ = 0.0;
	sumSquares_ = 0.0;
	numTrials_ = 0;
};

void Chromosome :: commonInit(const int sizeD, const int sizeB, const int sizeI, const int sizeC)

{
try {

	numDouble_ = sizeD;
	numBool_ = sizeB;
	numInt_ = sizeI;
	numChar_ = sizeC;

	boolStart_ = sizeD;
	intStart_ = boolStart_ + sizeB;
	charStart_ = intStart_ + sizeI;
	int size = charStart_ + sizeC;
	numGenes_ = size;
	
	if (sizeD != 0) geneDoubleArray_ = new GeneDouble[sizeD];
	else geneDoubleArray_ = 0;
	if (sizeB != 0) geneBoolArray_ = new GeneBool[sizeB];
	else geneBoolArray_ = 0;
	if (sizeI != 0) geneIntArray_ = new GeneInt[sizeI];
	else geneIntArray_ = 0;
	if (sizeC != 0) geneCharArray_ = new GeneChar[sizeC];
	else geneCharArray_ = 0;

} catch (std::bad_alloc) {
	cout << "bad_alloc!!!" << endl;
} catch (...) {
	cout << "not bad_alloc." << endl;
};
	clearScores();

};


void Chromosome :: replaceGenes(const GeneDouble* const & theD, const GeneBool* const & theB, 
					  const GeneInt* const & theI, const GeneChar* const & theC) {
					  
	int i;
	for (i = 0; i < numDouble_; ++i) {
		geneDoubleArray_[i].forceVal(theD[i].getVal()); 
	} 
	for (i = 0; i < numBool_; ++i) {
		geneBoolArray_[i].forceVal(theB[i].getVal()); 
	}	
	for (i = 0; i < numInt_; ++i) {
		geneIntArray_[i].forceVal(theI[i].getVal()); 
	}	
	for (i = 0; i < numChar_; ++i) {
		geneIntArray_[i].forceVal(theC[i].getVal()); 
	}
}


void Chromosome :: computeFitness() throw(Err)

{

	#if DEBUG_ON

		if (numTrials_ < 2) throw Err("Chromosome::computeFitness(), numTrials_ < 2");

	#endif



	double mean = sumScores_ / numTrials_;

	fitness_ = mean - sqrt((sumSquares_ - sumScores_ * mean) / (numTrials_ - 1));

};



void Chromosome :: mutateGene(int which, double howMuch) throw(Err)
{

	#if DEBUG_ON
		if (which >= numGenes_) throw Err("Chromosome::mutateGene(), invalid gene selected.");
		if (howMuch <= 1.0) throw Err("Chromosome::mutateGene(), howMuch must be > 1.0");
	#endif

	if (which < boolStart_) {
		geneDoubleArray_[which].mutate(howMuch); 
	} else { 
		if (which < intStart_) {
			geneBoolArray_[which - boolStart_].mutate(howMuch);
		} else {
			if (which < charStart_) {
				geneIntArray_[which - intStart_].mutate(howMuch);
			} else geneCharArray_[which - charStart_].mutate(howMuch);
		}
	}
};



void Chromosome :: mutateGene(int which, const double& loBound, const double& hiBound) throw(Err)

{

	#if DEBUG_ON

		if (which >= numGenes_) throw Err("Chromosome::mutateGene(), invalid gene selected.");

	#endif

	if (which < boolStart_) {
		geneDoubleArray_[which].mutate(loBound, hiBound); 
	} else {
		if (which < intStart_) {
			geneBoolArray_[which - boolStart_].mutate(loBound, hiBound);
		} else {
			if (which < charStart_) {
				geneIntArray_[which - intStart_].mutate(loBound, hiBound);
			} else geneCharArray_[which - charStart_].mutate(loBound, hiBound);
		}
	}
};

void Chromosome :: recordScore(double score)

{

	sumScores_ += score;

	sumSquares_ += score * score;

	++numTrials_;

};

void Chromosome :: saveChr(ofstream& outFile) {

	int i;
	for (i = 0; i < numDouble_; ++i) outFile << geneDoubleArray_[i].getVal() << ",";
	for (i = 0; i < numBool_; ++i) outFile << geneBoolArray_[i].getVal() << ",";
	for (i = 0; i < numInt_; ++i) outFile << geneIntArray_[i].getVal() << ",";
	for (i = 0; i < numChar_; ++i) outFile << geneCharArray_[i].getVal() << ",";

	outFile << endl;
}

void Chromosome :: saveChr(ofstream& outFile, const vector<long>& nNodes, const vector<int>& inputsToUse) {

	for (int iOut = 0; iOut < Btrue; ++iOut) {
		outFile << inputsToUse[iOut] << ",";
	}
	outFile << endl;
	for (int outLayer = 0; outLayer < 3; ++outLayer)
		outFile << nNodes[outLayer] << ",";
	outFile << endl;
	saveChr(outFile);
	outFile << "//" << endl;
}


void Chromosome :: saveChr(ofstream& outFile, const vector<long>& nNodes, const vector<int>& inputsToUse, const bool bias, const vector<string>& titles, const short OGAgen, const int whichOC) {

	outFile << OGAgen << " " << whichOC << endl;
	for (int iOut = 0; iOut < Btrue; ++iOut) {
		outFile << titles[inputsToUse[iOut]] << ",";
	}
	outFile << endl;
	outFile << ((bias) ? "BIAS" : "NOBIAS") << endl;
	for (int outLayer = 0; outLayer < 3; ++outLayer)
		outFile << nNodes[outLayer] << ",";
	outFile << endl;
	saveChr(outFile);
	outFile << "//" << endl; // Separator
}
