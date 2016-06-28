// classChromosome.h: header file for clasChromosome.cpp
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

// Modified from file 'classChromosome.h' created by Dr. Robert L. Charlebois, NeuroGadgets, Inc.

#include "classGene.h"

class Chromosome {
private:

	double fitness_;
	double sumScores_;	// used in evaluating fitness
	double sumSquares_;	// used in evaluating fitness
	int numTrials_;		// used in evaluating fitness
	int numGenes_;

	int boolStart_; // The first bool gene
	int intStart_;
	int charStart_;

	int numDouble_;
	int numBool_;
	int numInt_;
	int numChar_;
	
	int Btrue; // The number of true bool genes
	
	GeneDouble* geneDoubleArray_;
	GeneBool* geneBoolArray_;
	GeneInt* geneIntArray_;
	GeneChar* geneCharArray_;

	void commonInit(const int sizeD, const int sizeB, const int sizeI, const int sizeC);

		// Code in common to most of the constructors

	void mutateGene(int which, double howMuch) throw(Err);
		// Called from within Evolvable<T> functions,
		// modifying gene <which> by *= 1.0/howMuch to howMuch

		// e.g., if <howMuch> == 4.0, gene is multiplied by 0.25 to 4.0
		// <howMuch> must be > 1.0
		// Multiplying type <T> by a double must make sense here...

	void mutateGene(int which, const double& loBound, const double& hiBound) throw(Err);

		// More general version of above mutateGene()

		// Generates a random value between <loBound> and <hiBound>

		// Here, rnd(T) must make sense.
public:
	Chromosome() { };

	explicit Chromosome(const int sizeD, const int sizeB, const int sizeI, const int sizeC);

		// chromosome, initialized to zero 

	explicit Chromosome(const int sizeD, const int sizeB, const int sizeI, const int sizeC, const pStr& thePars);

		// chromosome, initialized to zero 

	Chromosome(const int sizeD, const int sizeB, const int sizeI, const int sizeC, const double& loBound, const
		double& hiBound, long& setSeed, const pStr& thePars);

		// chromosome, initialized to random value between <loBound> and <hiBound>

	Chromosome(const int sizeD, const int sizeB, const int sizeI, const int sizeC, const std::vector<double>& range, 
		long& setSeed, const pStr& thePars);

		// chromosome, initialized to random values defined in <theRange> for each gene

	Chromosome(const int sizeD, const int sizeB, const int sizeI, const int sizeC, std::ifstream& inFile) throw(FileErr);
		// chromosome, initial values from <inFile>
		// <size> is only used to check against file.

	Chromosome(const int num_double, const double * const the_links);

	Chromosome(const Chromosome& rhs);

	~Chromosome();
	Chromosome& operator=(const Chromosome& rhs) throw(Err);
//	void write(Writer& wrtr);
	
	void clearGenes();

	void itsBoolArray(std::vector<bool>& copyTo) 
	{	int bink = 0;
		for (int i = 0; i < numBool_; ++i) {
			copyTo[i] = (geneBoolArray_[i].getVal() == 1.0);
			if (geneBoolArray_[i].getVal() == 1.0) ++bink;
		}
	};

	void whichBoolTrue() {
		for (int i = 0; i < numBool_; ++i) {
			if (geneBoolArray_[i].getVal() == 1.0) {
				cout << i << ",";
			}
		}
		cout << endl;
	}

	double getVal(const int which) const { 
		if (which < boolStart_) {
			return geneDoubleArray_[which].getVal(); 
		} if (which < intStart_) {
			return geneBoolArray_[which - boolStart_].getVal();
		} if (which < charStart_) {
			return geneIntArray_[which - intStart_].getVal();
		} return geneCharArray_[which - charStart_].getVal();
	};

	int getTrueCount() const {  int n = 0; int i;
								for (i = 0; i < numBool_; ++i)
									if (geneBoolArray_[i].getVal() == 1.0) ++n; 
								return n;
							 };

	double itsFitness() const { return fitness_; };
	int itsNumGenes() const { return numGenes_; };

	int itsNumDoubleGenes() const { return numDouble_; };

	void geneValSet(const int which, const double val) { 
		if (which < boolStart_) {
			geneDoubleArray_[which].forceVal(val); 
		} else {
			if (which < intStart_) {
				geneBoolArray_[which - boolStart_].forceVal(val);
			} else {
				if (which < charStart_) {
					geneIntArray_[which - intStart_].forceVal(val);
				} else geneCharArray_[which - charStart_].forceVal(val);
			}
		}
	};


	int itsNumTrials() const { return numTrials_; };

	void modifyGeneVal(int which, const double& howMuch) { 
		if (which < boolStart_) {
			geneDoubleArray_[which].forceVal(geneDoubleArray_[which].getVal() + howMuch); 
		} else {
			if (which < intStart_) {
				geneBoolArray_[which - boolStart_].forceVal(geneBoolArray_[which - boolStart_].getVal() + howMuch);
			} else {
				if (which < charStart_) {
					geneIntArray_[which - intStart_].forceVal(geneIntArray_[which - intStart_].getVal() + howMuch);
				} else geneCharArray_[which - charStart_].forceVal(geneCharArray_[which - charStart_].getVal() + howMuch);
			}
		}
	};
	// Note: Permits a sign change of the value of <gene>

	void replaceGenes(const GeneDouble* const & theD, const GeneBool* const & theB, 
					  const GeneInt* const & theI, const GeneChar* const & theC);

	void clearScores();
		// Zeroes the scoring information

	void computeFitness() throw(Err);

		// fitness is the mean score minus the score's standard deviation
	void recordScore(double score);
	
	void saveChr(ofstream& outFile);
	void saveChr(ofstream& outFile, const vector<long>& nNodes, const vector<int>& inputsToUse);
	void saveChr(ofstream& outFile, const vector<long>& nNodes, const vector<int>& inputsToUse, const bool bias, const vector<string>& titles, const short OGAgen, const int whichOC);

	friend class Evolvable;

};

class CmpFitnesses {

public:

	int operator() (Chromosome* chr1,  Chromosome* chr2)

		{ return chr1->itsFitness() > chr2->itsFitness(); };

};
