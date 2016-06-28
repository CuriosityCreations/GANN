// classEvolvable.cpp defines the Evolvable objects used by the genetic algorithm
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

// Modified from file 'classEvolvable.cpp' created by Dr. Robert L. Charlebois, NeuroGadgets, Inc.

#include "classEvolvable.h"

#define DEBUG_ON false

Evolvable :: Evolvable(int numChr, int sizeD, int sizeB, int sizeI, int sizeC, const double& loBound,
	const double& hiBound, Chromosome* extraModifier, long& setSeed, const pStr& thePars)
{
	Btrue = elem(thePars.intGlob,"NUM_INPUTS");
	numChromosomes_ = numChr;
	boolStart_ = sizeD;
	intStart_ = boolStart_ + sizeB;
	charStart_ = intStart_ + sizeI;
	sizeChromosome_ = charStart_ + sizeC;
	currentSelection_ = 0;
	myChromosomes_.resize(numChr);

	for (int i = 0; i < numChr; ++i) {
		myChromosomes_[i] = new Chromosome(sizeD, sizeB, sizeI, sizeC, loBound, hiBound, setSeed, thePars);
	};
	commonInit(extraModifier, thePars);
};

Evolvable :: Evolvable(int numChr, int sizeD, int sizeB, int sizeI, int sizeC, const std::vector<double>& range,
	Chromosome* extraModifier, long& setSeed, const pStr& thePars)
{
	Btrue = elem(thePars.intGlob,"NUM_INPUTS");
	numChromosomes_ = numChr;
	boolStart_ = sizeD;
	intStart_ = boolStart_ + sizeB;
	charStart_ = intStart_ + sizeI;
	sizeChromosome_ = charStart_ + sizeC;
	currentSelection_ = 0;
	myChromosomes_.resize(numChr);

	for (int i = 0; i < numChr; ++i) {
		myChromosomes_[i] = new Chromosome(sizeD, sizeB, sizeI, sizeC, range, setSeed, thePars);
//		cout << "i = " << i << "   " << sizeD << " " << sizeB << " " << sizeI << " " << sizeC << endl;
	};
	commonInit(extraModifier, thePars);
};




Evolvable :: Evolvable(int numChr, int sizeD, int sizeB, int sizeI, int sizeC, std::ifstream& file_in,
	Chromosome* extraModifier, long& setSeed, const pStr& thePars) throw(Err)
{
	int fileNumChr, fileSizeChr, readRange;
	file_in >> fileNumChr;
	file_in >> fileSizeChr;

	#if DEBUG_ON
		if (fileSizeChr != sizeChr) throw Err("Evolvable::Evolvable(), fileSizeChr has unexpected value.");
	#endif

	numChromosomes_ = numChr;
	boolStart_ = sizeD;
	intStart_ = boolStart_ + sizeB;
	charStart_ = intStart_ + sizeI;
	sizeChromosome_ = charStart_ + sizeC;
	currentSelection_ = 0;

	myChromosomes_.resize(numChromosomes_);

	file_in >> readRange;
	if (readRange) {	// create random values from the ranges for each gene about to be read
		std::vector<double> theRange;
		theRange.resize(sizeChromosome_ * 2);	// lo-hi values for each gene
		for (int i = 0; i < sizeChromosome_ * 2; ++i) {
			file_in >> theRange[i];
		};
		for (int i = 0; i < numChromosomes_; ++i) {
			myChromosomes_[i] = new Chromosome(sizeD, sizeB, sizeI, sizeC, theRange, setSeed, thePars);
		}
	} else {	// read actual gene values from the file
		for (int i = 0; i < fileNumChr && i < numChromosomes_; ++i) {
			myChromosomes_[i] = new Chromosome(sizeD, sizeB, sizeI, sizeC, file_in);
		};
		for (int i = fileNumChr; i < numChromosomes_; ++i) {	// In case the file does not specify all.
			myChromosomes_[i] = new Chromosome(*myChromosomes_[i - fileNumChr]); // round-robin
			for (int g = 0; g < sizeChromosome_; ++g) {			// Mutate the replicated genes
				if (maybe()) {	// half the genes are mutated
					myChromosomes_[i]->mutateGene(g, 2.0);	// *= 0.5 to 2.0 
				}
			}
		}
	};
	commonInit(extraModifier, thePars);
};

Evolvable :: Evolvable(const Evolvable& rhs)
{
	numChromosomes_ = rhs.numChromosomes_;
	sizeChromosome_ = rhs.sizeChromosome_;
	currentSelection_ = rhs.currentSelection_;

	for (int i = 0; i < numChromosomes_; ++i) {
		myChromosomes_[i] = new Chromosome(*rhs.myChromosomes_[i]);
	}
//	myChromosomes_ = rhs.myChromosomes_;
	recombinants_ = rhs.recombinants_;
	if (rhs.intraModifier_) {
		intraModifier_ = new Chromosome(*(rhs.intraModifier_));
		myModifier_ = intraModifier_;
	} else {
		intraModifier_ = 0;
		myModifier_ = rhs.myModifier_;
	}
};

Evolvable :: ~Evolvable() {

	for (int i = 0; i < numChromosomes_; ++i) {
		delete myChromosomes_[i];
		myChromosomes_[i] = 0;
		delete recombinants_[i];
		recombinants_[i] = 0;
	}
	delete intraModifier_;
	intraModifier_ = 0;
}

Evolvable& Evolvable :: operator=(const Evolvable& rhs) throw(Err)

{

	#if DEBUG_ON

		if (numChromosomes_ != rhs.numChromosomes_ || sizeChromosome_ != rhs.sizeChromosome_)

			throw Err("Evolvable::operator=(), illegal assignment");

	#endif



	currentSelection_ = rhs.currentSelection_;

	for (int i = 0; i < numChromosomes_; ++i) {
		myChromosomes_[i] = rhs.myChromosomes_[i];
		recombinants_[i] = rhs.recombinants_[i];
	}

	if (rhs.intraModifier_) {

		*intraModifier_ = *rhs.intraModifier_;

		myModifier_ = intraModifier_;

	} else {

		intraModifier_ = 0;

		myModifier_ = rhs.myModifier_;

	};

	return *this;

};



void Evolvable :: commonInit(Chromosome* extraModifier, const pStr& thePars)
{
	
	int sizeD = boolStart_;
	int sizeB = intStart_ - boolStart_;
	int sizeI = charStart_ - intStart_;
	int sizeC = sizeChromosome_ - charStart_;

	recombinants_.resize(numChromosomes_);
	for (int i = 0; i < numChromosomes_; ++i) {
		recombinants_[i] = new Chromosome(sizeD, sizeB, sizeI, sizeC, thePars);
	};
	if (extraModifier) {
		intraModifier_ = 0;
		myModifier_ = extraModifier;
	} else {
		intraModifier_ = new Chromosome(sizeD, sizeB, sizeI, sizeC, thePars);
		intraModifier_->clearGenes();	// Zeroes the modification values
		myModifier_ = intraModifier_;
	}
};

/* void Evolvable :: write(Writer& wrtr)
{
	wrtr << numChromosomes_;
	wrtr << sizeChromosome_;
	wrtr << currentSelection_;
	for (int i = 0; i < numChromosomes_; ++i) {
		myChromosomes_[i]->write(wrtr);
	}
}; */

void Evolvable :: acceptImmigrants(const Evolvable* const & fromPopulation, 
	double fraction, bool bestOnly) throw(Err)
{
	#if DEBUG_ON

		if (fraction > 1.0) throw Err("Evolvable::acceptImmigrants(), fraction invalid.");

	#endif



	int numImmigrants = int(double(numChromosomes_) * fraction + 0.5);

	int count = -1;
	if (bestOnly) {
		for (int i = numChromosomes_ - numImmigrants; i < numChromosomes_; ++i) {
			*myChromosomes_[i] = fromPopulation->emigrateCopy(++count, sizeChromosome_);
			myChromosomes_[i]->clearScores();
		}
	} else {
		for (int i = numChromosomes_ - numImmigrants; i < numChromosomes_; ++i) {
			*myChromosomes_[i] = fromPopulation->emigrateCopy(count, sizeChromosome_);
			myChromosomes_[i]->clearScores();
		}
	}
};

const double Evolvable :: bestGene(int which) const throw(Err)
{

	#if DEBUG_ON
		if (which < 0 || which >= sizeChromosome_) throw Err("Evolvable::bestGene(), index out of range.");

	#endif

	double bG = myChromosomes_[0]->getVal(which);

	return bG;
};

bool Evolvable :: changePopulationSize(double factor) throw(Err)
{
	#if DEBUG_ON

		if (factor <= 0.0) throw Err("Evolvable<T>::changePopulationSize(), factor must be > 0.0");

	#endif



	int newPopSize = int(double(numChromosomes_) * factor + 0.5);

	if (newPopSize == numChromosomes_) return true;

	if (newPopSize == 0) return false;

	if (newPopSize < numChromosomes_) {
		for (int i = newPopSize; i < numChromosomes_; ++i) {
			myChromosomes_.pop_back();

			recombinants_.pop_back();
		}
	} else {
		// Create the new chromosomes as mutated quadriparental recombinants:
		int parent[4];
		Chromosome* toPush = new Chromosome(boolStart_, intStart_, charStart_, sizeChromosome_);
		for (int i = numChromosomes_; i < newPopSize; ++i) {
			recombinants_[i] = new Chromosome(boolStart_, intStart_, charStart_, sizeChromosome_);
			myChromosomes_[i] = toPush;
			for (int p = 0; p < 4; ++p) {
				parent[p] = rnd(numChromosomes_);
			};	// no checking for duplication in parentage
			for (int g = 0; g < sizeChromosome_; ++g) {
				myChromosomes_[i]->geneValSet(g, myChromosomes_[parent[rnd(4)]]->getVal(g));
				if (maybe()) {
					myChromosomes_[i]->mutateGene(g, 1.05);	// some mutation too (*= 0.9523 to 1.05) for good luck
				}
			}
		delete toPush;
		}
	};
	numChromosomes_ = newPopSize;

	return true;
};

void Evolvable :: chooseChromosome(int which) throw(Err) 
{

	#if DEBUG_ON
		if (which >= numChromosomes_) throw Err("Evolvable::chooseChromosome(), index out of range.");

	#endif


	if (which < 0) {
		currentSelection_ = rnd(numChromosomes_);
	} else {
		currentSelection_ = which;
	}
};

void Evolvable :: chooseChromosome(int which, double score) throw(Err)
{

	#if DEBUG_ON
		if (which >= numChromosomes_) throw Err("Evolvable::chooseChromosome(), index out of range.");

	#endif


	myChromosomes_[currentSelection_]->recordScore(score);
	if (which < 0) {
		currentSelection_ = rnd(numChromosomes_);
	} else {
		currentSelection_ = which;
	}
};

const Chromosome& Evolvable :: emigrateCopy(int which, int arraySize) const throw(Err)
{

	#if DEBUG_ON
		if (arraySize != sizeChromosome_) throw Err("Evolvable::emigrateCopy(), mismatched chromosome sizes.");
		if (which >= numChromosomes_) throw Err("Evolvable::emigrateCopy(), invalid index.");

	#endif


	if (which < 0) which = rnd(numChromosomes_);
	return *myChromosomes_[which];
};

const double Evolvable :: gene(int which) const throw(Err)
{
	#if DEBUG_ON

		if (which < 0 || which >= sizeChromosome_) throw Err("Evolvable::gene(), index out of range.");

	#endif



	return myChromosomes_[currentSelection_]->getVal(which) + myModifier_->getVal(which);
};

const double Evolvable :: geneVal(int which) const throw(Err)

{

	#if DEBUG_ON

		if (which < 0 || which >= sizeChromosome_) throw Err("Evolvable::geneVal(), index out of range.");

	#endif



	return myChromosomes_[currentSelection_]->getVal(which);

};



int Evolvable :: killClones(double minGeneFraction, double minAlleleDifference) throw(Err)
{
	int numDifferent = int(double(numChromosomes_) * minGeneFraction + 0.5);



	#if DEBUG_ON
		if (numDifferent < 1 || numDifferent > numChromosomes_)
			throw Err("Evolvable::killClones(), invalid minGeneFraction argument.");

	#endif


	double one, two;

	int count;
	int numReplaced = 0;
	int parent[4];
	bool different;
	minAlleleDifference += 1.0;	// saves computations within k-loop below
	for (int i = 0; i < numChromosomes_ - 1; ++i) {
		for (int j = i + 1; j < numChromosomes_; ++j) {
			count = 0;
			for (int k = 0; k < sizeChromosome_; ++k) {
				one = double(myChromosomes_[i]->getVal(k));
				two = double(myChromosomes_[j]->getVal(k));
				different = false;	// The following test is inelegant, but...
				if (one * two < 0.0) { 
					different = true;
				} else {	// same sign
					one = fabs(one);
					two = fabs(two);
					if (one < two) {
						if (one * minAlleleDifference < two) different = true;
					} else {
						if (two * minAlleleDifference < one) different = true;
					}
				};
				if (different) {
					if (++count >= numDifferent) break;
				}
			};
			if (count < numDifferent) {
				// Need a new chromosome j
				for (int p = 0; p < 4; ++p) {	// First, create a quadriparental recombinant
					parent[p] = rnd(numChromosomes_);
				};	// no checking for duplication in parentage
				for (int g = 0; g < sizeChromosome_; ++g) {
					myChromosomes_[j]->geneValSet(g, myChromosomes_[parent[rnd(4)]]->getVal(g));
					if (maybe()) {

					myChromosomes_[j]->mutateGene(g, 1.05);	// some mutation too (*= 0.9523 to 1.05) for good luck
					}
				};
				myChromosomes_[j]->clearScores();
				++numReplaced;
			}
		}
	};
	return numReplaced;
};

void Evolvable :: mutate(double fraction, double magnitude, double number, bool upTo) throw(Err)
{
	int num = int (number * double(sizeChromosome_));
	int numMutants = int(double(numChromosomes_) * fraction + 0.5);

	#if DEBUG_ON
		if (numMutants < 1 || numMutants > numChromosomes_) {
			throw Err("Evolvable::mutate(), fraction invalid.");

		};

		if (magnitude <= 1.0) throw Err("Evolvable::mutate(), magnitude invalid.");

	#endif


	if (upTo) {
		int howMany;
		for (int i = numChromosomes_ - numMutants; i < numChromosomes_; ++i) {
			howMany = rnd(num) + 1;
			for (int j = 0; j < howMany; ++j) {
				myChromosomes_[i]->mutateGene(rnd(sizeChromosome_), magnitude);
			};
			myChromosomes_[i]->clearScores();
		}
	} else {
		for (int i = numChromosomes_ - numMutants; i < numChromosomes_; ++i) {
			for (int j = 0; j < num; ++j) {
				myChromosomes_[i]->mutateGene(rnd(sizeChromosome_), magnitude);
			};
			myChromosomes_[i]->clearScores();
		}
	}
};

void Evolvable :: recombine(double parentFraction, double offspringFraction) throw(Err)
{
	int numParents = int(double(numChromosomes_) * parentFraction + 0.5);
	int numOffspring = int(double(numChromosomes_) * offspringFraction + 0.5);

	#if DEBUG_ON
		if ((numParents < 2) || (numParents > numChromosomes_) || (numOffspring < 1) ||
		  (numOffspring > numChromosomes_)) throw Err("Evolvable::recombine(), fraction invalid.");

	#endif

	int parent1, parent2;
	double proportionFrom1;
	for (int i = 0; i < numOffspring; ++i) {
		recombinants_[i]->clearGenes();
		parent1 = rnd(numParents);	// potential parents have index 0 to numParents-1
		parent2 = rnd(numParents);
		while (parent1 == parent2) {
			parent2 = rnd(numParents);
		};
		proportionFrom1 = rnd(0.9) + 0.05; // i.e., from 5% to 95%: arbitrary (for GeneDouble)
		for (int j = 0; j < boolStart_; ++j) {
			if (proportionFrom1 > rnd(1.0)) {
				recombinants_[i]->geneValSet(j, myChromosomes_[parent1]->getVal(j));
			} else {
				recombinants_[i]->geneValSet(j, myChromosomes_[parent2]->getVal(j));
			}
		};

		int oneTrue = 0;
		int twoTrue = 0;
		int currentBool = 0;
		int theVal = 0;
		int numBool = intStart_ - boolStart_;
		if (numBool != 0) {
			while (oneTrue < Btrue / 2) { // Collect 1/2 (rounded down) of inputs from parent 1
				currentBool = rnd(numBool) + boolStart_;
				theVal = (int) myChromosomes_[parent1]->getVal(currentBool);
				if ((theVal == 1) && recombinants_[i]->getVal(currentBool) == 0) {
					recombinants_[i]->geneValSet(currentBool, theVal);
					++oneTrue;
				}
			}
			while (twoTrue < Btrue - oneTrue) { // Collect 1/2 (rounded up) of inputs from parent 2
				currentBool = rnd(numBool) + boolStart_;
				theVal = (int) myChromosomes_[parent2]->getVal(currentBool);
				if ((theVal == 1) && recombinants_[i]->getVal(currentBool) == 0) {
					recombinants_[i]->geneValSet(currentBool, theVal);
					++twoTrue;
				}
			}
		}
		recombinants_[i]->clearScores();
	};
	for (int k = 0; k < numOffspring; ++k) {	// Separate loop in order to permit overlapping ranges
		for (int l = 0; l < sizeChromosome_; ++l) {
			myChromosomes_[numChromosomes_ - 1 - k]->geneValSet(l, recombinants_[k]->getVal(l));
		}
	}
};

bool Evolvable :: sortBySuccess(int minimumTrials)
{
	
	// Success is defined by the mean score minus the score's standard deviation
	// First, must check that numTrials_[i] is >= <minimumTrials> for all i:
	if (minimumTrials < 2) minimumTrials = 2;
	for (int i = 0; i < numChromosomes_; ++i) {
		if (myChromosomes_[i]->itsNumTrials() < minimumTrials) {
			return false;
		}
	};

	for (int i = 0; i < numChromosomes_; ++i) {
		myChromosomes_[i]->computeFitness();
	};

	// Sort (in descending order of success):
	std::sort(myChromosomes_.begin(), myChromosomes_.end(), CmpFitnesses());

	return true;
};
