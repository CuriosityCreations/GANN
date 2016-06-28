// classSelectableBase.h: header file for classSelectableBase.cpp
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

// Modified from file 'classSelectableBase.h' created by Dr. Robert L. Charlebois, NeuroGadgets, Inc.

#include "classEvolvable.h"

class SelectableBase {

private:

	SelectableBase(const SelectableBase& rhs);				// block copy construction

	SelectableBase& operator=(const SelectableBase& rhs);	// block assignment

protected:

	int numParamSets_;			// How many sets of parameters

	int numSuperSets_;			// How many sets of sets of parameters

	int currentSuperSet_;		// Which set of sets of parameters is currently active

	int requestedSuperSet_;		// Next SuperSet to try out; if < 0, keep currentSuperSet_

	int selectedParamSet_;		// Which set of parameters is currently active

	int Btrue;

	Chromosome* myModifier_;	// Modifies values of parameters in the Evolvable sets

	std::vector<Evolvable*> myParams_;		// Sets of parameters

public:

	SelectableBase(int numSets, int numSuperSets, int numChr, int sizeD, int sizeB, int sizeI, int sizeC, const
		double& loBound, const double& hiBound, long& setSeed, const pStr& thePars);

		// constructor that creates random chromosomes, whose values are linearly distributed

		// between <loBound> and <hiBound>

	SelectableBase(int numSets, int numSuperSets, int numChr, int sizeD, int sizeB, int sizeI, int sizeC, const
		std::vector<double>& range, long& setSeed, const pStr& thePars);

	SelectableBase(std::ifstream& theFile, int numChr, long& setSeed, const pStr& thePars);

		// constructor that initializes from a text file; <numChr> is ignored (distinguishes constructors)

		// chromosomes based on each gene having a random value within a range read from <theFile>

		// Calls: Chromosome::Chromosome(), Evolvable::Evolvable(), SelectableBase::commonConstructorCode()

			// Chromosome::clearGenes()

	~SelectableBase();

	int getTrueCount() const { return myParams_[selectedParamSet_]->getTrueCount(); }

	const double gene(int which) const { return myParams_[selectedParamSet_]->gene(which); };

	const double geneVal(int which) const { return myParams_[selectedParamSet_]->geneVal(which); };

	void geneValSet(int which, double toSet) { myParams_[selectedParamSet_]->geneValSet(which, toSet); };

	void itsBoolArray(std::vector<bool>& copyTo) { myParams_[selectedParamSet_]->itsBoolArray(copyTo); };

	double itsFitness() { return myParams_[selectedParamSet_]->itsFitness(); }

	const Chromosome& itsChromosome() const { return myParams_[selectedParamSet_]->theCurrentChromosome(); };

	Chromosome* itsModifierChrPointer() const { return myModifier_; };

	void replaceGenes(const GeneDouble* const & theD, const GeneBool* const & theB, 
					  const GeneInt* const & theI, const GeneChar* const & theC) { myParams_[selectedParamSet_]->replaceGenes(theD, theB, theI, theC); };

	void setNewSuperSet(int which) { requestedSuperSet_ = which; };

	void setParamSet(int the) { selectedParamSet_ = the; }

	void addSuperSet(std::ifstream& theFile, long& setSeed, const pStr& thePars);

		// Adds a SuperSet, described in <theFile>

	void chooseChromosome(int which) { myParams_[selectedParamSet_]->chooseChromosome(which); }
		
	void chooseChromosome(int which, double score) { myParams_[selectedParamSet_]->chooseChromosome(which, score); }

	virtual void optimizeParameters(const std::vector<double> GARange) = 0;

	void saveChr(ofstream& outFile) { myParams_[selectedParamSet_]->saveChr(outFile); }

	int whichChrSelected() { return myParams_[selectedParamSet_]->theCurrentChromosomeIndex(); };
};
