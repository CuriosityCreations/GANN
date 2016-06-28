// classSelectableBase defines the basic operations on Selectable objects in the genetic algorithm.
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

// Modified from file 'classSelectableBase.cpp' created by Dr. Robert L. Charlebois, NeuroGadgets, Inc.

#include "classSelectableBase.h"

#define DEBUG_ON false
// The constructor I use
SelectableBase :: SelectableBase(int numSets, int numSuperSets, int numChr, int sizeD, int sizeB, int sizeI, int sizeC,
	const double& loBound, const double& hiBound, long& setSeed, const pStr& thePars)
{
	Btrue = elem(thePars.intGlob,"NUM_INPUTS");
	numParamSets_ = numSets;
	numSuperSets_ = numSuperSets;
	currentSuperSet_ = 0;
	requestedSuperSet_ = -1;
	selectedParamSet_ = 0;
	myModifier_ = new Chromosome(sizeD, sizeB, sizeI, sizeC, thePars);
	myParams_.resize(numSets * numSuperSets);
	for (int i = 0; i < numSets * numSuperSets; ++i) {
		myParams_[i] = new Evolvable(numChr, sizeD, sizeB, sizeI, sizeC, loBound, hiBound, myModifier_, setSeed, thePars);
	};
	myModifier_->clearGenes();
};

SelectableBase :: SelectableBase(int numSets, int numSuperSets, 
	int numChr, int sizeD, int sizeB, int sizeI, int sizeC, const std::vector<double>& range, long& setSeed, const pStr& thePars)
{
	Btrue = elem(thePars.intGlob,"NUM_INPUTS");
	numParamSets_ = numSets;
	numSuperSets_ = numSuperSets;
	currentSuperSet_ = 0;
	requestedSuperSet_ = -1;
	selectedParamSet_ = 0;
	myModifier_ = new Chromosome(sizeD, sizeB, sizeI, sizeC, thePars);
	myParams_.resize(numSets * numSuperSets);
	for (int i = 0; i < numSets * numSuperSets; ++i) {
		myParams_[i] = new Evolvable(numChr, sizeD, sizeB, sizeI, sizeC, range, myModifier_, setSeed, thePars);
	};
	myModifier_->clearGenes();
};

SelectableBase :: SelectableBase(std::ifstream& theFile, int numChr, long& setSeed, const pStr& thePars)
{
	int sizeD, sizeB, sizeI, sizeC;

	theFile >> numParamSets_;
	theFile >> numSuperSets_;
	currentSuperSet_ = 0;

	theFile >> numChr;
	theFile >> sizeD;
	theFile >> sizeB;
	theFile >> sizeI;
	theFile >> sizeC;

	myModifier_ = new Chromosome(sizeD, sizeB, sizeI, sizeC);

	myParams_.resize(numParamSets_ * numSuperSets_);
	for (int i = 0; i < numParamSets_ * numSuperSets_; ++i) {
		myParams_[i] = new Evolvable(numChr, sizeD, sizeB, sizeI, sizeC, theFile, myModifier_, setSeed, thePars);
	};
	requestedSuperSet_ = -1;
	selectedParamSet_ = 0;
	myModifier_->clearGenes();
};

SelectableBase :: ~SelectableBase() {

	for (int i = 0; i < numParamSets_ * numSuperSets_; ++i) {
		delete myParams_[i];
		myParams_[i] = 0;
	}
	delete myModifier_;
	myModifier_ = 0;
}
/* void SelectableBase :: write(Writer& wrtr)
{
	wrtr << numParamSets_;
	wrtr << numSuperSets_;
	wrtr << currentSuperSet_;
	wrtr << myModifier_->itsNumGenes();
	for (int i = 0; i < myParams_.size(); ++i) {
		myParams_[i]->write(wrtr);
	}
}; */

void SelectableBase :: addSuperSet(std::ifstream& theFile, long& setSeed, const pStr& thePars)
{
	int checkNumParams, numChr, sizeChr, sizeD, sizeB, sizeI, sizeC;
	theFile >> checkNumParams;
	theFile >> numChr;
	theFile >> sizeD;
	theFile >> sizeB;
	theFile >> sizeI;
	theFile >> sizeC;
	sizeChr = sizeD + sizeB + sizeI + sizeC;

	#if DEBUG_ON
		if (checkNumParams != numParamSets_)
			throw Err("SelectableBase::addSuperSet(), file has wrong number of parameter sets.");
		if (numChr != (myParams_[0]->itsNumChromosomes())
			throw Err("SelectableBase::addSuperSet(), file has wrong number of chromosomes per set.")
		if (sizeChr != (myParams_[0]->itsChromosomeSize())
			throw Err("SelectableBase::addSuperSet(), file has wrong number of genes per chromosome.");
	#endif

	for (int i = 0; i < numParamSets_; ++i) {
		myParams_[i] = new Evolvable(numChr, sizeD, sizeB, sizeI, sizeC, theFile, myModifier_, setSeed, thePars);
	};
	++numSuperSets_;
};
