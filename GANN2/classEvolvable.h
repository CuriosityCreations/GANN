// classEvolvable.h: header file for classEvolvable.cpp
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

// Modified from file 'classEvolvable.h' created by Dr. Robert L. Charlebois, NeuroGadgets, Inc.

#include "classChromosome.h"

class Evolvable {
private:
	int numChromosomes_;				// How many type T parameter sets
	int sizeChromosome_;				// The number of type T parameters in each set
	int boolStart_;
	int intStart_;
	int charStart_;
	int Btrue;
	int currentSelection_;				// Which type T parameter set is currently selected
	std::vector<Chromosome*> myChromosomes_;	// Collection of type T parameter sets
	std::vector<Chromosome*> recombinants_;	// Collection of recombinant parameter sets (obviates need for many new/delete)
	Chromosome* intraModifier_;	// If no external modifier Chromosome is specified, intraModifier_ is
													// instantiated (and its parameters set to T(0))
	Chromosome* myModifier_;		// Points either to an external modifier Chromosome, or to *intraModifier_


	void commonInit(Chromosome* extraModifier, const pStr& thePars);
		// Supplies code in common to some of the constructors

public:

	const Chromosome getChr() { return *myChromosomes_[currentSelection_]; };

	Evolvable() { delete intraModifier_; };

	Evolvable(int sizeD, int sizeB, int sizeI, int sizeC, int sizeChr, const double& loBound, const double& hiBound,
		Chromosome* extraModifier, long& setSeed, const pStr& thePars);

	Evolvable(int numChr, int sizeD, int sizeB, int sizeI, int sizeC, const std::vector<double>& range, Chromosome*
		extraModifier, long& setSeed, const pStr& thePars);
		// Constructor that initializes gene data to random value between <loBound> and <hiBound>
	Evolvable(int numChr, int sizeD, int sizeB, int sizeI, int sizeC, std::ifstream& file_in, Chromosome* extraModifier, long& setSeed, const pStr& thePars) throw(Err);
		// Initializes gene data from the file <file_in>.

		// <sizeChr> is input only for verification against the file
		// If the file's numChr is less than <numChr>, the extra chromosomes are replicated from
			// those provided, in a round-robin fashion, and mutated.
			// If the file's numChr is greater, extra are ignored.
		// Calls: Chromosome::Chromosome(), Chromosome::clearGenes(), Chromosome::mutateGene()

	Evolvable(const Evolvable& rhs);
	~Evolvable();

	Evolvable& operator=(const Evolvable& rhs) throw(Err);
		
	int getTrueCount() const { return myChromosomes_[currentSelection_]->getTrueCount(); }

	int itsChromosomeSize() const { return sizeChromosome_; };
	int itsNumChromosomes() const { return numChromosomes_; };
	const Chromosome& theCurrentChromosome() const { return *myChromosomes_[currentSelection_]; };
	int theCurrentChromosomeIndex() const { return currentSelection_; };
	void replaceGenes(const GeneDouble* const & theD, const GeneBool* const & theB, 
					  const GeneInt* const & theI, const GeneChar* const & theC) { myChromosomes_[currentSelection_]->replaceGenes(theD, theB, theI, theC); };

	void itsBoolArray(std::vector<bool>& copyTo) { myChromosomes_[currentSelection_]->itsBoolArray(copyTo); };

//	bool isInput(const int which) const { return (myChromosomes_[currentSelection_]->isInput(which)); };

	void acceptImmigrants(const Evolvable* const& fromPopulation, double fraction, bool bestOnly) throw(Err);
		// The worst <fraction> of chromosomes is replaced by chromosomes from the <fromPopulation>.
		// These are randomly chosen, unless <bestOnly> is true.
		// This function checks that chromosome size is the same, but the user is responsible
		// for ensuring that the chromosomes are homologous.
	const double bestGene(int which) const throw(Err);
		// Returns the unmodified value of gene <which> from the currently most successful chromosome.

		// The chromosomes must be sorted by success first!
	bool changePopulationSize(double factor) throw(Err);
		// Grow or shrink the population. Must sortBySucess() first!
		// To grow, generate new chromosomes as does killClones().
		// To shrink, eliminate the worst chromosomes.

		// <factor> * <numChromosomes_> + 0.5 must be at least 1 (else returns false)
		// Calls: rnd(), Chromosome::Chromosome(), Chromosome::clearScores(), Chromosome::mutateGene()
	void chooseChromosome(int which) throw(Err);
		// Sets the selected chromosome <which> as the currently selected.
		// If <which> is less than zero, the chromosome is chosen randomly.
	void chooseChromosome(int which, double score) throw(Err);
		// Sets the selected chromosome <which> as the currently selected,
		// after updating the previous selection's score.
		// If <which> is less than zero, the chromosome is chosen randomly.
	const Chromosome& emigrateCopy(int which, int arraySize) const throw(Err);
		// Provides a copy of chromosome <which>, or of a random chromosome if <which> is negative.
		// <arraySize> is passed only in order to check that the chromosome sizes are compatible.

	const double gene(int which) const throw(Err);
		// Returns the possibly modified value of gene <which> from the currently selected chromosome.
	const double geneVal(int which) const throw(Err);

		// Returns the unmodified value of gene <which> from the currently selected chromosome.

	int killClones(double minGeneFraction, double minAlleleDifference) throw(Err);
		// This function finds and eliminates identical/similar copies of chromosomes, generating
		// replacements by random recombination and mutation.
		// Chromosomes are too similar if fewer than <minGeneFraction> of their genes
		// are at least <minAlleleDifference> different from one another
		// Returns the number of killed clones.
		// ****Must be way to convert type T into type double****

		// Calls: rnd(), Chromosome::clearScores(), Chromosome::mutateGene()
	void mutate(double fraction, double magnitude, double number, bool upTo) throw(Err);
		// Each of the worst <fraction> chromosomes suffers <number> mutations of up to <magnitude> size.
		// <magnitude> is > 1.0

		// If <upTo> is true, the number of mutations is between one and <number>.
		// Calls: Chromosome::mutateGene()
	void recombine(double parentFraction, double offspringFraction) throw(Err);
		// All genes are unlinked, and a random set is exchanged between random parents.
		// The best <parentFraction> of chromosomes are used to generate the <offspringFraction>,
		// replacing the <offSpringFraction> worst chromosomes.
		// <parentFraction> and <offspringFraction> may overlap.
	bool sortBySuccess(int minimumTrials = 4);
		// Called after the chromosomes are evaluated for fitness.
		// Returns false if there are not enough data available for evaluation.

		// <minimumTrials> cannot be set < 2 (it is forced to be >= 2)
	double itsFitness() const { return myChromosomes_[currentSelection_]->itsFitness(); };
	
	int itsNumGenes() const {return myChromosomes_[currentSelection_]->itsNumGenes();};

	int itsNumDoubleGenes() const {return myChromosomes_[currentSelection_]->itsNumDoubleGenes();};

	void saveChr(ofstream& outFile) {myChromosomes_[currentSelection_]->saveChr(outFile); };

	void saveChr(ofstream& outFile, const vector<long>& num_links, const vector<int>& inputsToUse) { myChromosomes_[currentSelection_]->saveChr(outFile, num_links, inputsToUse); };

	void geneValSet(const int which, const double val) {myChromosomes_[currentSelection_]->geneValSet(which, val); };
};
