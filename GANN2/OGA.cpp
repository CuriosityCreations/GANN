// OGA.cpp defines the usage of the Outer Genetic Algorithm object.
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

#include"OGA.h"

// CLASS OGA: The Outer Genetic Algorithm --------------------------------------------------------

OGA :: OGA(const int numSets, const int numSuperSets, const int numChr, const int sizeD, 
			const int sizeB, const int sizeI, const int sizeC, const vector<double>& range, 
			const bool isGA, const bool IVO, const string inFileName, long& setSeed, const pStr& thePars)
			: SelectableBase (numSets, numSuperSets, numChr, sizeD, sizeB, sizeI, sizeC, range, setSeed, thePars) {
		       
	numInputs = sizeB;

	// Read inputs from file
	setInputs(inFileName);

	string netName = elem(thePars.stringGlob,"NET_DEF");
	saveNetFile.open(netName.c_str());

	totalOC = elem(thePars.intGlob,"GA_CHR") * elem(thePars.intGlob,"GA_EVO");
	OGArounds = elem(thePars.intGlob,"OGA_TRAIN_ROUNDS");

	// Evolutionary parameters
	recRate = elem (thePars.doubleGlob,"OGA_REC_RATE");
	recRepl = elem (thePars.doubleGlob,"OGA_REC_REPL");
	mutRate = elem (thePars.doubleGlob,"OGA_MUT_RATE");
	mutAmt = elem (thePars.doubleGlob,"OGA_MUT_AMT");
	mutProp = elem (thePars.doubleGlob,"OGA_MUT_PROP");
	killFrac = elem (thePars.doubleGlob,"OGA_KILL_FRAC");
	killDiff = elem (thePars.doubleGlob,"OGA_KILL_DIFF");

	const int numInnerRounds = elem (thePars.intGlob,"NN_TRAIN_RUNS");

	// Create oCscore object
	if (!isGA) {
		chrScores = new BPscore(numInnerRounds);
		BPrepl = elem(thePars.intGlob,"REPLICATES");
		GAnumRec = 1;
		isGANN = false;
	} else {
		chrScores = new GAscore(numInnerRounds);
		BPrepl = 1;
		GAnumRec = elem(thePars.intGlob,"NUM_TO_RECORD");
		isGANN = true;
	}

	//input vector optimization?
	if (IVO) {
		vecOp = true;
	} else vecOp = false;
	
	// Build the exp lookup table
	for (int eBuild = 0; eBuild < ETABLESIZE * 2; ++eBuild) {
		eTable[eBuild] = 2.0/(1.0 + exp(-2 * (ETABLESIZE - eBuild) / 1000.0)) - 1;  // Equivalent to negative (eBuild - 10000.0)
	}
	
	bestTime = 1000000.0;
	worstTime = 0.0;

	best.resize(howManySets);
	bestNet.resize(howManySets);
	for (int i = 0; i < howManySets; ++i) {
		best[i] = 0;
		bestNet[i].resize(5);
	}
}

OGA :: ~OGA() {

	delete chrScores; 
	chrScores = 0;
	
	saveNetFile.close();

}

// setInputs reads the NN input values from the Perl-generated file, determines the sizes of all
// sets, and stores the input values in the OGA
void OGA :: setInputs (const string inFileName) {

	string header;
	string titleLine;
	string currentLine;

	int totalNum = 0;

	ifstream inFile(inFileName.c_str());
	getline(inFile, header);
	int numVars = atoi(header.c_str()); // Since the first line should be a single number

	int i;
	// Create the set storage
	getline(inFile, header);
	vector<string> catSizes = split(header, ",");
	howManySets = catSizes.size();
	howManyCategories = howManySets / 2;
	sets.resize(howManySets);
	for (i = 0; i < howManySets; ++i) {
		setSizes.push_back(atoi((catSizes[i]).c_str()));
		(sets[i]).resize(setSizes[i]);
		totalNum += setSizes[i];
		for (int j = 0; j < setSizes[i]; ++j) {
			(sets[i][j]).resize(numVars);
		}
	}

	// Store the titles for each category
	getline(inFile, titleLine);

	titles = split(titleLine, ",", 4); // Only use the fourth+ elements

	vector<int> storeCounter(howManySets, int(0));
	int whichSet = 0;

	// Store the input values in the appropriate arrays
	for (int createSet = 0; createSet < totalNum; ++createSet) {
		//if (createSet % 10 == 0) {
		//	cout << createSet << ",";
		//	if (createSet % 500 == 0) {
		//		cout << endl;
		//	}
		//}
		
		vector<string> chopLine; // Space for posNeg, trainTest, seqID
		getline(inFile, currentLine);
		chopLine = split(currentLine, ",");

		whichSet = atoi((chopLine[0]).c_str()) + howManyCategories * atoi((chopLine[1]).c_str());

		for (int sInputs = 3; sInputs < numVars + 3; ++sInputs) { // Start after PosNeg, TrainTest, SeqID
			sets[whichSet][storeCounter[whichSet]][sInputs - 3] = atof ((chopLine[sInputs]).c_str());
		}
		++storeCounter[whichSet];
	}
	inFile.close();
}

// Run an entire OGA training round
void OGA :: execute (const short OGAgen, const std::vector<double> GARange, const int numParams, const int cNum, const pStr& thePars) {
	
	int selRoundGA, evoRoundGA, chrRoundGA;

	chrScores->paramInit(); // Initialize next round

	selRoundGA = 0;

	while (selRoundGA < elem(thePars.intGlob, "GA_SEL")) {
		evoRoundGA = 0;		
		while (evoRoundGA < elem(thePars.intGlob, "GA_EVO")) {
			setParamSet(evoRoundGA);
			chrRoundGA = 0;
			while (chrRoundGA < cNum) {
				// Choose the next set of GA parameters
				chooseChromosome(chrRoundGA);

				// Get the relevant neural network parameters from the Chromosome's genes
				vector<double> currentGAPar;
				currentGAPar.resize(numParams);
				for (int copyParams = 0; copyParams < numParams; ++copyParams) {
					currentGAPar[copyParams] = gene(copyParams);
				}
				
				// Recreate the object to score the neural network with the current Chr and Evo size
				chrScores->nextConfig(currentGAPar, chrRoundGA, OGAgen, elem(thePars.intGlob,"NUM_INPUTS"));

				OCrun(OGAgen, currentGAPar, numParams, thePars); // Train the outer Chromosome

				chrRoundGA++; // Next outer Chromosome
			}
			evoRoundGA++; // Next outer Evolvable
		}
		selRoundGA++; // Next outer Selectable
	}	

	optimizeParameters(GARange); // Evolve the OGA

} // Execute

// Run an outer Chromosome of the OGA
void OGA :: OCrun (const short OGAgen, std::vector<double>& currentGAPar, const int numParams, const pStr& thePars) {
					
	bool saveThisNet = false;
	
	// Create the ANN and set the input vectors
	vector<int> inputsToUse(elem (thePars.intGlob, "NUM_INPUTS"));
	createNN(currentGAPar, numParams, OGAgen, inputsToUse, thePars);
	
	// How many of the best scores to record
	const int scoresToRecord = (isGANN) ? GAnumRec : BPrepl;

	// Record the best generalization scores
	vector<double> bestGen;
	bestGen.resize(scoresToRecord);
	for (int qw = 0; qw < scoresToRecord; ++qw) bestGen[qw] = 0.0;

	long seed = 0 - rand(); // Seed to randomize neural networks

	// Timers (not for selection purposes)
	clock_t beginnn, enddd;
	beginnn = clock(); 
	// Measure the time taken per run
	clock_t startEl, stopEl;
	long elapsedTime = 0;
		
	if (!isGANN) { chrScores->chrScoresInit(); }

	int replNum = 0;
	const int repStop = BPrepl;

	for (replNum = 0; replNum < repStop; ++replNum) {
	
		// Randomize weights for DNA_net
		neuralNet->Reset(seed);
		chrScores->reset();

		cout << "Replicate " << replNum << "  ";
	
		bool willStop = false;
		int generation = 0;

		if (isGANN) { chrScores->chrScoresInit(); }

		// Run the network training phase, with phase termination dependent on the WHILE loop
		while (!willStop) {

			startEl = clock();
			neuralNet->NNexec(generation, false, chrScores, eTable, setSizes);

			// Calculate elapsed and recorded time
			stopEl = clock();

			elapsedTime = (long) (stopEl - startEl) / CLOCKS_PER_SEC;		

			// If this is a testing round, test
			neuralNet->NNexec(generation, true, chrScores, eTable, setSizes);

			// Check to see whether best scores have been exceeded
			if (chrScores->recordScores(elem(thePars.doubleGlob,"MIN_GEN"))) { // If TRUE, then new best gen score has been reached
				neuralNet->copyTo(toSave);
				saveThisNet = true;
			}
			++generation;
			willStop = neuralNet->stopConditions(generation, chrScores);
			
			chrScores->ifNewBestTime(elapsedTime);
			if (isGANN) chrScores->chrScoresOutput(false);
			chrScores->nextInner();

		} // Next generation				

		// BestGen = the best test FP score
		for (int i = 0; i < GAnumRec; ++i) {
			bestGen[replNum + i] = chrScores->getBestScores(2, i); // Get the ith best score (0 to NUM_TO_RECORD)
		}

		cout << bestGen[replNum];

		if (!isGANN) chrScores->chrScoresOutput(false);

		cout << endl;

	} // Next replicate

	if (saveThisNet) {
		vector<long> nNodes(3);
		nNodes[0] = elem(thePars.intGlob,"NUM_INPUTS");
		nNodes[1] = (int) currentGAPar[XHIDNODE];
		nNodes[2] = (int) currentGAPar[XOUTNODE];
		toSave->saveChr(saveNetFile, nNodes, inputsToUse, neuralNet->isBiased(), titles, OGAgen, whichChrSelected());
	}

	enddd = clock();

	std::sort(bestGen.begin(), bestGen.end());
	cout << "   Best generalization score: " << bestGen[scoresToRecord - 1] << endl;
	cout << "   Time for one outer Chromosome training run: " << (double)(enddd - beginnn) / CLOCKS_PER_SEC << endl;	

	const double worstScore = elem(thePars.doubleGlob,"WORST_SCORE");
	// Correct the scores and record them
	for (int reap = 0; reap < scoresToRecord; ++reap) {
		bestGen[reap] -= worstScore; // Subtract the base score
		bestGen[reap] /= (1 - worstScore); // Force range between 0 and 1
		if (bestGen[reap] < 0.0) bestGen[reap] = 0.0; // If less than 0, make it 0
		chooseChromosome(whichChrSelected(), bestGen[reap]); // Add the average scores to the appropriate oC
	}

	for (int bestS = 0; bestS < howManyCategories; ++bestS) {
		if (chrScores->getBestScores(bestS, 0) > best[bestS]) {
			best[bestS] = chrScores->getBestScores(bestS, 0);
		}
	}

	chrScores->endOfRound(currentGAPar, inputsToUse);

	delete toSave;
	toSave = 0;
	delete neuralNet;
	neuralNet = 0;

} // OCrun

void OGA :: createNN (std::vector<double>& currentGAPar, const int numParams, const int OGAgen, vector<int>& inputsToUse, const pStr& thePars) {

	bool netIsBias = (currentGAPar[XISBIAS] >= 1.0);

	const int currentChr = whichChrSelected();
	const int inNodes = inputsToUse.size();
	const int hidNodes = (int) currentGAPar[XHIDNODE];
	const int outNodes = (int) currentGAPar[XOUTNODE];
	const int netLinks = inNodes * hidNodes + hidNodes * outNodes + ((netIsBias) ? hidNodes + outNodes : 0); // The links from the bias nodes

	std::vector<int> nodeSizes;
	nodeSizes.resize(3);
	nodeSizes[0] = inNodes;
	nodeSizes[1] = hidNodes;
	nodeSizes[2] = outNodes;

	// Create the genetic algorithm or backpropagation neural network
	if (isGANN) {
		const int thisNumEvo = (int) currentGAPar[XNUMEVO];
		const int thisNumChr = (int) currentGAPar[XNUMCHR] / thisNumEvo;		
		currentGAPar[XNUMCHR] /= thisNumEvo;
		const int totalChr = thisNumEvo * thisNumChr;
		long setSeed = rand();
		neuralNet = new GANet(currentGAPar, numParams, thisNumEvo, 1, thisNumChr, 
								netLinks, 0, 0, 0, -1.0, 1.0, 3, nodeSizes, setSizes, setSeed, thePars); // New GANN
		toSave = new Chromosome(netLinks,0,0,0);
	} else {
		neuralNet = new BPNet(currentGAPar, numParams, 3, nodeSizes, setSizes, thePars); // New BPNN
		toSave = new Chromosome(netLinks,0,0,0);
	}

	// Build an array to indicate which inputs will be used (all if no input vector optimization)
	std::vector<bool> inputArray;
	inputArray.resize(numInputs);

	itsBoolArray(inputArray);	

	int tally = 0;
	for (int makeArray = 0; makeArray < numInputs; ++makeArray) {
		if (inputArray[makeArray]) {
			inputsToUse[tally] = makeArray;
			++tally;
		}
	}

	// Get train and test input vectors
	neuralNet->buildFragments(inputsToUse, sets); 

	cout << "\nGA Parameter Optimization Round: " << OGAgen + 1 << " of " << OGArounds <<
		        "   outer Chromosome #" << chrScores->getChrID() + 1 << " of " << totalOC << endl;
	cout << "Input nodes: " << inNodes << "  Hidden Nodes: " << hidNodes << "  Output nodes: " << outNodes << endl;

}

// Evolve the OGA
void OGA :: optimizeParameters (const std::vector<double> GARange) {

		int r;
		for (int i = 0; i < numParamSets_; ++i) {
			if (myParams_[i]->sortBySuccess()) {
				myParams_[i]->recombine(recRate, recRepl);	// Might consider optimizing the rate of evolution as well
				myParams_[i]->mutate(mutRate, mutAmt, mutProp, false);	// though this might be dangerous
				if (numParamSets_ > 1 && rnd(1.0) < 0.05) {
					r = rnd(numParamSets_);
					while (r == i) r = rnd(numParamSets_);
					myParams_[i]->acceptImmigrants(myParams_[r + currentSuperSet_ * numParamSets_], 0.05, maybe());
						// The immigrants replace some of the mutants
				};
//				myParams_[i]->killClones(killProp, killDiff);
			}
		};
	int numOfGenes = myParams_[0]->itsNumGenes();

	// Ensure that GA params are within the specified boundaries
	for (int checkParams = 0; checkParams < myParams_[0]->itsNumChromosomes(); ++checkParams) {
		chooseChromosome(checkParams);
		for (int checkGenes = 0; checkGenes < myParams_[0]->itsNumDoubleGenes(); ++checkGenes) {
			if (myParams_[0]->geneVal(checkGenes) < GARange[checkGenes * 2]) myParams_[0]->geneValSet(checkGenes, GARange[checkGenes * 2]);
			else if (myParams_[0]->geneVal(checkGenes) > GARange[checkGenes * 2 + 1]) myParams_[0]->geneValSet(checkGenes, GARange[checkGenes * 2 + 1]);
		}
	}
}

// Write information from (stuff). Don't know if I'll use this
// void OGA :: write(Writer& wrtr) {}

// Output the start time and headers for master file
void OGA :: startRun(ofstream& outFile) {

	struct tm *newtime;
	time_t aclock;
	
   	time( &aclock );                
    newtime = localtime( &aclock ); 
    
	outFile << "Optimizing Genetic Algorithm Neural Networks using " << ((isGANN) ? "Genetic Algorithm" : "Backpropagation") << endl;
	outFile << "Run started " << asctime(newtime) << endl;
	outFile << "---------------------------------------\n" << endl;

	outFile << "Generation, bestTestFPscore, bestDistTestScore, " << endl;
} // OGA :: startRun

// Output the generation parameters
void OGA :: genOutput(ofstream& outFile, const int generation) {

	outFile << generation << ",";
	for (int i = 0; i < best.size(); ++i) {
		outFile << best[i] << ",";
	}
	outFile << endl;
} // OGA :: genOutput
