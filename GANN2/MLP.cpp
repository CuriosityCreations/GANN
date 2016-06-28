// MLP.cpp defines the multi-layer perceptron objects used in the GANN core software
// GANN CORE
// Copyright (C) 1999-2004 Robert G. Beiko

// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

#include"MLP.h"

// ABSTRACT CLASS MLP: The Multi-layer perceptron ------------------------------------------------
MLP :: MLP (const std::vector<double>& GAP, const int numParams, const int numLayers, const std::vector<int>& layerSizes, const vector<int> ssize, const pStr& thePars) {

	GAprms.resize(numParams);
	for (int copyPar = 0; copyPar < numParams; ++copyPar) {
		GAprms[copyPar] = GAP[copyPar];
	}

	if (GAP[XISBIAS] < 1.0) { isBias = false; } else { isBias = true; }

	totalRounds = elem(thePars.intGlob,"NN_TRAIN_RUNS");
	roundCheckScore = elem(thePars.intGlob,"CHECK_SCORE");
	minScoreDif = elem(thePars.doubleGlob,"NO_SCOREDIF");

	node_sum = 0;
	num_layers = numLayers;
	layersMinusOne = num_layers - 1;
	
	numSets = ssize.size();
	numCategories = numSets / 2;
	
	num_nodes.resize(numLayers); 

	// Store the number of nodes in each layer, keeping a total running count
	for (int layerCreate = 0; layerCreate < numLayers; ++layerCreate) {
		
		num_nodes[layerCreate] = layerSizes[layerCreate];
		node_sum += num_nodes[layerCreate];
	}

	inSize = num_nodes[0] + (int) isBias; // (+1)
	hidSize = num_nodes[1] + (int) isBias; // (+1)
	outSize = num_nodes[2];

	// Create node value arrays for each layer
	in_node_vals = new double[inSize];
	hid_node_vals = new double[hidSize];
	out_node_vals = new double[outSize];

	if (isBias) { // Set the bias node values, since they will never change
		in_node_vals[inSize - 1] = 1.0;
		hid_node_vals[hidSize - 1] = 1.0;
		inSet = inSize - 1;
		hidSet = hidSize - 1;
	} else {
		inSet = inSize;
		hidSet = hidSize;
	}
	
	// Create node error value arrays
	hid_node_errs = new double[hidSize];
	out_node_errs = new double[outSize];	
	
	num_links.resize(numLayers - 1);
	link_base.resize(numLayers);

	num_links[0] = inSize * hidSet; num_links[1] = hidSize * outSize;
	link_base[0] = 0; link_base[1] = num_links[0]; link_base[2] = num_links[0] + num_links[1];
	link_sum = link_base[2];

	setInputs = new double**[numSets];
	// Create space for the training and test sets
	for (int createType = 0; createType < numSets; ++createType) {
		int sz = ssize[createType];
		sizes.push_back(sz);
		setInputs[createType] = new double*[sz];
		for (int createRow = 0; createRow < sz; ++createRow) {
			setInputs[createType][createRow] = new double[inSet];
		}
	}
}

// Destructor for Multi-Layer Perceptron
MLP :: ~MLP() {

	// Destroy the NN inputs 
	for (int killType = 0; killType < numSets; ++killType) {
		int ksz = (int) sizes[killType];
		for (int killRow = 0; killRow < ksz; ++killRow) {
			delete [] setInputs[killType][killRow];
			setInputs[killType][killRow] = 0;
		}
		delete [] setInputs[killType];
		setInputs[killType] = 0;
	}
	delete [] setInputs;
	setInputs = 0;

	delete [] in_node_vals;
	in_node_vals = 0;
	delete [] hid_node_vals;
	hid_node_vals = 0;
	delete [] out_node_vals;
	out_node_vals = 0;
	
	delete [] hid_node_errs;
	hid_node_errs = 0;
	delete [] out_node_errs;
	out_node_errs = 0;
}

// buildFragments copies the subset of input values that will be used into the MLP array
void MLP :: buildFragments(const vector<int>& inputArray, const vector<vector<vector<double> > >& sets) {

	int sz;

	for (int setType = 0; setType < numSets; ++setType) {
		sz = (int) sizes[setType];
		for (int eachMember = 0; eachMember < sz; ++eachMember) {
			for (int eachInput = 0; eachInput < inSet; ++eachInput) {
				setInputs[setType][eachMember][eachInput] = sets[setType][eachMember][inputArray[eachInput]];
			}
		}
	}
}

void MLP :: NNround (const bool isGA, const bool isTest, vector<double>& FPSums, vector<double>& distSums, const double eTable[EARRAYSIZE]) {

	const int linkNum = sumLinks();

	std::vector<double> thisNetLinks;
	thisNetLinks.resize(linkNum);
	for (int copyLink = 0; copyLink < linkNum; ++copyLink) {
		thisNetLinks[copyLink] = getLinkVal(copyLink);
	}

	double FPScore = 0.0, distScore = 0.0; // Scores
	
	int outNodes = (int) GAprms[XOUTNODE]; // Number of output nodes

	int arrayRound;

	int i;
	if (!isTest) { // training set, therefore use trainNet
		for (i = 0; i < numCategories; ++i) {
			for (arrayRound = 0; arrayRound < sizes[i]; ++arrayRound) {
				trainNet(i, arrayRound, thisNetLinks, FPSums, distSums, eTable);
			} 
		}
	} else { // Test set, therefore use runNet
		for (i = numCategories; i < numSets; ++i) {
			for (arrayRound = 0; arrayRound < sizes[i]; ++arrayRound) {	
				runNet(i, arrayRound, thisNetLinks, FPSums, distSums, eTable);
			}
		}
	} 
} // NNround

// Universal 'run' function that does not train the neural network
void MLP :: runNet(const int whichGroup, const int pos, const std::vector<double>& linkVals, vector<double>& FPSums, vector<double>& distSums, 
	const double eTable[EARRAYSIZE]) {

	double scoreNet = 0.0;

	Set_Value(whichGroup, pos, true);

	Run(linkVals, eTable, false);

	scoreNet = thePrediction(whichGroup, FPSums, distSums);
// cout << out_node_vals[0] << "," << flush;
} // MLP :: runNet

// thePrediction calculates the neural network's prediction and returns the distance score
double MLP :: thePrediction(const int whichGroup, vector<double>& FPSums, vector<double>& distSums) const {

	double FP;

	if (outSize == 1) { // A single output node for two groups
		double target = (whichGroup % 2 == 1) ? 1.0 : -1.0;
		FP = scoreOneFP(target);
		FPSums[whichGroup] += FP;
		distSums[whichGroup] += scoreOneDist(target);
	} else { // Corresponds to number of output nodes
		FP = scoreMultFP(whichGroup);
		FPSums[whichGroup] += FP;
		distSums[whichGroup] += scoreMultDist(whichGroup);
	}
	return FP; // The score used for selection
} // MLP :: thePrediction()

// Score a single output node using the distance between predicted and expected
double MLP :: scoreOneDist (const double target) const {
	double pred = out_node_vals[0];
	double theScore = ((pred > target) ? (2 - pred + target) : (2 - target + pred));
	return theScore;	
}

// Score a single output node using a cutoff
double MLP :: scoreOneFP (const double target) const {
	double pred = out_node_vals[0];
	return ((target > 0 && pred > 0) || (target <= 0 && pred <= 0)) ? 1.0 : 0.0;
}

// Score > 1 nodes using a distance measure
double MLP :: scoreMultDist (const int target) const {

	double sumDist = 0.0;

	for (int i = 0; i < numCategories; ++i) {
		sumDist += ((i == target % numCategories) 
			? (1.0 - out_node_vals[i]) * (numCategories - 1) // This will weight the positive case equally with all of the negative cases
			: 1.0 + out_node_vals[i]);
	}
	double theScore = sumDist / numCategories * 2 - 2;
	return theScore;
}

// Score > 1 nodes by finding the max element
double MLP :: scoreMultFP(const int target) const { 
	
	int max_i = 0;
	double max_val = -1.0;
	
	for (int i = 0; i < numCategories; ++i) {
		if (out_node_vals[i] > max_val) {
			max_i = i;
			max_val = out_node_vals[i];
		}
	}
	return ((max_i == target % numCategories) ? 1.0 : 0.0); // This one doesn't need to be weighted
}

// Set_Value copies the values from a pattern object into the input layer of the 
// neural network.
void MLP :: Set_Value(const int whichGroup, const int whichPattern, const bool isTest) {

	for (int trainRun = 0; trainRun < inSet; ++trainRun) {
		in_node_vals[trainRun] = setInputs[whichGroup][whichPattern][trainRun];
	}
} // Set_Value

// Run propagates node values forward through the network. The node values in a given 
// layer correspond to the the product of each node in the previous layer with the 
// connection weight leading from that node, summed together, then passed through a 
// 'transfer function' to keep the node values between 0.0 and 1.0.

// NOTE: Run for backprop differs from Run for GA. For backprop, the hidden node values must be set,
// whereas for GA (and for testing), the variables are merely stored in an array.

void MLP :: Run(const std::vector<double>& linkArr, const double eTable[EARRAYSIZE], const bool setHid) {
	
	double runningTotal = 0.0;  // The sum of (node value * connection weight)
	const int inBase = link_base[0];
	const int hidBase = link_base[1];
	int runningBase = 0;
	int nodeValSet = 0;	
	int inNode = 0;
	int arrayVal = 0; // $$$
	int div = 0;
	int signDiv = 0;
		
	// Determine the value for each node in the hidden layer...
	for (nodeValSet = 0; nodeValSet < hidSet; ++nodeValSet) {

		runningTotal = 0.0;
		
		runningBase = inBase + nodeValSet;
		
		// ...By iterating each node in the input layer
		for (inNode = 0; inNode < inSize; ++inNode) {

			// Add the (node value * connection weight) to runningTotal
			runningTotal += in_node_vals[inNode] * linkArr[runningBase];
			runningBase += hidSize; 
		}
		// When finished with this node in the current layer, pass its value through
		// the transfer function

		arrayVal = (int) (runningTotal * EDIVTEN);
		if (arrayVal > EMINUSONE) arrayVal = EMINUSONE;
		else if (arrayVal < NEGEMONE) arrayVal = NEGEMONE;
		hid_node_vals[nodeValSet] = eTable[arrayVal + ETABLESIZE]; // Set the hidden node value
	}
	
	// Determine the value for each node in the output layer...
	for (nodeValSet = 0; nodeValSet < outSize; ++nodeValSet) {

		runningTotal = 0.0;
		
		runningBase = hidBase + nodeValSet;
		
		// ...By iterating each node in the hidden layer
		for (inNode = 0; inNode < hidSize; ++inNode) {

			// Add the (node value * connection weight) to runningTotal
			runningTotal += hid_node_vals[inNode] * linkArr[runningBase];
			runningBase += outSize; 
		}


		// When finished with this node in the current layer, pass its value through
		// the transfer function

		arrayVal = (int) (runningTotal * EDIVTEN);
		if (arrayVal > EMINUSONE) arrayVal = EMINUSONE;
		else if (arrayVal < NEGEMONE) arrayVal = NEGEMONE;
		out_node_vals[nodeValSet] = eTable[arrayVal + ETABLESIZE]; // Set node value
	}
} // Run

// DERIVED CLASS BPNET: The Backpropagation Neural Network ---------------------------------------
BPNet :: BPNet (const std::vector<double>& GAP, const int numParams, const int numLayers, const std::vector<int>& layerSizes, const vector<int> ssize, const pStr& thePars) 
							: MLP (GAP, numParams, numLayers, layerSizes, ssize, thePars) {

	long ran_link = 145;
 
	the_links = new BPLink[link_sum];
	for (int q = 0; q < link_sum; ++q) the_links[q].randomize_weights(ran_link);
	
	if (GAP[XBATCH] >= 1.0) {
		isBatch = true;
		batch_errVal = new double[link_sum];
		for (int sZ = 0; sZ < link_sum; ++sZ)
			batch_errVal[sZ] = 0.0;
	}
	else isBatch = false;
	
	minLR = elem(thePars.doubleGlob,"LR_TOOLOW");
	checkLR = elem(thePars.intGlob,"LR_CHECKROUND");
	minWT = elem(thePars.doubleGlob,"TINY_WEIGHT");

	learnRate = GAP[XLRNRATE];
	momentum = GAP[XMOMENT];

	lrnDec = 1.0 - GAP[XLRNDECAY];
	wtDec = 1.0 - GAP[XWEIGHTDECAY];

	lrnDecStart = (int) ((double) totalRounds * GAP[XLRNDECAYSTART]); // Start LR decay
	wtDecStart = (int) ((double) totalRounds * GAP[XWTSTART]); // Start weight decay
}

// Read a predefined neural network from a file
BPNet :: BPNet(ifstream& inNet, const int numEntries, const vector<string>& names, const vector<vector<double> >& entries) 
							: MLP (inNet, numEntries, names, entries) {

	isBatch = false;

	string inLine;
	// First line: ignored
	getline (inNet, inLine);
	// Second line: index names
	getline (inNet, inLine);
	vector<string> theNames = split(inLine, ",");

	vector<int> whereNames(theNames.size());
	for (int i = 0; i < theNames.size(); ++i) {
		bool thisFound = false;
		string searchFor = theNames[i];
		for (int j = 0; !thisFound; ++j) {
			if (searchFor == names[j]) {
				whereNames[i] = j;
				thisFound = true;
			} else {
				if (j == names.size()) {
					cout << "Name " << searchFor << " not found in index file." << endl;
					exit (-1);
				}
			}
		}
	}
		
	// Next line: bias nodes or not
	getline (inNet, inLine);
	isBias = (inLine == "BIAS") ? true : false;
	
	// Next line: set sizes
	getline (inNet, inLine);
	vector<string> nSizes = split(inLine, ",");
	inSize = atoi((nSizes[0]).c_str()); hidSize = atoi((nSizes[1]).c_str()); outSize = atoi((nSizes[2]).c_str());

	if (theNames.size() != inSize) {
		cout << "Mismatch in ANN file between input names and # inputs." << endl;
		exit(-1);
	}

	in_node_vals = new double[inSize];
	hid_node_vals = new double[hidSize];
	out_node_vals = new double[outSize];

	// Create node error value arrays (for destructor compatibility more than anything else)
	hid_node_errs = new double[hidSize];
	out_node_errs = new double[outSize];

	if (isBias) { // Set the bias node values, since they will never change
		in_node_vals[inSize - 1] = 1.0;
		hid_node_vals[hidSize - 1] = 1.0;
		inSet = inSize - 1;
		hidSet = hidSize - 1;
	} else {
		inSet = inSize;
		hidSet = hidSize;
	}
	
	int numLayers = 3;

	num_links.resize(numLayers - 1);
	link_base.resize(numLayers);

	num_links[0] = inSize * hidSet; num_links[1] = hidSize * outSize;
	link_base[0] = 0; link_base[1] = num_links[0]; link_base[2] = num_links[0] + num_links[1];
	link_sum = link_base[2];

	numSets = 1;

	sizes.push_back(numEntries);

	setInputs = new double**[numSets]; // Only a single set (we don't know the classification)
	// Create space for the entries and set the input values
	setInputs[0] = new double*[numEntries];
	for (int createRow = 0; createRow < numEntries; ++createRow) {
		setInputs[0][createRow] = new double[inSet];
		for (int putVal = 0; putVal < inSet; ++putVal) { // Set the values for each row
			setInputs[0][createRow][putVal] = entries[createRow][whereNames[putVal]];
		}
	}
	
	// Last line: link values
	the_links = new BPLink[link_sum];

	getline (inNet, inLine);
	vector<string> linkText = split(inLine, ",");
	for (int setLink = 0; setLink < link_sum; ++setLink) {
		(the_links[setLink]).set_weight(atof((linkText[setLink]).c_str()));
	}
}

BPNet :: ~BPNet() {

	delete [] the_links;
	the_links = 0;
	if (isBatch) {
		delete [] batch_errVal;
		batch_errVal = 0;
	}
}


void BPNet :: NNexec(const int round, const bool isTest, oCscore*& chrScore, const double eTable[EARRAYSIZE], const vector<int>& ssize) {

	const double n = (!isTest) ? ssize[0] + ssize[1] : ssize[2] + ssize[3];
	vector<double> FPSums(numSets, double(0.0)); // Vector for holding sums
	vector<double> distSums(numSets, double(0.0));

	NNround(false, isTest, FPSums, distSums, eTable);

	// Record average scores and standard deviations in the oNetScore object
	double FPScore = 0.0;
	double distScore = 0.0;

	int start = (isTest) ? numCategories : 0;  
	int finish = start + numCategories;

	for (int i = start; i < finish; ++i) {
		FPScore += FPSums[i] / sizes[i]; // Correct for set size
		distScore += distSums[i] / sizes[i];
	}
	FPScore /= numCategories;
	distScore /= numCategories;

	chrScore->setScoreDev(start, round, FPScore);
	chrScore->setScoreDev(start + 1, round, distScore);

	if (ifBatch()) batchLearn();

	chrScore->findBestScore(isTest);

	if (round > lrnDecStart) learnDecrease();  // Decrease the learning rate if appropriate
	if (round > wtDecStart) weightDecrease();  // Weight decay

} // BPNet :: NNexec

// Train a backpropagation neural network
void BPNet :: trainNet(const int whichGroup, const int pos, std::vector<double>& linkVals, vector<double>& FPSums, vector<double>& distSums, 
	const double eTable[EARRAYSIZE]) {

	double scoreNet = 0.0;

	for (int copyLink = 0; copyLink < link_sum; ++copyLink) {
		linkVals[copyLink] = getLinkVal(copyLink);
	}

	Set_Value(whichGroup, pos, false); 

	Run(linkVals, eTable, true);

	Set_Error(pos, false, whichGroup);

	if (!ifBatch()) {
		Learn();  // Do not pass linkVals, since they will be modified
	} else {
		batchAdd();
	}

	scoreNet = thePrediction(whichGroup, FPSums, distSums);

} // BPNet :: trainNet

void BPNet :: predictSet(vector<vector<double> >& thePredictions) {

	// Build the exp lookup table
	double eTable[ETABLESIZE * 2];
	for (int eBuild = 0; eBuild < ETABLESIZE * 2; ++eBuild) {
		eTable[eBuild] = 2.0/(1.0 + exp(-2 * (ETABLESIZE - eBuild) / 1000.0)) - 1;  // Equivalent to negative (eBuild - 10000.0)
	}

	// Set the link values
	std::vector<double> linkVals(link_sum);
	for (int copyLink = 0; copyLink < link_sum; ++copyLink) {
		linkVals[copyLink] = getLinkVal(copyLink);
	}
	// Run the net for each 
	for (int i = 0; i < thePredictions.size(); ++i) {
		Set_Value(0, i, true);
		Run(linkVals, eTable, false);
		(thePredictions[i]).resize(outSize);
		for (int setOut = 0; setOut < outSize; ++setOut) { // Store predictions in the big vector of vectors
			thePredictions[i][setOut] = out_node_vals[setOut];
		}
	}
}

// Stop the neural network if any of the stopping conditions are met
bool BPNet :: stopConditions(const int round, oCscore*& chrScore) {

	// Progress bar
	int tick = (int) ((double) totalRounds / 10.0 + 0.9);
	if (totalRounds >= tick) {
		if ((round + 1) % tick == 0) { cout << "*"; cout.flush(); }
	}

	// Stop the run if the last generation has been reached
	if (round >= totalRounds) return true;

	// Stop the run if the learning rate has decreased below a threshold
	if ((getLrnRate() < minLR) && (round > checkLR)) { chrScore->setWhyStop(LR_CRASH); return true; }

	// Stop the run if NN weights have crashed
	int weight_check = 0;
	bool all_zero = true;
	while (weight_check < howManyLinks()) {
		if (fabs(getLinkVal(weight_check)) > minWT) {
			all_zero = false;
			break;
		}
		++weight_check;
	}
	if (all_zero) {
	chrScore->setWhyStop(WT_CRASH); cout << "Crash!" << endl; return true; 
	}

	// Stop the run if no generalization improvement has occurred in the last 20 training rounds
	if (round > roundCheckScore) {
		if (fabs(chrScore->getScoreDev(TESTFPSCORE, round) - chrScore->getScoreDev(TESTFPSCORE, round - roundCheckScore)) < minScoreDif) { 
			chrScore->setWhyStop(NO_LEARN); return true;
		}
	}

	// If the net has made it this far, it can go another round
	return false;
}

// Save the best links into an array
void BPNet :: saveBestLinks(double*& linkSet, const oCscore* const & source, const int whichNet) {

	for (int i = 0; i < link_sum; ++i) {
		linkSet[i] = the_links[i].get_weight();
	}
}

// For batch learning, calculate the error * value term but do not modify connection weights
void BPNet :: batchAdd() {

	const int lBase = link_base[layersMinusOne - 1];

	// Add the error * nodeVal to batch_errVal for each link between input and hidden layer nodes
	for (int whichLink = 0; whichLink < lBase; ++whichLink) {
		batch_errVal[whichLink] += in_node_vals[whichLink / hidSize] * hid_node_errs[whichLink % hidSize];
	}
	
	int lFix = 0;
	// Add for each link between hidden and output layer nodes
	for (int wLink = lBase; wLink < link_sum; ++wLink) {
		lFix = wLink - lBase;
		batch_errVal[wLink] += hid_node_vals[lFix / outSize] * out_node_errs[lFix % outSize];
	}
}

// Change connection weights for batch learning, and reset batch_errVal to 0
void BPNet :: batchLearn() {

	for (int linkIt = 0; linkIt < link_sum; ++linkIt) {
		(the_links[linkIt]).set_weight(learnRate, momentum, batch_errVal[linkIt]);
	}

	for (int backToZero = 0; backToZero < link_sum; ++backToZero) {
		batch_errVal[backToZero] = 0.0;
	}

}

// Show the initial configuration
void BPNet :: showInit()  {

	for (int i = 0; i < link_sum; ++i) {
		cout << i << " " << the_links[i].get_weight() << endl;
	}
	cout << "------------------=====-=-==-=-=-=-=-=--=-=-======------------------" << endl; 
} // showInit()

// Set_Error uses the desired output pattern from a Pattern object to calculate the error
// in the output layer nodes.
void BPNet :: Set_Error (const int whichPattern, const bool isTest, const int whichGroup) {

	const int layersMinusTwo = layersMinusOne - 1;
	const int lBase = link_base[layersMinusTwo];	
	double theValue = 0.0;
	
	vector<double> target(outSize, double(-1.0));
	if (outSize == 1) { 
		if (whichGroup == 1) { target[0] = 1.0; } 
	} else {	
		target[whichGroup] = 1.0;
	}

	if (!isTest) {
		// Set output layer error
		for (int errRun = 0; errRun < outSize; ++errRun) {
			theValue = out_node_vals[errRun];
			out_node_errs[errRun] = (1 - theValue * theValue) * (target[errRun] - theValue);
		}
	} else {
		// Set output layer error
		for (int errRun = 0; errRun < outSize; ++errRun) {
			theValue = out_node_vals[errRun];
			out_node_errs[errRun] = (1 - theValue * theValue) * (target[errRun] - theValue);
		}
	}
	// Set hidden layer error for all nodes
	for (int hErrRun = 0; hErrRun < hidSet; ++hErrRun) {
		double hNodeVal = hid_node_vals[hErrRun];
		double sum_err = 0;
		int loh = lBase + hErrRun * outSize;
		for (int outN = 0; outN < outSize; ++outN) {
			sum_err += the_links[loh + outN].get_weight() *
				out_node_errs[outN];
		}
		sum_err *= (1 - hNodeVal * hNodeVal);

		hid_node_errs[hErrRun] = sum_err;
// if (hid_node_errs[hErrRun] > 100.0) cout << "Bighid! " << hNodeVal << " " << sum_err << " ";
	}
	 
}  // Set_Error

// Do the dirty work. Propagate errors back through the network, changing the appropriate link values
// (Output node error was already calculated by calc_error)
void BPNet :: Learn() {

	const int lBase = link_base[layersMinusOne - 1];
	// Change the input-hidden links based on the formula in BPLink :: set_weight
	for (int linkChange = 0; linkChange < lBase; ++linkChange) {
		(the_links[linkChange]).set_weight(learnRate, momentum, hid_node_errs[linkChange % 
			hidSize] * in_node_vals[linkChange / hidSize]);
	}

	int lFix = 0;
	// Change the hidden-output links based on the same formula
	for (int lChange = lBase; lChange < link_sum; ++lChange) {
		lFix = lChange - lBase;
		(the_links[lChange]).set_weight(learnRate, momentum, out_node_errs[lFix %
			outSize] * hid_node_vals[lFix / outSize]);
	}
} // Learn()

// This function does not currently do much of anything.
void BPNet :: loadNet(ifstream& inNet, vector<int>& inputsToUse) {

	char inputBuf[80];
	int shouldBuf[3] = {0};
	int shouldBeLinks = 0;

	char* toss;
	
	int inLinkCount = 0;
	
	std::vector<double> replacing;
	replacing.resize(link_sum);

	inNet.getline(inputBuf, 80, '\n');
	shouldBuf[0] = strtol(strtok(inputBuf, " "), &toss, 10);
	shouldBuf[1] = strtol(strtok(NULL, " "), &toss, 10);
	shouldBuf[2] = strtol(strtok(NULL, " "), &toss, 10);
	
	shouldBeLinks = shouldBuf[0] * shouldBuf[1] + shouldBuf[0] * shouldBuf[1];
		
	while (inLinkCount < shouldBeLinks) {
		inNet.getline(inputBuf, 80, ' ');
		replacing[inLinkCount] = strtod(inputBuf, &toss);
		inLinkCount++;
	}
} //loadNet

// Save the neural network to outFile
void BPNet :: saveNet(ofstream& outFile, const vector<int>& inputsToUse) {

	for (int nOut = 0; nOut < num_layers; ++nOut) {
		outFile << num_nodes[nOut] << ",";
	}
	outFile << endl;

	for (int lOut = 0; lOut < link_sum; ++lOut) {
		outFile << getLinkVal(lOut) << ",";
	}
	outFile << endl << endl;
} // saveNet

// Reset randomizes the connection weights and restores the initial parameters of an existing neural network
void BPNet :: Reset(long &seed) {

	for (int i = 0; i < link_sum; ++i) the_links[i].randomize_weights(seed);

	learnRate = GAprms[XLRNRATE];
	momentum = GAprms[XMOMENT];
	lrnDec = 1.0 - GAprms[XLRNDECAY];
} // Reset

// Copy to the 'best' chromosome
void BPNet :: copyTo(Chromosome*& dump) {

	delete dump;
	dump = 0;

	double * allVals = new double[link_sum];
	for (int i = 0; i < link_sum; ++i) {
		allVals[i] = getLinkVal(i);
	}

	dump = new Chromosome(link_sum, allVals);

	delete [] allVals;
	allVals = 0;
}


// DERIVED CLASS GANET: The Genetic Algorithm Neural Network -------------------------------------
GANet :: GANet (const std::vector<double>& GAP, const int numParams, const int numSets, const int numSuperSets, const int numChr, 
					const int sizeD, const int sizeB, const int sizeI, const int sizeC, 
					const double& loBound, const double& hiBound, const int numLayers, const std::vector<int>& layerSizes, const vector<int> ssize, long& setSeed, const pStr& thePars) 
					: MLP (GAP, numParams, numLayers, layerSizes, ssize, thePars) 
					, SelectableBase (numSets, numSuperSets, numChr, sizeD, sizeB, sizeI, sizeC,
					loBound, hiBound, setSeed, thePars) {

	bestLrnScore.resize(totalRounds); // Or intGlob["MAX_GA_ROUNDS?"]

	nChr = numChr;
	nEvo = numSets;
}

// Do the dirty work. Interpret the scores obtained with the neural nets, and perform X-over, mutation, 
// and migration
void GANet :: optimizeParameters (const std::vector<double> GARange) {

	int r;
	for (int i = 0; i < numParamSets_; ++i) {
		if (myParams_[i]->sortBySuccess()) {
			double mutRate=GAprms[XMUTRATE];
			double mutAmt=GAprms[XMUTAMT];
			double mutProp=GAprms[XMUTPROP];
			double recRate=GAprms[XRECRATE];
			double recRepl=GAprms[XRECREPL];
			double migRate=GAprms[XMIGRATE];
			double killFrac=GAprms[XKILLFRAC];
			double killDiff=GAprms[XKILLDIFF];
			
			myParams_[i]->recombine(recRate, recRepl);	// Might consider optimizing the rate of evolution as well
			myParams_[i]->mutate(mutRate, mutAmt, (int) (mutProp * (double) myParams_[i]->itsChromosomeSize()), false); // though this might be dangerous
			if (numParamSets_ > 1 && rnd(1.0) < migRate) {
				r = rnd(numParamSets_);
				while (r == i) r = rnd(numParamSets_);
				myParams_[i]->acceptImmigrants(myParams_[r + currentSuperSet_ * numParamSets_], 0.05, maybe());
					// The immigrants replace some of the mutants
			};
			//myParams_[i]->killClones(killFrac, killDiff);
		}
	};
}

// NNexec iterates through the inner Chromosomes of the genetic algorithm, preparing them for scoring
void GANet :: NNexec(const int round, const bool isTest, oCscore*& chrScore, const double eTable[EARRAYSIZE], const vector<int>& ssize) {

	const double n = (!isTest) ? ssize[0] + ssize[1] : ssize[2] + ssize[3];

	int netSel = 0;
	int netEvol = 0;
	int netChr = 0;
	int numEvo = (int) GAprms[XNUMEVO];
	int numChr = (int) GAprms[XNUMCHR];
	
	double maxFP = 0.0;
	int bestChr = 0;

	while (netSel < 1) {

		while (netEvol < numEvo) {
		
			setParamSet(netEvol);
		
			while (netChr < numChr) {

				chooseChromosome(netChr);

				const double n = (!isTest) ? ssize[0] + ssize[1] : ssize[2] + ssize[3];
				vector<double> FPSums(numSets, double(0.0)); // Vector for holding sums
				vector<double> distSums(numSets, double(0.0));

				NNround(true, isTest, FPSums, distSums, eTable);
				
				int index = netEvol * numChr + netChr;
				
				double FPScore = 0.0;
				double distScore = 0.0;

				int start = (isTest) ? numCategories : 0;  
				int finish = start + numCategories;
				
				for (int i = start; i < finish; ++i) {
					FPScore += FPSums[i] / sizes[i]; // Correct for set size
					distScore += distSums[i] / sizes[i];
				}
				FPScore /= numCategories;
				distScore /= numCategories;
								
				// Record average scores and standard deviations in the oNetScore object
				chrScore->setScoreDev(start++, index, FPScore);
				chrScore->setScoreDev(start++, index, distScore);

				if (FPScore > maxFP) { maxFP = FPScore; bestChr = netChr; }

				netChr++;
			}
			netEvol++;
			netChr = 0;
		}
		netSel++;
		netEvol = 0;
	}
			
	// Record scores
	chrScore->findBestScore(isTest);
	chooseChromosome(bestChr);

/*	if (isTest) { // If this is a testing round
		chrScore->doTheMaxMatrices(round); // Score the outer Chromosome for this training round
		chrScore->resultsOutput(true);
	} else { // If this is NOT a testing round
		if ((round + 1) % 10 != 0) {
			chrScore->doTheTrainMatrices(round);
			chrScore->resultsOutput(false);
		}
	} */

	const std::vector<double> BPRange;

	// NOW evolve (Chromosomes are rearranged, must perform this last)
	if (!isTest) optimizeParameters(BPRange);
} // GANet :: NNexec

// Score a genetic algorithm neural network
void GANet :: trainNet(const int whichGroup, const int pos, std::vector<double>& linkVals, vector<double>& FPSums, vector<double>& distSums, 
	const double eTable[EARRAYSIZE]) {

	double scoreNet = 0.0;

	Set_Value(whichGroup, pos, false);

	Run(linkVals, eTable, false); 

	scoreNet = thePrediction(whichGroup, FPSums, distSums);
	
	chooseChromosome(whichChrSelected(), scoreNet); 
} // GANet :: trainNet

// Determine whether the training run should be stopped (return true if so)
bool GANet :: stopConditions (const int round, oCscore*& chrScore) {

	// Progress bar
	//if ((round + 1) % (totalRounds / 10) == 0) { cout << "*"; cout.flush(); }

	int tick = (int) ((double) totalRounds / 10.0 + 0.9);
	if (totalRounds >= tick) {
		if ((round + 1) % tick == 0) { cout << "*"; cout.flush(); }
	}

	// Stop the run if the last generation has been reached
	if (round >= totalRounds) return true;
	
	bestLrnScore[round] = chrScore->getBestScore(0);
	
	// Stop the run if training score has not changed substantially in the past 20 rounds
	if (round > roundCheckScore) {
		if (fabs(bestLrnScore[round] - bestLrnScore[round - roundCheckScore]) < minScoreDif) { 
			chrScore->setWhyStop(NO_LEARN); return true; 
		}
	}
	
	// Stop the run if NN weights have crashed
	int weight_check = 0;
	bool all_zero = true;
	while (weight_check < sumLinks()) {
		if (getLinkVal(weight_check) > 0.001) {
			all_zero = false;
			break;
		}
		++weight_check;
	}
	if (all_zero) { chrScore->setWhyStop(WT_CRASH); return true; }

	// If the net has made it this far, it can go another round
	return false;

} // GANet :: stopConditions

// Save the neural network in a standard format
void GANet :: saveNet (ofstream& outFile, const vector<int>& inputsToUse) {

	outFile << "BEGIN" << endl; // Start of record
	for (int i = 0; i < inSet; ++i) {
		outFile << inputsToUse[i] << ","; // Which inputs from file are being used
	}
	outFile << endl;
	myParams_[selectedParamSet_]->saveChr(outFile, num_links, inputsToUse); // Connection weights
	outFile << "END" << endl; // End of record
}

// Save the best links into an array
void GANet :: saveBestLinks(double*& linkSet, const oCscore* const & source, const int whichNet) {
/*
	int bestIndex = 0;

	if (whichNet < 2) { // distance scores
		if (whichNet == 0) {
			bestIndex = (int) source->getMaxTrainDist(3);
		} else {
			bestIndex = (int) source->getMaxTestDist(3);
		}
	} else {
		if (whichNet == 2) {
			bestIndex = (int) source->getMaxTrainFP(3);
		} else {
			bestIndex = (int) source->getMaxTestFP(3);
		}
	}
	setParamSet(bestIndex / nChr);
	chooseChromosome(bestIndex % nChr);
	
	for (int i = 0; i < link_sum; ++i) {
		linkSet[i] = gene(i);
	}	*/
}

// Show the initial configuration
void GANet :: showInit() {

	for (int zeke = 0; zeke < numParamSets_; ++zeke) {
		setParamSet(zeke);
		for (int jed = 0; jed < myParams_[zeke]->itsNumChromosomes(); ++jed) {
			myParams_[zeke]->chooseChromosome(jed);
			cout << "A  " << geneVal(0) << endl;
		}
		cout << "------------------=====-=-==-=-=-=-=-=--=-=-======------------------" << endl;
	}
}

// Set_Error uses the desired output pattern from a Pattern object to calculate the error
// in the output layer nodes.
void GANet :: Set_Error (const int whichPattern, const bool isTest, const int whichGroup) {

	double theValue;

	vector<double> target(outSize, double(-1.0));
	if (outSize == 1) { 
		if (whichGroup == 1) { target[0] = 1.0; } 
	} else {	
		target[whichGroup] = 1.0;
	}
	
	if (!isTest) {
		// Set output layer error
		for (int errRun = 0; errRun < outSize; ++errRun) {
			theValue = out_node_vals[errRun];
			out_node_errs[errRun] = theValue * (1 - theValue) * (target[errRun] - theValue);
		}
	} else {
		// Set output layer error
		for (int errRun = 0; errRun < outSize; ++errRun) {
			theValue = out_node_vals[errRun];
			out_node_errs[errRun] = theValue * (1 - theValue) * (target[errRun] - theValue);
		}
	}	
	
}  // Set_Error 

void GANet :: loadNet(ifstream& inNet, vector<int>& inputsToUse) {

/*	char inputBuf[80];
	int shouldBuf[3] = {0};
	int shouldBeLinks = 0;

	char* toss;
	
	int inLinkCount = 0;
	
	std::vector<GeneBase*> replacing;
	replacing.resize(link_sum);

	for (int i = 0; i < N_GA_PARAMS; ++i) {
		replacing[i] = new GeneDouble();
	}
	int tSize = N_GA_PARAMS + NUM_INPUTS;
	for (int j = N_GA_PARAMS; j < tSize; ++j) {
		replacing[j] = new GeneBool();
	}
	inNet.getline(inputBuf, 80, '\n');
	shouldBuf[0] = strtol(strtok(inputBuf, " "), &toss, 10);
	shouldBuf[1] = strtol(strtok(NULL, " "), &toss, 10);
	shouldBuf[2] = strtol(strtok(NULL, " "), &toss, 10);
	
	shouldBeLinks = shouldBuf[0] * shouldBuf[1] + shouldBuf[0] * shouldBuf[1];
		
	while (inLinkCount < shouldBeLinks) {
		inNet.getline(inputBuf, 80, ' ');
		replacing[inLinkCount]->forceVal(strtod(inputBuf, &toss));
		inLinkCount++;
	}
	
	replaceGenes(replacing);
	
	for (int k = 0; k < tSize; ++k) {
		delete replacing[k];
		replacing[k] = 0;
	} */
} // loadNet

// Randomize the populations
void GANet :: Reset(long &seed) {
	// for (int i = 0; i < link_sum; ++i) the_links[i].randomize_weights(seed);
} // Reset

// Write information from (stuff). Don't know if I'll use this
// void GANet :: write(Writer& wrtr) {}

void GANet :: copyTo(Chromosome*& dump) {

	delete dump;
	dump = 0;
	dump = new Chromosome(myParams_[selectedParamSet_]->getChr());
}
