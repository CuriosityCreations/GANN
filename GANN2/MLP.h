// MLP.h: header file for MLP.cpp
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

#include"oCscore.h"

// CLASS BPNET ----------------------------------------------------

/* BPNet is an implementation of the backpropagation neural network. A neural network 
has many nodes organized into layers. Each layer of nodes is completely connected to 
the adjacent layers of nodes. An input pattern is copied into the 'input layer' of the 
network (Set_Value), and these values propagate through the network, modified by the 
weights of the connections between layers (Run). The last layer of the neural network is
the output layer, and the node values in this layer represent the neural network's prediction.
During the training phase, these values are compared with a set of known values that 
correspond to the actual value of the input pattern, to determine the error in prediction 
(Set_Error). This error is then used to modify the connection weights in the network, to 
improve its predictions (LearnAll).

There are three types of nodes: input layer nodes, hidden layer nodes, and output layer nodes.
The base class BPNode is the generic backprop network node, and is used to represent the input
layer neurons, which have a single node value associated with them. The hidden and output layer 
neurons, on the other hand, have a node value and an error value, and the method of error 
calculation differs between these two types of neurons. Therefore, these types of nodes are 
derived from the BPNode class.

Finally, a link between two node in adjacent layers is represented by an object of the class 
BPLink. BPLink objects have two values: a connection strength between -1.0 and 1.0, and a 
value to 'remember' their last change, for momentum calculations.

*/

/////////////////////////////////////////////////////////////////

class BPNode {

public:

	explicit BPNode () { node_val = 0.0; this_error = 0.0; };  // Null constructor
	~BPNode() {}; // Destructor
		
	virtual void calc_error(const double param) = 0;  // Pure virtual function, calculate or set node error
	virtual double get_error() const = 0;

	double getNodeVal() const { return node_val; };  // Return the node value
	void setNodeVal(const double setMe) { node_val = setMe; };  // Set the node value

	void setNodeErr(const double setErr) { this_error = setErr; };  // Explicitly set the node error to setErr

protected:

	double node_val;  // The current value of this node
	double this_error;

};

///////////////////


class BPInNode : public BPNode {

public:

	BPInNode() {};
	~BPInNode() {};
	
	void calc_error (const double param) {};
	double get_error () const { return 0.0; };

};

///////////////////////////////////


class BPHidNode : public BPNode {

public:

	BPHidNode() {};
	~BPHidNode() {};

	void calc_error(const double param) { this_error = param; 
						// if (this_error > 1) 
						cout << "Error: " << this_error << " "; };
	double get_error() const { return this_error; };

};

//////////////////////////////////


class BPOutNode : public BPNode {

public:

	BPOutNode() {};
	~BPOutNode() {};
	
	void calc_error(const double param) { this_error = node_val * (1 - node_val) * (param - node_val); 
						// if (this_error > 1) 
						cout << "Error: " << this_error << " "; };
	double get_error() const { return this_error; };

};

////////////////////////////////////////////////////////////////


class BPLink {

public:

	BPLink() { weight = 0.0; old_delta = 0.0; };
	~BPLink() {};
	
	double get_weight() const { return weight; };
	double get_old_delta() const { return old_delta; };

	void randomize_weights (long& linkSet) { weight = ran1(&linkSet) * 2.0 - 1.0; };

	void set_weight(const double toSet) { weight = toSet; };

	void set_weight(const double lrnRate, const double moment, const double errorTimesVal) {
		double new_delta = lrnRate * errorTimesVal;
		weight += new_delta + (moment * old_delta);
		old_delta = new_delta;
	};

	void decay_weight(const double factor) { weight *= factor; };
	
private:

	double weight;
	double old_delta;

}; // Class BPLink

////////////////////////////////////////////////////////////////////
// CLASS MLP -------------------------------------------------------
//
//

class MLP {

public:

	MLP(const std::vector<double>& GAP, const int numParams, const int numLayers, const std::vector<int>& layerSizes, const vector<int> ssize, const pStr& thePars);
	
	MLP (ifstream& inNet, const int numEntries, const vector<string>& names, const vector<vector<double> >& entries) { }; 

	virtual ~MLP();

	// Create the input/output vectors for training and testing
	void buildFragments(const vector<int>& inputArray, const vector<vector<vector<double> > >& sets);

	// Determine the scaling to apply to input vectors
	void Input_Scales();

	// Reset the neural network parameters
	virtual void Reset(long& seed) = 0;

	// Prepare the neural network for training
	virtual void NNexec(const int round, const bool isTest, oCscore*& chrScore, const double eTable[EARRAYSIZE], const vector<int>& ssize) = 0;

	// Run a round of neural network training
	void NNround(const bool isGA, const bool isTest, vector<double>& FPSums, vector<double>& distSums, const double eTable[EARRAYSIZE]);

	// Run the neural network without training it
	void runNet(const int whichGroup, const int pos, const std::vector<double>& linkVals, vector<double>& FPSums, vector<double>& distSums, 
	const double eTable[EARRAYSIZE]);	

	// Train the neural network
	virtual void trainNet(const int whichGroup, const int pos, std::vector<double>& linkVals, vector<double>& FPSums, vector<double>& distSums, 
	const double eTable[EARRAYSIZE]) = 0;

	// Interpret and score the neural network's predicitons
	double thePrediction(const int whichGroup, vector<double>& FPSums, vector<double>& distSums) const;

	// scoreOne takes a single set of predictions and calculates the score that will be passed to the chromosome.
	// Try to favour false negatives (1 < -1) over false positives (1 > -1)
	double scoreOneDist (const double target) const;

	// scoreOneFP returns 1 if the prediction is on the right side of the cutoff, or 0 if not
	double scoreOneFP(const double target) const;

	// Score > 1 nodes using a distance measure
	double scoreMultDist (const int target) const;

	// Score > 1 nodes by finding the max element
	double scoreMultFP(const int target) const;

	// See if conditions for stopping the neural network are met
	virtual bool stopConditions(const int round, oCscore*& chrScore) = 0;

	// Return the number of nodes in a given layer
	int itsNumNodes(const int whichLayer) const { return num_nodes[whichLayer]; };

	// Return the number of links in the neural network
	long sumLinks() const { return link_sum; };

	// Save the best links into an array
	virtual void saveBestLinks(double*& linkSet, const oCscore* const & source, const int whichNet) = 0;

	// Set the input node values from a Pattern object
	void Set_Value(const int whichGroup, const int whichPattern, const bool isTest);
	
	// Get the value of a given output node
	double Get_Value(const int position) const { return out_node_vals[position]; };

	// Does the ANN use bias nodes?
	bool isBiased() const { return isBias; };

	// Get the value of a given link
	virtual double getLinkVal(const int which) const = 0;

	// Propagate values forward through the network
	void Run(const std::vector<double>& linkArr, const double eTable[EARRAYSIZE], const bool setHid);

	virtual void Set_Error(const int whichPattern, const bool isTest, const int whichGroup) = 0;

	// Show the initial state of the network
	virtual void showInit() = 0;

	// Save the current chromosome to disk
	virtual void saveNet(ofstream& outFile, const vector<int>& inputsToUse) = 0;
	
	// Copy link values from MLP to a vector of double
	virtual void copyTo(Chromosome*& dump) = 0;

protected:
	
	int totalRounds;
	int roundCheckScore; // The round to start checking scores
	double minScoreDif; // Minimum score difference

	vector<int> num_nodes;  // Number of nodes in each layer
	int node_sum;	 // Number of nodes in all layers

	vector<long> num_links; // Number of links between adjacent layers
	long link_sum; // Number of links in the entire neural network
	
	vector<long> link_base; // Number of links between all previous layers

	int num_layers;  // Number of layers incl. input, hidden, and output
	int layersMinusOne; // Above, minus one

	vector<double> GAprms;  // Parameters of the wrapping GA

	double * in_node_vals; // Node values for different layers 
	double * hid_node_vals; 
	double * out_node_vals;

	double * hid_node_errs; // Error values for the last two layers
	double * out_node_errs;
		
	bool iScale; // Whether or not to scale inputs
	
	int inSize; // Layer sizes to avoid array lookup in num_nodes[]
	int hidSize;
	int outSize;
	int inSet; // If (isBias), then only set (inSize - 1) node values
	int hidSet; // Ditto for hidNodes

	// Is there input vector optimization?
	bool IVO;

	// Bias nodes or no
	bool isBias;

	int numCategories; // = # of distinct sets
	int numSets; // = numCategories * 2

	// The training and test patterns
	double*** setInputs;
	vector<double> sizes;
};

////////////////////////////////////////////////////////////////////
// CLASS BPNET -----------------------------------------------------
//
//

class BPNet : public MLP {

public:

	// Constructors
	BPNet(const std::vector<double>& GAP, const int numParams, const int layers, const std::vector<int>& layerSizes, 
		const vector<int> ssize, const pStr& thePars); // Specified dimensions

	// Constructor, read from file
	BPNet(ifstream& inNet, const int numEntries, const vector<string>& names, const vector<vector<double> >& entries);

	// Destructor
	~BPNet();
	
	double getLrnRate() const { return learnRate; }; // Return the learning rate
	int howManyLinks() const { return link_sum; };	// Return the number of links in the network

	// Start the training run
	void NNexec(const int round, const bool isTest, oCscore*& chrScore, const double eTable[EARRAYSIZE], const vector<int>& ssize);

	// Train the BP net
	void trainNet(const int whichGroup, const int pos, std::vector<double>& linkVals, vector<double>& FPSums, vector<double>& distSums,
	const double eTable[EARRAYSIZE]);
		
	// Execute a loaded BPNet with a new data set
	void predictSet(vector<std::vector<double> >& thePredictions);

	// Should the neural network be stopped prematurely?
	bool stopConditions(const int round, oCscore*& chrScore);
	
	// Save the best links into an array
	void saveBestLinks(double*& linkSet, const oCscore* const & source, const int whichNet);

	void learnDecrease() { learnRate *= lrnDec; }; // Decrease the learning rate
	void weightDecrease() { for (int i = 0; i < link_sum; ++i) the_links[i].decay_weight(wtDec); };	// Decay weights

	void batchAdd(); // Store the errors and node values for later propagation
	void batchLearn(); // Batch Learn()ing

	// Does the network use batch learning?
	bool ifBatch() const { return isBatch; };

	// Show the initial configuration
	void showInit();
	
	// Get the value of a particular link
	double getLinkVal(const int which) const { return the_links[which].get_weight(); };

	// Calculate the errors for the neural network
	void Set_Error(const int whichPattern, const bool isTest, const int whichGroup);

	// Propagate errors back through the network, and change connection weights
	void Learn();
	
	// Load BPNet from disk
	void loadNet(ifstream& inNet, vector<int>& inputsToUse);

	// Save the current BPNet to disk
	void saveNet(ofstream& outFile, const vector<int>& inputsToUse);
	
	// Reset and restore the weights and parameters
	void Reset(long &seed);
	
	void copyTo(Chromosome*& dump);
	
private:

	BPLink* the_links; // All the links in the network

	bool isConj; // Learn by conjugate-gradient algorithm?
	bool isBatch; // Batch learning?
	double learnRate; // The learning rate
	double momentum; // The momentum

	double minLR; // Minimum learning rate
	int checkLR; // Round to start checking learning rate
	double minWT; // Minimum weight values

	double lrnDec; // The learning rate decrement
	double lrnDecStart; // The learning rate decrement start
	double wtDec; // The weight decrement start
	double wtDecStart; // The weight decrement start

	double* batch_errVal; // error storage for a batch network (one value per link)

};
////////////////////////////////////////////////////////////////////
// CLASS GANET -----------------------------------------------------
//
//

class GANet : public MLP, public SelectableBase {

public:

	// Constructors
	GANet(const std::vector<double>& GAP, const int numParams, const int numSets, const int numSuperSets, const int numChr, 
		const int sizeD, const int sizeB, const int sizeI, const int sizeC, 
		const double& loBound, const double& hiBound, const int numLayers, const std::vector<int>& layerSizes, 
		const vector<int> ssize, long& setSeed, const pStr& thePars); // Specified dimensions

	GANet(ifstream& inNet, const int numEntries);

	// Destructor
	~GANet() { };
	
	// Interpret scores and mess with populations of chromosomes
	virtual void optimizeParameters(const std::vector<double> GARange);
	
	// Write information from (stuff)
//	virtual void write(Writer& wrtr);

	// Start the training run
	void NNexec(const int round, const bool isTest, oCscore*& chrScore, const double eTable[EARRAYSIZE], const vector<int>& ssize);

	// Run the GA net and save the scores for later optimization
	void trainNet(const int whichGroup, const int pos, std::vector<double>& linkVals, vector<double>& FPSums, vector<double>& distSums, 
	const double eTable[EARRAYSIZE]);
		
	// Should the neural network be stopped prematurely?
	bool stopConditions(const int round, oCscore*& chrScore);
		
	// Save the best links into an array
	void saveBestLinks(double*& linkSet, const oCscore* const & source, const int whichNet);

	// Get the value of a particular link
	double getLinkVal(const int which) const { return geneVal(which); };

	// Show the initial configuration
	void showInit();

	// Calculate the errors for the neural network
	void Set_Error(const int whichPattern, const bool isTest, const int whichGroup);
	
	// Load GANet from disk
	void loadNet(ifstream& inNet, vector<int>& inputsToUse);
	
	// Save the current chromosome to disk
	void saveNet(ofstream& outFile, const vector<int>& inputsToUse);

	// Reset and restore the weights and parameters
	void Reset(long &seed);
	
	void copyTo(Chromosome*& dump);
	
private:

	vector<double> bestLrnScore; // Size = MAX_GA_ROUNDS

	int nChr; // Number of Chromosomes
	int nEvo; // Number of Evolvables
}; 

// Copy arrays
inline void copyDoubleArray(const double from[], double to[], const int size) {
	for (int i = 0; i < size; ++i) to[i] = from[i];
}

// Calculate a node value
inline int calcArr(const int div, const int signDiv, const int arrayVal) {

	return (arrayVal - (((div * div - 1) >> INTSHIFT) + 1) * ((div * ETABLESIZE) - 
		   (signDiv * ETABLESIZE - arrayVal % ETABLESIZE)) + ETABLESIZE);		   
}

