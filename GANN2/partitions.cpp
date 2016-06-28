// partitions.cpp reads the input files to GANN, and creates the OGA object.
// GANN CORE
// Copyright (C) 1998-2004 Robert G. Beiko

// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

#include"OGA.h"

// ---------------------------------- partitions.cpp -----------------------------------------

/*
Partitions.cpp reads an input file containing all members of the positive and negative members of the training and test sets.
These sets are represented by a number of values, describing such properties as nucleotide composition, helical properties, 
etc. of DNA. The program will set up an outer genetic algorithm (OGA) which will specify a numbre of training instances, each 
representing 20 out of all possible neural network inputs. The best predictors among all inputs will recombine, leading to
a population of neural networks that use effective input values to classify DNA sequences as regulatory or non-regulatory.

*/

// Other functions -----------------------------------------------------------------


void detect(const bool NNisGA, const string inFileName, const string inCommand);
void runOnce(const string& NNfile, const string& inputfile, const string& outfile);
void readInFile(ifstream& NNfile, vector<string>& names, vector<string>& seqIDs, vector<vector<double> >& entries);
void writePreds(ofstream& outPreds, const vector<string>& seqIDs, const vector<vector<double> >& outNodePredictions);
void readCommand(const string inCommand, pStr& thePars);
void setDefaults(pStr& thePars);

// ------------------------------------- MAIN ----------------------------------------

int main(int argc, char* argv[]) {

	bool runSuccess = false;

	if ((argc == 4) && ((argv[1][1] == 'b') || (argv[1][1] == 'g'))) { // Check # of params and ID of second param
		detect(
			(argv[1][1] == 'g') ? true : false, // GA or NN
			argv[2], // Input file name
			argv[3]); // Command file name
			runSuccess = true;
	} else {
		if ((argc == 5) && (argv[1][1] == 'r')) {
			runOnce(
				argv[2], // Trained NN
				argv[3], // Input indices
				argv[4]); // Output file name
				runSuccess = true;
		}
	}
	if (!runSuccess) {
		cout << U_MESSAGE << argv[0] << U_FLAGS << endl;
	} 
	cout << "Done!" << endl;
	return 0;
} // MAIN


// ------------------------- Function Definitions -----------------------

// Detect initializes and runs the outer genetic algorithm for a number of generations. 
void detect (const bool NNisGA, const string inFileName, const string inCommand) {

	pStr thePars;

	srand((unsigned) time(NULL));
	long setSeed = 0 - rand();

	// Initialize the random number generator
	ran1(&setSeed);

	short numParams;
	if (NNisGA) numParams = N_GA_PARAMS;
	else numParams = N_BP_PARAMS;

	// Read the command file to set global constants
	readCommand(inCommand, thePars); 

	// Define allowable range of GA parameters
	vector<double> GARange;
	GARange.resize(numParams * 2);
	setGAVec(GARange, NNisGA, thePars);
	
	// Read the number of inputs from the input file
	string header;
	ifstream inFile;
	inFile.open(inFileName.c_str());

	if (!inFile.is_open()) {
		cout << "Can't open index file." << endl;
		exit(-1);
	}
	getline(inFile, header);
	inFile.close();

	// Create the set storage
	int numInputs = atoi (header.c_str());	

	// The 'outer' GA: Generate truly random parameters for a population of 'inner' GAs
	OGA GAPar(thePars.intGlob["GA_EVO"], thePars.intGlob["GA_SEL"], thePars.intGlob["GA_CHR"], numParams, numInputs, 0, 0, GARange, NNisGA, thePars.boolGlob["IVO"], inFileName, setSeed, thePars);

	// Open output file
	ofstream oGAout(thePars.stringGlob["OGA_DEF"].c_str());
	GAPar.startRun(oGAout);
	
	// Run the 'outer' GA for a number of generations
	short OGAgen = 0;
	while (OGAgen < thePars.intGlob["OGA_TRAIN_ROUNDS"]) {
		GAPar.execute(OGAgen, GARange, numParams, thePars.intGlob["GA_CHR"], thePars);
		GAPar.genOutput(oGAout, OGAgen);
		OGAgen++;
	}

	oGAout.close();
	
}

// Read in a trained NN file and a set of indices, then save each ID with the full set of predictions
void runOnce(const string& NNfile, const string& inputfile, const string& outfile) {

	// Read the index file into a data structure
	ifstream inData(inputfile.c_str());
	if (!inData) {
		cout << "Cannot open index file. " << endl;
		exit (-1);
	}
	vector<string> names;
	vector<string> seqIDs;
	vector<vector<double> > entries;
	readInFile(inData, names, seqIDs, entries);
	
	int numEntries = entries.size();

	// Open the ifstream to read the neural network file
	ifstream inNet(NNfile.c_str());
	if (!inNet) { 
		cout << "Cannot open NN file. " << endl; 
		exit (-1);
	}
	
	// Call the NNfile constructor with the ifstream
	BPNet theNet(inNet, numEntries, names, entries); // Since training is not an issue, use BPNet because it's simpler
	
	// Run the NN with the data
	vector<vector<double> > outNodePredictions;
	outNodePredictions.resize(numEntries);
	theNet.predictSet(outNodePredictions);
	
	// Write the NN predictions to a file
	ofstream outPreds(outfile.c_str());
	writePreds(outPreds, seqIDs, outNodePredictions);

}

// Read the input file for a single NN pass
void readInFile(ifstream& NNfile, vector<string>& names, vector<string>& seqIDs, vector<vector<double> >& entries) {

	string inLine = "";
	vector<string> nameLine;

	int colStart;

	// Get names from the file
	while (inLine.find("SeqID") == string::npos) {
		getline(NNfile, inLine);
	}
	nameLine = split(inLine, ",");
	bool seqFound = false;
	int stopSearch = nameLine.size() - 1;
	for (int i = 0; !seqFound; ++i) {
		if ((nameLine[i]).find("SeqID") != string::npos) { // If the string 'SeqID' is found
			seqFound = true;
			colStart = i;
		} 
		++i;
		if (i == stopSearch) {
			cout << "Header 'SeqID' not found in index file. " << endl;
			exit (-1);
		}
	}

	vector<vector<string> > tableEntries;
	// Read the rest of the file
	while (!NNfile.eof()) {
		getline(NNfile, inLine);
		if (inLine != "") { tableEntries.push_back(split(inLine, ",", colStart + 1)); } // "SeqID" will be the first entry stored
	}
	
	entries.resize(tableEntries.size());
	int lastCol = stopSearch - colStart;
	// Copy everything into the appropriate vector
	for (int j = colStart + 1; j <= stopSearch; ++j) {
		names.push_back(nameLine[j]);
	}
	for (int k = 0; k < tableEntries.size(); ++k) {
		seqIDs.push_back(tableEntries[k][0]);
		(entries[k]).resize(lastCol);
		for (int thisCopy = 1; thisCopy <= lastCol; ++thisCopy) { // Copy the entries from one table to the other
			entries[k][thisCopy - 1] = atof((tableEntries[k][thisCopy]).c_str());
		}
	}
}

// Write the NN predictions to a file
void writePreds(ofstream& outPreds, const vector<string>& seqIDs, const vector<vector<double> >& outNodePredictions) {

	int outSize = seqIDs.size();
	int outNodes = (outNodePredictions[0]).size();
	for (int i = 0; i < outSize; ++i) {
		outPreds << seqIDs[i] << ",";
		for (int j = 0; j < outNodes; ++j) {
			outPreds << outNodePredictions[i][j] << ",";
		}
		outPreds << endl;
	}
}


// Read the command file
void readCommand(const string inCommand, pStr& thePars) {

	setDefaults(thePars);

	ifstream cmdFile(inCommand.c_str());
	string thisLine, globCat, globName;

	if (!cmdFile.is_open()) {
		cout << "Can't open command file.\n";
		exit(-1);
	} 
	// Read a line from the file
	getline(cmdFile, thisLine);
	while (!(cmdFile.eof())) {
		vector<string> thisVec = split(thisLine, " ");
		string vName = thisVec[1];
		string vVal = thisVec[2];
		if (thePars.doubleGlob.find(vName) != thePars.doubleGlob.end()) { // Search for name in the different vectors
			thePars.doubleGlob[vName] = atof (vVal.c_str());
		} else {
			if (thePars.intGlob.find(vName) != thePars.intGlob.end()) {
				thePars.intGlob[vName] = atoi (vVal.c_str());
			} else {
				if (thePars.boolGlob.find(vName) != thePars.boolGlob.end()) {
					thePars.boolGlob[vName] = (bool) (atoi (vVal.c_str()));
				} else {
					if (thePars.stringGlob.find(vName) != thePars.stringGlob.end()) {
						thePars.stringGlob[vName] = vVal;
					} else { // Doesn't go anywhere!
						cout << "Name " << vName << " not found. You may be using an unexpected default value. " << endl;
					}
				}
			}
		}
		getline(cmdFile, thisLine);
	}	
	cmdFile.close(); 
}

// Set default values for the global variables
void setDefaults(pStr& thePars) {

	// Set shared defaults
	thePars.doubleGlob["NHIDNODEA"] = 2.0; thePars.doubleGlob["NHIDNODEB"] = 50.0; 
	thePars.doubleGlob["NOUTNODEA"] = 1.0; thePars.doubleGlob["NOUTNODEB"] = 1.0; // This value may be overridden by the number of categories
	thePars.doubleGlob["NISBIASA"] = 0.00; thePars.doubleGlob["NISBIASB"] = 1.99;

	// GABP
	thePars.doubleGlob["NLRNRATEA"] = 0.0100; thePars.doubleGlob["NLRNRATEB"] = 0.1000;
	thePars.doubleGlob["NMOMENTA"] = 0.0; thePars.doubleGlob["NMOMENTB"] = 0.9;
	thePars.doubleGlob["NWEIGHTDECAYA"] = 0.001; thePars.doubleGlob["NWEIGHTDECAYB"] = 0.05;
	thePars.doubleGlob["NWTSTARTA"] = 0.1; thePars.doubleGlob["NWTSTARTB"] = 1.0;
	thePars.doubleGlob["NLRNDECAYA"] = 0.001; thePars.doubleGlob["NLRNDECAYB"] = 0.100;
	thePars.doubleGlob["NLRNDECAYSTARTA"] = 0.1; thePars.doubleGlob["NLRNDECAYSTARTB"] = 1.0;
	thePars.doubleGlob["NBATCHA"] = 0.0; thePars.doubleGlob["NBATCHB"] = 0.0;
	
	// GAGA
	thePars.doubleGlob["NMUTRATEA"] = 0.1; thePars.doubleGlob["NMUTRATEB"] = 0.9;
	thePars.doubleGlob["NMUTAMTA"] = 1.10; thePars.doubleGlob["NMUTAMTB"] = 1.90;
	thePars.doubleGlob["NMUTPROPA"] = 0.05; thePars.doubleGlob["NMUTPROPB"] = 1.00;
	thePars.doubleGlob["NRECRATEA"] = 0.05; thePars.doubleGlob["NRECRATEB"] = 0.90;
	thePars.doubleGlob["NRECREPLA"] = 0.05; thePars.doubleGlob["NRECREPLB"] = 0.90;
	thePars.doubleGlob["NMIGRATEA"] = 0.5; thePars.doubleGlob["NMIGRATEB"] = 0.5;
	thePars.doubleGlob["NNUMCHRA"] = 10.0; thePars.doubleGlob["NNUMCHRB"] = 100.0;
	thePars.doubleGlob["NNUMEVOA"] = 1.0; thePars.doubleGlob["NNUMEVOB"] = 1.0;
	thePars.doubleGlob["NKILLFRACA"] = 0.01; thePars.doubleGlob["NKILLFRACB"] = 0.90;
	thePars.doubleGlob["NKILLDIFFA"] = 0.05; thePars.doubleGlob["NKILLDIFFB"] = 0.90;
	
	// Score check
	thePars.doubleGlob["LR_TOOLOW"] = 0.001;
	thePars.intGlob["LR_CHECKROUND"] = 50;
	thePars.doubleGlob["TINY_WEIGHT"] = 0.001;
	thePars.intGlob["CHECK_SCORE"] = 20;
	thePars.doubleGlob["NO_SCOREDIF"] = 0.0001;
	
	// OGA
	thePars.doubleGlob["OGA_REC_RATE"] = 0.5;
	thePars.doubleGlob["OGA_REC_REPL"] = 0.5;
	thePars.doubleGlob["OGA_MUT_RATE"] = 0.5;
	thePars.doubleGlob["OGA_MUT_AMT"] = 0.5;
	thePars.doubleGlob["OGA_MUT_PROP"] = 0.5;
	thePars.doubleGlob["OGA_KILL_PROP"] = 0.5;
	thePars.doubleGlob["OGA_KILL_DIFF"] = 0.5;
	thePars.intGlob["GA_CHR"] = 100; 
	thePars.intGlob["GA_EVO"] = 1;
	thePars.intGlob["GA_SEL"] = 1;
	thePars.intGlob["NUM_INPUTS"] = 10; 
	thePars.intGlob["OGA_TRAIN_ROUNDS"] = 30;
	thePars.intGlob["MAX_GA_ROUNDS"] = 500;
	
	// SCORE
	thePars.doubleGlob["WORST_SCORE"] = 0.45;
	thePars.doubleGlob["MIN_GEN"] = 0.9;
	thePars.boolGlob["IVO"] = true;
	thePars.boolGlob["ONODEAVG"] = true;

	// FILENAME
	thePars.stringGlob["NET_DEF"] = "trained.net";
	thePars.stringGlob["OGA_DEF"] = "ogastats.csv";
	
	// NNRUNS
	thePars.intGlob["REPLICATES"] = 1;
	thePars.intGlob["NUM_TO_RECORD"] = 10;
	thePars.intGlob["NN_TRAIN_RUNS"] = 100; 

}
