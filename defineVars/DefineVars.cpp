// DefineVars.cpp creates a list of variable values for input to the GANN core software.
// Copyright (C) 2002 Robert G. Beiko

// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

#include"DefineVars.h"

// A particularly useless main()
int main() {

	mainMenu();
	return 0;
}

void mainMenu() {

	// Create and size the vector of vectors
	vector<vector<string> > SizeVec(MENU_ITEMS);
	for (int i = 0; i < MENU_ITEMS; ++i) {
		(SizeVec[i]).resize(sizes[i]);
	}
	
	// Initialize vector with default parameters
	initVec(SizeVec);

/*	// Determine the start points within the vector for each set of variables
	int cumulCount = 0;
	int varBase[8] = { 0 };
	for (int j = 0; j <= 8; ++j) {
		varBase[j] = cumulCount;
		cumulCount += sizes[j];
	} */

	char selection = 'X';

	while (toupper(selection) != 'Q') {

		printMainMenu(selection);
	
		switch (toupper(selection)) {
			case 'N':
				editVec(NNRUNS, NNRUNS_DESC, SizeVec[C_NNRUNS]);
				break;
			case 'F':
				editVec(FILENAME, FILENAME_DESC, SizeVec[C_FILENAME]);
				break;
			case 'S':
				editVec(SCOREPARAM, SCOREPARAM_DESC, SizeVec[C_SCOREPARAM]);
				break;
			case 'O':
				editVec(OGAPARAM, OGAPARAM_DESC, SizeVec[C_OGAPARAM]);
				break;
			case 'T':
				editVec(STOPCOND, STOPCOND_DESC, SizeVec[C_STOPCOND]);
				break;
			case 'R':
				editVec(NNRANGESHARED, NNRANGESHARED_DESC, SizeVec[C_NNRANGESHARED]);
				break;
			case 'B':
				editVec(NNBPRANGE, NNBPRANGE_DESC, SizeVec[C_NNBPRANGE]);
				break;
			case 'G':
				editVec(NNGARANGE, NNGARANGE_DESC, SizeVec[C_NNGARANGE]);
				break;
			case 'A':
				saveVector(SizeVec);
				break;
			case 'Q':
				break;
			default:
				cout << "Invalid character. Try again." << endl;
		}
	}
}

void printMainMenu(char& selection) {

	cout << "*** Main Menu \n\n" << endl;
	cout << "Please select which category of variables to modify:\n" << endl;
	cout << "<N>eural network run parameters\n" << endl;
	cout << "<F>ile names\n" << endl;
	cout << "<S>coring parameters\n" << endl;
	cout << "<O>GA parameters\n" << endl;
	cout << "S<T>op conditions\n" << endl;
	cout << "NN shared <R>ange parameters\n" << endl;
	cout << "<B>P parameters\n" << endl;
	cout << "<G>A parameters\n" << endl;
	cout << "S<A>ve my settings\n" << endl;
	cout << "<Q>uit this program\n" << endl;
	cout << "Well? ";
	cin >> selection;
	cin.ignore(100,'\n');
}

// Edit the vector of parameter values. The first three parameters all have the same cardinality as the fourth
int editVec (const string varNames[], const string desc[], vector<string>& paramEdit) {

	int numQuantities = paramEdit.size(); // Number of entries in the vector
	int selInt;
	string sel;

	while (1) {
		cout << "Choose which quantity to modify: \n" << endl;
		for (int i = 0; i < numQuantities; ++i) {
			cout << i << "  " << desc[i] << "   (" << paramEdit[i] << ")" << endl;
		}
		cout << "99  Return to main menu" << endl;
		cout << "\nEnter your selection: ";
		getline(cin, sel, '\n');
		selInt = atoi(sel.c_str());
		if (selInt == 99) { return 1; }
		if (selInt >= paramEdit.size()) { cout << "Invalid selection.\n"; }
		else { 
			cout << "Enter new value for variable " << varNames[selInt] << ": ";
			getline(cin, sel, '\n');
			paramEdit[selInt] = sel; 
		}
	}
}

// This is probably the ugliest function I have ever written. It sets the default values for the entire parameter set
void initVec (vector<vector<string> >& paramInit) {

	paramInit[C_NNRUNS][0] = "100"; // NN_TRAIN_RUNS, 
	paramInit[C_NNRUNS][1] = "1"; // REPLICATES, 
	paramInit[C_NNRUNS][2] = "10"; // NUM_TO_RECORD 

	paramInit[C_FILENAME][0] = "trained.net"; // NET_DEF
	paramInit[C_FILENAME][1] = "ogastats.csv"; // OGA_DEF
	
	paramInit[C_SCOREPARAM][0] = "0.45"; // WORST_SCORE
	paramInit[C_SCOREPARAM][1] = "0.9"; // MIN_GEN
	paramInit[C_SCOREPARAM][2] = "1.99"; // IVO
	paramInit[C_SCOREPARAM][3] = "1.99"; // ONODEAVG

	paramInit[C_OGAPARAM][0] = "0.3"; // OGA_REC_RATE
	paramInit[C_OGAPARAM][1] = "0.3"; // OGA_REC_REPL
	paramInit[C_OGAPARAM][2] = "0.7"; // OGA_MUT_RATE
	paramInit[C_OGAPARAM][3] = "1.1"; // OGA_MUT_AMT
	paramInit[C_OGAPARAM][4] = "100"; // OGA_MUT_PCT
	paramInit[C_OGAPARAM][5] = "0.1"; // OGA_KILL_PROP
	paramInit[C_OGAPARAM][6] = "0.1"; // OGA_KILL_DIFF
	paramInit[C_OGAPARAM][7] = "150"; // GA_CHR
	paramInit[C_OGAPARAM][8] = "1"; // GA_EVO
	paramInit[C_OGAPARAM][9] = "1"; // GA_SEL
	paramInit[C_OGAPARAM][10] = "10"; // NUM_INPUTS
	paramInit[C_OGAPARAM][11] = "30"; // OGA_TRAIN_ROUNDS
	paramInit[C_OGAPARAM][12] = "500"; // MAX_GA_ROUNDS

	paramInit[C_STOPCOND][0] = "0.001"; // LR_TOOLOW
	paramInit[C_STOPCOND][1] = "50"; // LR_CHECKROUND
	paramInit[C_STOPCOND][2] = "0.001"; // TINY_WEIGHT
	paramInit[C_STOPCOND][3] = "20"; // CHECK_SCORE
	paramInit[C_STOPCOND][4] = "0.0001"; // NO_SCOREDIF

	paramInit[C_NNRANGESHARED][0] = "2.0"; // NHIDNODEA
	paramInit[C_NNRANGESHARED][1] = "50.0"; // NHIDNODEB
	paramInit[C_NNRANGESHARED][2] = "1.0"; // NOUTNODEA
	paramInit[C_NNRANGESHARED][3] = "1.0"; // NOUTNODEB
	paramInit[C_NNRANGESHARED][4] = "0.00"; // NISBIASA
	paramInit[C_NNRANGESHARED][5] = "1.99"; // NISBIASB

	paramInit[C_NNBPRANGE][0] = "0.01"; // NLRNRATEA
	paramInit[C_NNBPRANGE][1] = "0.1"; // NLRNRATEB
	paramInit[C_NNBPRANGE][2] = "0.0"; // NMOMENTA
	paramInit[C_NNBPRANGE][3] = "0.9"; // NMOMENTB
	paramInit[C_NNBPRANGE][4] = "0.001"; // NWEIGHTDECAYA
	paramInit[C_NNBPRANGE][5] = "0.05"; // NWEIGHTDECAYB
	paramInit[C_NNBPRANGE][6] = "0.1"; // NWTSTARTA
	paramInit[C_NNBPRANGE][7] = "1.0"; // NWTSTARTB
	paramInit[C_NNBPRANGE][8] = "0.001"; // NLRNDECAYA
	paramInit[C_NNBPRANGE][9] = "0.1"; // NLRNDECAYB
	paramInit[C_NNBPRANGE][10] = "0.1"; // NLRNDECAYSTARTA
	paramInit[C_NNBPRANGE][11] = "1.0"; // NLRNDECAYSTARTB
	paramInit[C_NNBPRANGE][12] = "0.0"; // NBATCHA
	paramInit[C_NNBPRANGE][13] = "0.0"; // NBATCHB
	
	paramInit[C_NNGARANGE][0] = "0.1"; // NMUTRATEA
	paramInit[C_NNGARANGE][1] = "0.9"; // NMUTRATEB
	paramInit[C_NNGARANGE][2] = "1.1"; // NMUTAMTA
	paramInit[C_NNGARANGE][3] = "1.9"; // NMUTAMTB
	paramInit[C_NNGARANGE][4] = "0.05"; // NMUTPROPA
	paramInit[C_NNGARANGE][5] = "1.00"; // NMUTPROPB
	paramInit[C_NNGARANGE][6] = "0.05"; // NRECRATEA
	paramInit[C_NNGARANGE][7] = "0.9"; // NRECRATEB
	paramInit[C_NNGARANGE][8] = "0.05"; // NRECREPLA
	paramInit[C_NNGARANGE][9] = "0.9"; // NRECREPLB
	paramInit[C_NNGARANGE][10] = "0.5"; // NMIGRATEA
	paramInit[C_NNGARANGE][11] = "0.5"; // NMIGRATEB
	paramInit[C_NNGARANGE][12] = "10"; // NNUMCHRA
	paramInit[C_NNGARANGE][13] = "100"; // NNUMCHRB
	paramInit[C_NNGARANGE][14] = "1.0"; // NNUMEVOA
	paramInit[C_NNGARANGE][15] = "1.0"; // NNUMEVOB
	paramInit[C_NNGARANGE][16] = "0.01"; // NNKILLFRACA
	paramInit[C_NNGARANGE][17] = "0.9"; // NNKILLFRACB
	paramInit[C_NNGARANGE][18] = "0.05"; // NNKILLDIFFA
	paramInit[C_NNGARANGE][19] = "0.9"; // NNKILLDIFFB

}

void saveVector(const vector<vector<string> >& paramsToSave) {

	cout << "Enter command file name: " << endl;
	string outName;
	getline(cin, outName);
	
	ofstream outF(outName.c_str());
	
	for (int i = 0; i < MENU_ITEMS; ++i) {
		for (int j = 0; j < sizes[i]; ++j) {
			outF << C_NAMES[i] << " " << ALLNAMES[i][j] << " " << paramsToSave[i][j] << endl;
		}
	}
}
