// DefineVars.h: header file for DefineVars.cpp
// Copyright (C) 2002 Robert G. Beiko

// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

#include <iostream>
#include <vector>
#include <ctype.h>
#include <fstream>

using namespace std;

const int MENU_ITEMS = 8;

const int C_NNRUNS = 0;
const int C_FILENAME = 1;
const int C_SCOREPARAM = 2;
const int C_OGAPARAM = 3;
const int C_STOPCOND = 4;
const int C_NNRANGESHARED = 5;
const int C_NNBPRANGE = 6;
const int C_NNGARANGE = 7;

const int sizes[MENU_ITEMS] = { 3, 2, 4, 13, 5, 6, 14, 20 };

const string C_NAMES[MENU_ITEMS] = { "NNRUNS", "FILENAME", "SCOREPARAM", "OGAPARAM", "STOPCOND", "NNRANGESHARED", "NNBPRANGE", "NNGARANGE" };

// The variable names to be saved
const string NNRUNS[3] = { "NN_TRAIN_RUNS", "REPLICATES", "NUM_TO_RECORD" };
const string FILENAME[2] = { "NET_DEF", "OGA_DEF"};
const string SCOREPARAM[4] = { "WORST_SCORE", "MIN_GEN", "IVO", "ONODEAVG" };
const string OGAPARAM[13] = { "OGA_REC_RATE", "OGA_REC_REPL", "OGA_MUT_RATE", "OGA_MUT_AMT", "OGA_MUT_PROP", "OGA_KILL_PROP", "OGA_KILL_DIFF", "GA_CHR", "GA_EVO", "GA_SEL", "NUM_INPUTS", "OGA_TRAIN_ROUNDS", "MAX_GA_ROUNDS" };
const string STOPCOND[5] = { "LR_TOOLOW", "LR_CHECKROUND", "TINY_WEIGHT", "CHECK_SCORE", "NO_SCOREDIF" };
const string NNRANGESHARED[6] = { "NHIDNODEA","NHIDNODEB","NOUTNODEA","NOUTNODEB","NISBIASA","NISBIASB" };
const string NNBPRANGE[14] = { "NLRNRATEA","NLRNRATEB","NMOMENTA","NMOMENTB","NWEIGHTDECAYA","NWEIGHTDECAYB","NWTSTARTA","NWTSTARTB","NLRNDECAYA","NLRNDECAYB","NLRNDECAYSTARTA","NLRNDECAYSTARTB","NBATCHA","NBATCHB" };
const string NNGARANGE[20] = { "NMUTRATEA","NMUTRATEB","NMUTAMTA","NMUTAMTB","NMUTPROPA","NMUTPROPB","NRECRATEA","NRECRATEB","NRECREPLA","NRECREPLB","NMIGRATEA","NMIGRATEB","NNUMCHRA","NNUMCHRB","NNUMEVOA","NNUMEVOB","NKILLFRACA","NKILLFRACB","NKILLDIFFA","NKILLDIFFB" };

const string* ALLNAMES[8] = { NNRUNS, FILENAME, SCOREPARAM, OGAPARAM, STOPCOND, NNRANGESHARED, NNBPRANGE, NNGARANGE };

// The variable descriptions and one-letter abbrev.
const string NNRUNS_DESC[3] = { "Number of runs to train neural network" , "Number of replicates to use ", "Number of scores to consider in fitness calculation " };

const string FILENAME_DESC[2] = { "Output file for trained neural networks ", "OGA progress file " };

const string SCOREPARAM_DESC[4] = { "Low cutoff of scoring range for selection ", "Minimum generalization score required to save NN ", "Input vector optimization (bool) ", "Average output node values to determine score (useless) "};

const string OGAPARAM_DESC[13] = { "Proportion of chromosomes that recombine ", "Proportion of chromosomes that are replaced by recombinants ", "Mutation rate ", "Magnitude of mutation ", "Percent of loci subject to mutation ",
								   "Proportion of chromosomes to kill due to similarity ", "Maximum difference between chromosomes to be killed ", "OGA Chromosomes ", "OGA Evolvables ", "OGA Selectables ", "Number of inputs per neural network ",
								   "OGA training rounds ", "Maximum GA training rounds " };

const string STOPCOND_DESC[5] = { "For BP, minimum learning rate allowed before training is stopped ", "Round to start checking learning rate decay ", "Minimum connection weights to continue training ", 
				 				  "Use the score from this many previous rounds to compare vs. current score ",
								  "Minimum score improvement required between these two rounds to continue training " };

const string NNRANGESHARED_DESC[6] = { "Hidden nodes MIN ", "Hidden nodes MAX ", "Output nodes MIN ", "Output nodes MAX ", "Bias nodes (bool) MIN ", "Bias nodes (bool) MAX " };

const string NNBPRANGE_DESC[14] = { "Learning rate MIN ", "Learning rate MAX ", "Momentum MIN ", "Momentum MAX ", "Weight decay (bool) MIN ", "Weight decay (bool) MAX ", "Round to begin weight decay MIN ", "Round to begin weight decay MAX ", 
									 "Learning rate decay (bool) MIN ", "Learning rate decay (bool) MAX ", "Round to begin LR decay MIN ", "Round to begin LR decay MAX ", "Batch learning (bool) MIN ", "Batch learning (bool) MAX " };

const string NNGARANGE_DESC[20] = { "Mutation rate MIN ", "Mutation rate MAX ", "Magnitude of mutation MIN ", "Magnitude of mutation MAX ", "Proportion of loci subject to mutation MIN ", "Proportion of loci subject to mutation MAX ",
									"Proportion of IGA chromosomes that recombine MIN ", "Proportion of IGA chromosomes that recombine MAX ", "Proportion of chromosomes that are replaced by recombinants MIN ", "Proportion of chromosomes that are replaced by recombinants MAX ", 
									"Migration rate MIN ", "Migration rate MAX ", "INN Chromosomes MIN ", "INN Chromosomes MAX ", "INN Evolvables MIN ", "INN Evolvables MAX ", 
									"Proportion of IGA chromosomes to kill due to similarity MIN ", "Proportion of IGA chromosomes to kill due to similarity MAX ", "Maximum difference between chromosomes to be killed MIN ", "Maximum difference between chromosomes to be killed MAX "



 };

// The main menu
void mainMenu();
void printMainMenu(char& selection);
void saveVector(const vector<vector<string> >& paramsToSave);
int editVec (const string varNames[], const string desc[], vector<string>& paramEdit);
void initVec (vector<vector<string> >& paramInit);