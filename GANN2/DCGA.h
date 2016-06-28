// DCGA.h defines constants and global variables for use by GANN
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

#ifndef DCGAdef
#define DCGAdef

using namespace std;

#include <iostream>
#include <fstream>
#include <time.h>
#include <cmath>

#include <string>
#include <functional>
#include <vector>
#include <algorithm>
#include <map>

// #pragma warning (disable : 4244)

const double OC_SENTRY = -100.0;

const string U_MESSAGE = "Invalid argument. Usage: ";
const string U_FLAGS = " ([-b | -g] indexFile commandFile OR -r trainednet indexfile outfile)";

const double PI = 3.14159265358979;

// The e^x shortcut table
const int ETABLESIZE = 4096; // Used to be 8192
const int EARRAYSIZE = 2 * ETABLESIZE - 1;
const int EMINUSONE = ETABLESIZE - 1;
const int NEGEMONE = EMINUSONE * -1;
const double EDIVTEN = ETABLESIZE / 10.0;
const int LOGET = (int)( log((double) ETABLESIZE) / log(2.0));
const int INTSHIFT = sizeof(int) * 8 - 1;

// Shared plus specific
const short N_BP_PARAMS = 10;
const short N_GA_PARAMS = 13;

//////////// Shared OGA parameters
const int XHIDNODE = 0, XOUTNODE = 1, XISBIAS = 2;

//////////////// GABP specific
const int XLRNRATE = 3, XMOMENT = 4, XWEIGHTDECAY = 5, XWTSTART = 6, XLRNDECAY = 7, 
			 XLRNDECAYSTART = 8, XBATCH = 9;

//////////////// GAGA specific
const int XMUTRATE = 3, XMUTAMT = 4, XMUTPROP = 5, XRECRATE = 6, XRECREPL = 7, XMIGRATE = 8, 
			 XNUMCHR = 9, XNUMEVO = 10, XKILLFRAC = 11, XKILLDIFF = 12;

//int TOTAL_BPSIZE; 
// = N_BP_PARAMS + NUM_INPUTS;
//int TOTAL_GASIZE; 
// = N_GA_PARAMS + NUM_INPUTS;

const vector<double> BP_PARAMETERS;
const vector<double> GA_PARAMETERS;

const string BP_PARAM_RANGES[N_BP_PARAMS * 2] = {"NHIDNODEA","NHIDNODEB","NOUTNODEA","NOUTNODEB","NISBIASA","NISBIASB","NLRNRATEA","NLRNRATEB","NMOMENTA","NMOMENTB","NWEIGHTDECAYA","NWEIGHTDECAYB","NWTSTARTA", "NWTSTARTB","NLRNDECAYA","NLRNDECAYB","NLRNDECAYSTARTA","NLRNDECAYSTARTB","NBATCHA","NBATCHB" };
const string GA_PARAM_RANGES[N_GA_PARAMS * 2] = { "NHIDNODEA","NHIDNODEB","NOUTNODEA","NOUTNODEB","NISBIASA","NISBIASB","NMUTRATEA","NMUTRATEB","NMUTAMTA","NMUTAMTB","NMUTPROPA","NMUTPROPB","NRECRATEA","NRECRATEB","NRECREPLA","NRECREPLB","NMIGRATEA","NMIGRATEB","NNUMCHRA","NNUMCHRB","NNUMEVOA","NNUMEVOB","NKILLFRACA","NKILLFRACB","NKILLDIFFA","NKILLDIFFB" };

const string BP_PARAM_NAMES[N_BP_PARAMS] = { "hidNodes","outNodes","biasYN","lrnRate","momentum","weightDecay","wtDecayStart","lrnDecay","lrnDecayStart","BatchYN" };
const string GA_PARAM_NAMES[N_GA_PARAMS] = { "hidNodes","outNodes","biasYN","mutRate","mutAmt","mutProp","recRate","recRepl","migRate","numChr","numEvo","killFrac","killDiff" };

const string OGADEFAULT = "ogastats.csv";

const string OCOUTFILES[5] = { "nntrainfpscores.csv",
					  "nntraindistscores.csv",
					  "nntestfpscores.csv",
					  "nntestdistscores.csv",
					  "nnparams.csv" };

// For scoring --------------------------------------
enum { NO_LEARN = 1, WT_CRASH = 2, LR_CRASH = 3 };

enum { TRAINFPSCORE = 0, TRAINDIFFSCORE = 1, TESTFPSCORE = 2, TESTDIFFSCORE = 3 };


struct pars {
	map<string,double> doubleGlob;
	map<string,int> intGlob;
	map<string,bool> boolGlob;
	map<string,string> stringGlob;
};

typedef struct pars pStr;

// NNRUNS ---------------------------------------------------------------------

const int C_NNRUNS = 3;

//int NN_TRAIN_RUNS; 
//int REPLICATES;
//int NUM_TO_RECORD;

// FILENAME --------------------------------------------------------------------

const int C_FILENAME = 2;

//string NET_DEF;
//string OGA_DEF;

// SCOREPARAM ------------------------------------------------------------------

const int C_SCOREPARAM = 4;

//double WORST_SCORE;
//double MIN_GEN;
//bool IVO;
//bool ONODEAVG;

// OGAPARAM ---------------------------------------------------------------------

const int C_OGAPARAM = 13;

//double OGA_REC_RATE;
//double OGA_REC_REPL;
//double OGA_MUT_RATE;
//double OGA_MUT_AMT;
//double OGA_MUT_PROP;
//double OGA_KILL_PROP;
//double OGA_KILL_DIFF;
//int GA_CHR; 
//int GA_EVO;
//int GA_SEL;
//int NUM_INPUTS; 
//int OGA_TRAIN_RUNS;
//int MAX_GA_ROUNDS;

// STOPCOND ----------------------------------------------------------------------

const int C_STOPCOND = 5;

//double LR_TOOLOW;
//int LR_CHECKROUND;
//double TINY_WEIGHT;
//int CHECK_SCORE;
//double NO_SCOREDIF;

// NNRANGE -------------------------------------------------------------------------

const int C_NNRANGE = 40;
/*									
// Shared OGA
double NHIDNODEA; double NHIDNODEB;
double NOUTNODEA; double NOUTNODEB;  
double NISBIASA; double NISBIASB;

// GABP

double NLRNRATEA; double NLRNRATEB;
double NMOMENTA; double NMOMENTB;
double NWEIGHTDECAYA; double NWEIGHTDECAYB;
double NWTSTARTA; double NWTSTARTB;
double NLRNDECAYA; double NLRNDECAYB;
double NLRNDECAYSTARTA; double NLRNDECAYSTARTB;
double NBATCHA; double NBATCHB;

// GAGA
double NMUTRATEA; double NMUTRATEB;
double NMUTAMTA; double NMUTAMTB;
double NMUTPROPA; double NMUTPROPB;
double NRECRATEA; double NRECRATEB;
double NRECREPLA; double NRECREPLB;
double NMIGRATEA; double NMIGRATEB;
double NNUMCHRA; double NNUMCHRB;
double NNUMEVOA; double NNUMEVOB;
double NKILLFRACA; double NKILLFRACB;
double NKILLDIFFA; double NKILLDIFFB;
*/
// The number of variables that must be read from an external file
const int EXTERN_VARS = C_NNRUNS + C_FILENAME + C_SCOREPARAM + C_NNRANGE + C_OGAPARAM + C_STOPCOND; 

// General function declarations
void setGAVec(std::vector<double>& paramVec, const bool NNisGA, const pStr& thePars); 

// Split a line of text into a vector, using any of the provided delimeters to split text into vector elements
std::vector<string> split(const string& source, const string& delims, const int startAt = 1);

// Get map element in a const-friendly way (unlike [])
template <class T> T elem(const map<string, T>& mapFind, const string searchParam) { return (*(mapFind).find(searchParam)).second; }

double rnd(double random);
	// Returns a random double between 0.0 and <random>, with 32-bit resolution
	// <random> may be positive or negative.

inline unsigned long rnd(unsigned long random) {
	unsigned long myResult = long(rnd(double(random)));
	return (myResult < random) ? myResult : myResult - 1;
};

inline unsigned short rnd(unsigned short random) {
	unsigned short myResult = short(rnd(double(random)));
	return (myResult < random) ? myResult : myResult - (unsigned short) 1;
};

inline unsigned int rnd(unsigned int random) {
	unsigned int myResult = int(rnd(double(random)));
	return (myResult < random) ? myResult : (unsigned int) myResult - 1;
};

inline long rnd(long random) {
	long myResult = long(rnd(double(random)));
	return (myResult == random) ? myResult - 1 + 2 * long(random < 0) : myResult;
};

inline short rnd(short random) {
	short myResult = short(rnd(double(random)));
	return (myResult == random) ? myResult - 1 + 2 * short(random < 0) : myResult;
};

inline int rnd(int random) {
	int myResult = int(rnd(double(random)));
	return (myResult == random) ? myResult - 1 + 2 * int(random < 0) : myResult;
};
	// The above return a random integer between 0 and <random> - 1, with 32-bit resolution
	// (For negative <random>, return value is between <random> + 1 and 0.)
	// They call: double rnd(double random)

inline bool maybe() { return (rnd(1.0) < 0.5); };

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)

// ran1 function declaration
double ran1(long *);

#endif

