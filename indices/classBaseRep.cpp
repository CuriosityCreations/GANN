#include"classBaseRep.h"

int degenBaseRep :: countTrue() {

	int yes = 0;

	for (int i = 1; i <= 8; i *= 2) {
		if (theBase & i) { ++yes; }
	}
	return yes;
}

void degenBaseRep :: printComponent(ofstream& outfile) const {

	outfile << BASE_ARRAY[theBase];

}

// ----------------- WEIGHTEDBASEREP ---------------------

// Create from a string in the format 'An Cn Gn Tn' DO NOT SET TO ONE, EVER!!!!
weightedBaseRep :: weightedBaseRep(const string& createFrom) {

	vector<string> eachBase = split(createFrom, " ");
	for (int i = 0; i < 4; ++i) {
		baseScores[i] = atof(eachBase[i].substr(1, eachBase[i].length() - 1).c_str());
	}
}

double weightedBaseRep :: isCompatibleWith(const char baseTest) const { 

	switch (baseTest) {
		case 'A':
		return baseScores[0]; 			
		case 'C':
		return baseScores[1]; 			
		case 'G':
		return baseScores[2];
		case 'T':
		return baseScores[3];
	}
	return -1000000; // Degen bases shouldn't be passed
};


// Split a line of text into a vector of string
std::vector<string> split(const string& source, const string& delims, const int startAt) {
                            
    std::vector<string> destination;   
    int numDelims = delims.length();
    int foundMin = 0;
    int foundStart = 0, foundFinish;
    int totalLength = source.length();
    for (int entry = 1; ; ++entry) { // We'll return out of this loop, so no finishing condition
        foundMin = totalLength; // == string::npos
        for (int i = 0; i < numDelims; ++i) { // Iterate over all delimeters
            int foundThis = (source.substr(foundStart)).find(delims[i]);
            if (foundThis != string::npos && foundThis < foundMin) {
                    foundMin = foundThis;
            }
        }
        foundFinish = foundMin + foundStart; // Minimum element: this is the next delimiter
        if (foundFinish >= totalLength) { // Then found was == string::npos
                string lastPush = source.substr(foundStart, totalLength - foundStart);
                if (lastPush.length() != 0) {
                        destination.push_back(lastPush);
                }
                return destination;
        } else {
            string toPush = source.substr(foundStart, foundFinish - foundStart);
            if (toPush.length() != 0) {
				if (entry >= startAt) {
	                destination.push_back(toPush);
				}
            } else { // We just found a pair of delimiters; don't count the zero-width space between them
            	--entry;
            }
            foundStart = foundFinish + 1;
        }
    }
}

double rnd(double random)
{
	static const double invTwo32 = 1.0 / pow(2.0, 32.0);
	static unsigned long idum = 0;
	idum = 1664525L * idum + 1013904223L;	// From Press et al. 1992
	return idum * random * invTwo32;
};

int firstPosToNumber(const char toChange) {

	char toEmit = toChange;
	switch(toEmit) {
		case 'X':
			return 0;
		case 'A':
			return 1;
		case 'C':
			return 2;
		case 'M':
			return 3;
		case 'G':
			return 4;
		case 'R':
			return 5;
		case 'S':
			return 6;
		case 'V':
			return 7;
		case 'T':
			return 8;
		case 'W':
			return 9;
		case 'Y':
			return 10;
		case 'H':
			return 11;
		case 'K':
			return 12;
		case 'D':
			return 13;
		case 'B':
			return 14;
		case 'N':
			return 15;
	}
	return 0;
}