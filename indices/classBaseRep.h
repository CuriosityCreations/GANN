#include<iostream>
#include<map>
#include<string>
#include<vector>
#include<fstream>
// #include<strstream>
#include<math.h>
#include<sstream>

// class degenBaseRep represents a single base as a combination of A, C, G, and T.

using namespace std;

int firstPosToNumber(const char toChange);

const char BASE_ARRAY[16] = { 'X', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N' };

class baseRep { // Abstract base class

public:
	
	virtual baseRep* clone() const = 0;

	virtual void printComponent(ofstream& outfile) const = 0;

	virtual double isCompatibleWith(const char baseTest) const = 0;

private:

};

class degenBaseRep : public baseRep {

public:
	
	degenBaseRep(const string& createFrom) { theBase = firstPosToNumber(createFrom[0]); countTrue(); };

	degenBaseRep(const degenBaseRep& rhs) { theBase = rhs.theBase; numTrue = rhs.numTrue; };	

	~degenBaseRep() { };
	
	virtual degenBaseRep* clone() const { return new degenBaseRep(*this); };

	int countTrue(); // How many of the four bases are represented?

	virtual double isCompatibleWith(const char baseTest) const { return (!(~theBase & firstPosToNumber(baseTest))) ? 1.0 : 0.0; }; // See whether baseTest is a legal member of the Rep
//	virtual double isCompatibleWith(const char baseTest) const { return (!(~theBase & BASE_ARRAY[baseTest])) ? 1.0 : 0.0; }; // See whether baseTest is a legal member of the Rep
		
	virtual void printComponent(ofstream& outfile) const;

private:

	int theBase;
	int numTrue;

};

class weightedBaseRep : public baseRep {
public:

	weightedBaseRep(const string& createFrom);

//	weightedBaseRep(const weightedBaseRep& rhs) { for (int i = 0; i < 4; ++i) baseScores[i] = rhs.itsBaseScore(i); };

	~weightedBaseRep() { };
	virtual weightedBaseRep* clone() const { return new weightedBaseRep(*this); };

	virtual double isCompatibleWith(const char baseTest) const;

	virtual void printComponent(ofstream& outfile) const { outfile << "A" << baseScores[0] << " C" << baseScores[1] << " G" << baseScores[2] << " T" << baseScores[3]; };

private:

	float baseScores[4]; // Scores associated with A, C, G, T

};

typedef vector<baseRep*>::const_iterator baseIt;

double rnd(double random);

inline void die(const string& message) { cout << message << "\n"; exit(-1); }

std::vector<string> split(const string& source, const string& delims, const int startAt = 1);

