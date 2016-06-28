// classGene.h: header file for classGene.cpp
// GANN CORE
// Copyright (C) 2000-2002 Robert G. Beiko

// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

#include"classBaseRep.h"

class GeneBase {

public:

	virtual ~GeneBase() {};
	
	virtual GeneBase* clone() = 0; // Create a copy of the original

	virtual bool isCompatible(const string& seq, const int startPt, 
		const int lengthToCheck, vector<map<string, double> >& ruleVec, 
		map<string, int>& nameToNum) const = 0; // Is the motif found at a given point in the sequence?
	virtual void fixIfIllegal() { }; // This is for motifGenes
	virtual void fixIfIllegal(const int, const int) { }; // This is for posGenes

	virtual int itsMinLength() const = 0;
	virtual int itsMaxLength() const = 0;

	virtual void printComponent(ofstream& outfile) = 0;

	virtual double itsThreshold() const { return -1; };
};

// -----------------------------------------
// Gene that identifies the first and final permissible starting positions of a motif, scaled between 0 and seqlen
class posGene : public GeneBase {

public:

	explicit posGene(const int seqLen);
	posGene(const int start, const int stop) { firstPos = start; lastPos = stop; }
	posGene(const posGene& rhs) { firstPos = rhs.firstPos; lastPos = rhs.lastPos; };
	virtual ~posGene() {};

	virtual posGene* clone() { return new posGene(*this); };

	void fixIfIllegal(const int maxMotif, const int seqLen); // Fix firstPos and lastPos if motif is too large
	virtual bool isCompatible(const string&, const int, 
		const int, vector<map<string, double> >&, 
		map<string, int>&) const { return true; }; // Is the motif found at a given point in the sequence?
	virtual void printComponent(ofstream& outfile);

	virtual int itsMinLength() const { return lastPos; };
	virtual int itsMaxLength() const { return firstPos; };
	
private:

	int firstPos, lastPos;

};

// -----------------------------------------

class motifGene : public GeneBase {

public:

	virtual ~motifGene() {};
	
	virtual motifGene* clone() = 0;

	virtual bool isCompatible(const string& seq, const int startPt, 
		const int lengthToCheck, vector<map<string, double> >& ruleVec, 
		map<string, int>& nameToNum) const = 0; // Is the motif found at a given point in the sequence?
	virtual void fixIfIllegal() { };
	virtual void printComponent(ofstream& outfile) = 0;

	virtual int itsMinLength() const { return minLength; };
	virtual int itsMaxLength() const { return maxLength; };

protected:

	void setLength(const int minLeng, const int maxLeng) { minLength = minLeng; maxLength = maxLeng; };

	int minLength; // Number of bases covered by this motif component
	int maxLength;

};

// -----------------------------------------

class discreteGene : public motifGene {

public:

	discreteGene(const string& createFrom); // Create from string rep
//	discreteGene(const discreteGene& rhs);
	virtual ~discreteGene();

	virtual discreteGene* clone() { return new discreteGene(*this); };

	virtual bool isCompatible(const string& seq, const int startPt, 
		const int lengthToCheck, vector<map<string, double> >& ruleVec, 
		map<string, int>& nameToNum) const;
	virtual void fixIfIllegal();
	virtual void printComponent(ofstream& outfile);

	virtual double itsThreshold() const { return threshold; };

private:

	vector<baseRep*> baseSeq;
	double threshold; // Threshold score for accepting isCompatible()

};

// -----------------------------------------

class continuousGene : public motifGene {

public:

	continuousGene() { };
	continuousGene(const string& createFrom);

	continuousGene(const continuousGene& rhs) { minLength = rhs.minLength; maxLength = rhs.maxLength; };
	virtual ~continuousGene() { };
	
	virtual continuousGene* clone() { return new continuousGene(*this); };

	virtual bool isCompatible(const string&, const int, 
		const int, vector<map<string, double> >&, 
		map<string, int>&) const { 
//		cout << "C" << flush;
		return true; 
		
		}; // Spacer is always compatible
	virtual void fixIfIllegal();
	virtual void printComponent(ofstream& outfile) { outfile << "N" << minLength << "-" << maxLength << "|"; };
};

// ------------------------------------------

class rangedContinuousGene : public continuousGene {

public:

	rangedContinuousGene(const string& createFrom, vector<map<string, double> >& ruleVec, vector<string>& ruleNames);
//	rangedContinuousGene(const rangedContinuousGene& rhs) { minLength = rhs.minLength; maxLength = rhs.maxLength; ruleType = rhs.ruleType; ruleName = rhs.ruleName; ruleSize = rhs.ruleSize; minVal = rhs.minVal; maxVal = rhs.maxVal; };
	virtual ~rangedContinuousGene() {};
	
	virtual rangedContinuousGene* clone() { return new rangedContinuousGene(*this); };

	virtual bool isCompatible(const string& seq, const int startPt, 
		const int lengthToCheck, vector<map<string, double> >& ruleVec, 
		map<string, int>& nameToNum) const;
	virtual void fixIfIllegal();
	virtual void printComponent(ofstream& outfile);

private:

	int ruleType;
	string ruleName;
	int ruleSize; // The number of bases specified in the rule
	double minVal, maxVal; // The range of the continuous character
	double mutateVal;
	
};
