// classGene.h: header file for classGene.cpp
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

#include"classErr.h"

class GeneBase {

public:

	GeneBase() { geneVal = 0.0; };
	~GeneBase() {};

	virtual void ranStart() = 0;
	virtual void rangeStart(const double& lo, const double& hi) = 0;

	virtual int mutate(double howMuch) = 0;
	virtual int mutate(const double& lo, const double& hi) = 0;

	virtual void forceVal(double toForce) = 0; 

	void copy (const GeneBase& src) { geneVal = src.getVal(); };
	double getVal() const { return geneVal; };

protected:

	double geneVal;
};

// -----------------------------------------

class GeneBool : public GeneBase {

public:

	explicit GeneBool() {};
	~GeneBool() {};

	void ranStart();
	void rangeStart(const double& lo, const double& hi) { this->ranStart(); }	

	int mutate(double howMuch);
	int mutate(const double& lo, const double& hi);

	void forceVal(double toForce) { geneVal = (toForce >= 1.0); };

};

// -----------------------------------------

class GeneInt : public GeneBase {

public:

	explicit GeneInt() {};
	~GeneInt() {};

	void ranStart() { geneVal = rnd(30000.0); };
	void rangeStart(const double& lo, const double& hi) { geneVal = (int) (lo + rnd(hi - lo)); }	

	int mutate(double howMuch);
	int mutate(const double& lo, const double& hi);

	void forceVal(double toForce) { geneVal = (int) toForce; };

};
// -----------------------------------------

class GeneChar : public GeneBase {

public:

	explicit GeneChar() {};
	~GeneChar() {};

	void ranStart() { geneVal = rnd(256); };
	void rangeStart(const double& lo, const double& hi) { geneVal = (char) (lo + rnd(hi-lo)); }	

	int mutate(double howMuch);
	int mutate(const double& lo, const double& hi);

	void forceVal(double toForce) { geneVal = (char) toForce; };

};
// -----------------------------------------

class GeneDouble : public GeneBase {

public:

	explicit GeneDouble() {};
	~GeneDouble() { };

	void ranStart() { geneVal = rnd(30000.0); };
	void rangeStart(const double& lo, const double& hi) { geneVal = lo + rnd (hi - lo); };	

	int mutate(double howMuch);
	int mutate(const double& lo, const double& hi);

	void forceVal(double toForce) { geneVal = toForce; };
};
