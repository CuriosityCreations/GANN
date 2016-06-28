// classGene.cpp defines various types of Gene objects used by Chromosomes in the genetic algorithm.
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

#include"classGene.h"

// Class GeneBool -----------------------------
void GeneBool :: ranStart() { 
	geneVal = 0.0;
	if (maybe()) {
		if (maybe()) {
			geneVal = 1.0; 
		}
	}		
}

// DO NOT MUTATE BOOL GENES!!!
int GeneBool :: mutate(double howMuch) {
	return (int) geneVal;
}

int GeneBool :: mutate(const double&lo, const double& hi) {
	return mutate(0.0);
}

// Class GeneInt ------------------------------
int GeneInt :: mutate(double howMuch) {

	if (maybe()) {
                howMuch = 1.0 / howMuch;
                geneVal *= howMuch + rnd(1.0 - howMuch);
        } else {
                geneVal *= 1.0 + rnd(howMuch - 1.0);
        }
	return 0;
}

int GeneInt :: mutate(const double& lo, const double& hi) {
	geneVal = lo + rnd(hi - lo);
	return 0;
}

// Class GeneChar -----------------------------
int GeneChar :: mutate(double howMuch) { this->ranStart(); return 0; }

int GeneChar :: mutate(const double& lo, const double& hi) {
	geneVal = lo + rnd(hi - lo);
	return 0;
}

// Class GeneDouble ---------------------------

int GeneDouble :: mutate(double howMuch) {

	if (maybe()) {
                howMuch = 1.0 / howMuch;
                geneVal *= howMuch + rnd(1.0 - howMuch);
        } else {
                geneVal *= 1.0 + rnd(howMuch - 1.0);
        }
	return 0;
}

int GeneDouble :: mutate(const double& lo, const double& hi) {

	geneVal = lo + rnd(hi - lo);
	return 0;
}
