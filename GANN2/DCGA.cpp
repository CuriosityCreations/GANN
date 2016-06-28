// DCGA.cpp defines a couple of useful functions used by the GANN core.
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

#include "DCGA.h"

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
                        }
                        foundStart = foundFinish + 1;
                }
        }
} 

void setGAVec(std::vector<double>& paramVec, const bool isGA, const pStr& thePars) {

	if (isGA) {
        for (int i = 0; i < N_GA_PARAMS * 2; ++i)
                paramVec[i] = elem(thePars.doubleGlob,GA_PARAM_RANGES[i]);
	} else {
		for (int i = 0; i < N_BP_PARAMS * 2; ++i)
                paramVec[i] = elem(thePars.doubleGlob,BP_PARAM_RANGES[i]); 
	}
}

// From Numerical Recipes in C
double ran1(long *idum) {

	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <=0 || !iy) {
		if (-(*idum) < 1) *idum = 1;
		else *idum = -(*idum);
		for (j=NTAB+7; j>=0; j--) {
			k = (*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k= (*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j = iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

double rnd(double random)
{
	static const double invTwo32 = 1.0 / pow(2.0, 32.0);
	static unsigned long idum = 0;
	idum = 1664525L * idum + 1013904223L;	// From Press et al. 1992
	return idum * random * invTwo32;
};
