/*
*  Created by Jan Poleszczuk
*  Last modified November, 2017 by Jan Poleszczuk
*/

#pragma once

#include "stdafx.h"
#include "SIMparameters.h"
#include <algorithm>

class ChemotaxisMap
{

private:

	struct point {//defining a point
		int x;
		int y;
		int z;
	};

	float* map;
	float* kernel;

	int kernelDim;
	int kernelWidth;
	int kernelNum;
	
	char* changes;

	
	unsigned int N1, N2, N3;
	unsigned int numSpots;

	bool is3D;

	SIMparameters* params;



public:

	ChemotaxisMap();
	~ChemotaxisMap();

	void initialize(SIMparameters*);
	void updateMap();

	float getValue(unsigned int pos){ return map[pos]; };

	void addSource(unsigned int);
	void removeSource(unsigned int);

	mxArray* getState();
	
	void printChemoMap() {
		//only 2D
		for (unsigned int i = 0; i < N1; ++i) {
			for (unsigned int j = 0; j < N2; ++j)
				printf("%.2f ", map[i + j*N1]);
			printf("\n");
		}
	}
};
