/*
*  Created by Jan Poleszczuk
*  Last modified November, 2017 by Jan Poleszczuk
*/

#pragma once

#include "stdafx.h"
#include "SIMparameters.h"
#include "ChemotaxisMapSuiteS.h"
//#include "NecrosisMap.h"
#include "NecrosisMapSuiteS.h"

class Environment
{

private:

	struct point{//defining a point
		int x;
		int y;
		int z;
	};

	unsigned int *indcNeigh;//neighborhood
	bool *lattice;
	bool *latticeOnlyTU;
	bool *necrosis;
	bool *fibrosis;

	unsigned int numOccupiedSpots;
	unsigned int numNecroticSpots;
	unsigned int numFibroticSpots;

	float *AdjuMap;

	ChemotaxisMap ChtaxMap;
	NecrosisMap NecroMap;
	

	unsigned int N1, N2, N3;
	unsigned int numSpots;

	CRandomMersenne* randGen; //pointer to random number generator
	SIMparameters* params; //pointer to parameters

	bool is3D;
	int numNeigh;
	unsigned int* neigh;
	float* chemoGrad;

	std::vector<point> adjuMask, smoothMask;

public:
	unsigned int returnEmptyPlace(unsigned int);
	unsigned int returnEmptyPlaceImmune(unsigned int, float);
	
	unsigned int findTarget(unsigned int);

	float getAdjuValue(unsigned int pos){ return ChtaxMap.getValue(pos); };
	double getNecrosisValue(unsigned int pos){ return NecroMap.getValue(pos); };

	Environment();
	~Environment();
	void initialize(SIMparameters*, CRandomMersenne*);
	void readState(const mxArray*);

	unsigned int getCenter();
	void newTUcell(unsigned int);
	void deleteTUcell(unsigned int);
	void setOccupancy(unsigned int, bool);
	void markNecrosis(unsigned int pos){ necrosis[pos] = true; };

	void modulateAdjuvanticity(std::vector<unsigned int>::iterator, std::vector<unsigned int>::iterator, float);
	void modulateFibrosis(std::vector<unsigned int>::iterator, std::vector<unsigned int>::iterator);
	void decayAdjuvanicityMap();

	void updateNecroMap() { NecroMap.updateMap(); };
	void updateChemoMap() { ChtaxMap.updateMap(); };
    

	unsigned int* generatePostions(int*);

	mxArray* getState();
};
