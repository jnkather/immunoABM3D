/*
*  Created by Jan Poleszczuk
*  Last modified November, 2017 by Jan Poleszczuk
*/

#pragma once

#include "stdafx.h"
#include "SIMparameters.h"
#include "Environment.h"

class Macrophages
{
private:
	struct Macrophage {//defining a tumor cell
		unsigned int place;
		unsigned char p;
		unsigned char state;
	};

	std::vector<Macrophage> cells; //vector containing all cells present in the system
	Environment* env; //pointer to the current environment
	CRandomMersenne* randGen; //pointer to random number generator
	SIMparameters* params; //pointer to parameters

public:
	size_t NumMacrophages();

	Macrophages(){};
	~Macrophages(){};

	void initialize(SIMparameters*, Environment*, CRandomMersenne*);
	void initializeFromState(const mxArray*, SIMparameters*, Environment*, CRandomMersenne*);

	void influx();
	void influxInput(const mxArray*);

	void action();

	mxArray* getState();
};

