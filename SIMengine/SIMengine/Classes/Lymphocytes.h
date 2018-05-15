/*
*  Created by Jan Poleszczuk
*  Last modified November, 2017 by Jan Poleszczuk
*/

#pragma once

#include "stdafx.h"
#include "SIMparameters.h"
#include "Environment.h"
#include "TumorCells.h"

class Lymphocytes
{
private:
	struct Lymphocyte {//defining a tumor cell
		unsigned int place;
		unsigned char p;
		unsigned char kcap;
		unsigned char engaged;
	};

	std::vector<Lymphocyte> cells; //vector containing all cells present in the system
	Environment* env; //pointer to the current environment
	CRandomMersenne* randGen; //pointer to random number generator
	SIMparameters* params; //pointer to parameters
	TumorCells* TUcells;

public:
	size_t NumLymphocytes();

	Lymphocytes(){};
	~Lymphocytes(){};

	void initialize(SIMparameters*, Environment*, CRandomMersenne*, TumorCells*);
	void initializeFromState(const mxArray*, SIMparameters*, Environment*, CRandomMersenne*, TumorCells*);

	void influx();
	void influxInput(const mxArray*);

	void fibrosify();

	void action();

	mxArray* getState();
};

