/*
*  Created by Jan Poleszczuk
*  Last modified November, 2017 by Jan Poleszczuk
*/

#pragma once

#include "stdafx.h"
#include "SIMparameters.h"
#include "Environment.h"

class TumorCells
{

public:
	struct TumorCell {//defining a tumor cell
		unsigned int place;
		unsigned char p;
		bool is_stem;
		float Antigen;
		unsigned char damage;
	};

	void go_grow_die();
	void modulateAdjuvanticity();

	void getStateForImmune(std::map<unsigned int, std::vector<TumorCells::TumorCell>::iterator>*);

	size_t NumTUcells();

	TumorCells(){};
	~TumorCells(){};

	void initialize(SIMparameters*, Environment*, CRandomMersenne*);
	void initializeFromState(const mxArray*, SIMparameters*, Environment*, CRandomMersenne*);
	void readState(mxArray*);

	mxArray* getState();

private:
	

	std::vector<TumorCell> cells; //vector containing all cells present in the system
	std::vector<unsigned int> dyingCells; //vector containing information about dying cells in previous iteration

	Environment* env; //pointer to the current environment
	CRandomMersenne* randGen; //pointer to random number generator
	SIMparameters* params; //pointer to parameters

};

