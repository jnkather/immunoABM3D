/*
*  Created by Jan Poleszczuk
*  Last modified September, 2017 by Jan Poleszczuk
*/


#include "SIMcore.h"

SIMcore::SIMcore(){};
SIMcore::~SIMcore(){};

void SIMcore::initialize(const mxArray* sysTempl, const mxArray* cnst) {

	//initializing parameters
	params.initialize(sysTempl, cnst);
	//params.listSIMParameters();
    
	//intitializing random number generator
	if (mxIsNaN(params.initSeed)) {//no seed specified
        //mexPrintf("No initial seed specified!\n");
		randGen.RandomInit(time(NULL));
	}
	else {
        //mexPrintf("Initial seed specified!\n");
		randGen.RandomInit((int)params.initSeed);
		randGen.LastInterval = 0;
	}

	//initializing environment (important to be before agents initialization)
	env.initialize(&params, &randGen);

	//initializing cells
	TUcells.initialize(&params, &env, &randGen);
	lymphocytes.initialize(&params, &env, &randGen, &TUcells);
	macrophages.initialize(&params, &env, &randGen);
}

void SIMcore::initializeMexFree() {

	params.listSIMParameters();

	//intitializing random number generator
	randGen.RandomInit((int)params.initSeed);
	randGen.LastInterval = 0;

	//initializing environment
	env.initialize(&params, &randGen);

	//initializing cells
	TUcells.initialize(&params, &env, &randGen);
	lymphocytes.initialize(&params, &env, &randGen, &TUcells);
	macrophages.initialize(&params, &env, &randGen);
}

void SIMcore::initializeFromState(const mxArray* sysTempl, const mxArray* cnst) {

	//initializing parameters
	params.initialize(sysTempl, cnst);
	//params.listSIMParameters();

	//intitializing random number generator
     if (mxIsNaN(params.initSeed)) {//no seed specified
        //mexPrintf("No initial seed specified!\n");
        randGen.RandomInit(time(NULL));
    }
    else {
        //mexPrintf("Initial seed specified!\n");
        randGen.RandomInit((int)params.initSeed);
        randGen.LastInterval = 0;
    }

	//initializing environment (important to be before agents initialization)
	env.initialize(&params, &randGen);
	env.readState(sysTempl);

	//initializing cells
	TUcells.initializeFromState(sysTempl, &params, &env, &randGen);

	lymphocytes.initializeFromState(sysTempl, &params, &env, &randGen, &TUcells);
	
	macrophages.initializeFromState(sysTempl, &params, &env, &randGen);

	//updating maps
	env.updateChemoMap();
	env.updateNecroMap();
}

void SIMcore::readState(const mxArray *state) {
	mxArray *TUc = mxGetField(state, 0, "TUcells");
	TUcells.readState(TUc);
}



size_t SIMcore::NumTUcells() {
	return TUcells.NumTUcells();
}

void SIMcore::getState(mxArray*** theOutput) {
	//this function outputs the inner data to MATLAB

	//1 - creating the most outer structure
	const char *field_names[] = { "env", "TUcells","Lymphocytes","Macrophages" };
	mwSize dims[2] = { 1, 1 };
	(*theOutput)[0] = mxCreateStructArray(1, dims, 4, field_names);

	//exporting environment
	mxSetField((*theOutput)[0],0,"env",env.getState());

	//exporting TUcells
	mxSetField((*theOutput)[0],0,"TUcells",TUcells.getState());

	//exporting TUcells
	mxSetField((*theOutput)[0], 0, "Lymphocytes", lymphocytes.getState());

	//exporting TUcells
	mxSetField((*theOutput)[0], 0, "Macrophages", macrophages.getState());
	
}

void SIMcore::TUcellsNum(mxArray*** theOutput) {
	//this function outputs the inner data to MATLAB
	mwSize dims[3] = { 1, 1, 1 };
	//1 - creating the most outer structure
	(*theOutput)[0] =  mxCreateNumericArray(1, dims, mxINT32_CLASS, mxREAL);
	int *ptr = (int*)mxGetPr((*theOutput)[0]);
	ptr[0] = (int)TUcells.NumTUcells();		
}
