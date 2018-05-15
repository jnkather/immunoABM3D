/*
 *  Created by Jan Poleszczuk
 *  Last modified September, 2017 by Jan Poleszczuk
 */


#include "ClassHandle.h"

#include "SIMengine/Classes/SIMcore.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Get the command string
	char cmd[64];
	if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
		mexErrMsgTxt("First input should be a command string less than 64 characters long.");

	// New
	if (!strcmp("new", cmd)) {
		// Check parameters
		if (nlhs != 1)
			mexErrMsgTxt("New: One output expected.");
		// Return a handle to a new C++ instance
		plhs[0] = convertPtr2Mat<SIMcore>(new SIMcore);
		return;
	}

	// Check there is a second input, which should be the class instance handle
	if (nrhs < 2)
		mexErrMsgTxt("Second input should be a class instance handle.");

	// Delete
	if (!strcmp("delete", cmd)) {
		// Destroy the C++ object
		destroyObject<SIMcore>(prhs[1]);
		// Warn if other commands were ignored
		if (nlhs != 0 || nrhs != 2)
			mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
		return;
	}
	
	
	// Get the class instance pointer from the second input
	SIMcore *SIMcore_instance = convertMat2Ptr<SIMcore>(prhs[1]);

	if (!strcmp("initialize", cmd)) {
		SIMcore_instance->initialize(prhs[2], prhs[3]);
		return;
	}
    
    
    if (!strcmp("initializeFromState", cmd)) {
		SIMcore_instance->initializeFromState(prhs[2], prhs[3]);
		return;
	}
    
    
    if (!strcmp("readState", cmd)) {
		SIMcore_instance->readState(prhs[2]);
		return;
	}

	if (!strcmp("getState", cmd)) {
		SIMcore_instance->getState(&plhs);
		return;
	}

	if (!strcmp("TU_go_grow_die", cmd)) {
		SIMcore_instance->TU_go_grow_die();
		return;
	}

	if (!strcmp("modulateAdjuvanticity", cmd)) {
		SIMcore_instance->modulateAdjuvanticity();
		return;
	}

	if (!strcmp("decayAdjuvanicityMap", cmd)) {
		SIMcore_instance->decayAdjuvanicityMap();
		return;
	}

	
	if (!strcmp("IMinflux", cmd)) {
        if (nrhs > 2) {
            SIMcore_instance->IMinfluxInput(prhs[2]);
        } else {
            SIMcore_instance->IMinflux();
        }
		return;
	}
	
	if (!strcmp("MPinflux", cmd)) {
        if (nrhs > 2) {
            SIMcore_instance->MPinfluxInput(prhs[2]);
        } else {
            SIMcore_instance->MPinflux();
        }
		return;
	}

	if (!strcmp("lymphocytesAct", cmd)) {
		SIMcore_instance->lymphocytesAct();
		return;
	}

	if (!strcmp("macrophagesAct", cmd)) {
		SIMcore_instance->macrophagesAct();
		return;
	}
	

	if (!strcmp("updateNecroMap", cmd)) {
		SIMcore_instance->updateNecroMap();
		return;
	}
    
    	if (!strcmp("updateChemoMap", cmd)) {
		SIMcore_instance->updateChemoMap();
		return;
	}
    
	if (!strcmp("seedFibrosis", cmd)) {
		SIMcore_instance->seedFibrosis();
		return;
	}
	
	if (!strcmp("TUcellsNum", cmd)) {
		SIMcore_instance->TUcellsNum(&plhs);
		return;
	}
            
	// Got here, so command not recognized
	mexErrMsgTxt("Command not recognized.");
}
