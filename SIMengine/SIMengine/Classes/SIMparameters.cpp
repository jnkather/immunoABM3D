/*
*  Created by Jan Poleszczuk
*  Last modified September, 2017 by Jan Poleszczuk
*/


#include "SIMparameters.h"


SIMparameters::SIMparameters()
{
	// START GENERAL SYSTEM PROPERTIES -------------------------------------------
	initSeed = 0;
	seedUnderneath = true;
    CellDieAtProl = 0;
	// END SYSTEM PROPERTIES -------------------------------------------

	// START INITIALIZE TUMOR CELLS -------------------------------------------
	TUpmax = 10; //divisions before proliferative capacity exhaustion
	TUpprol = 1. / 24.; //division probability
	TUpdeath = 0.05; //spontaneous death probability
	TUps = 0.3; //probability of symmetric division
	TUpmig = 10. / 24.; //probability of migration
	TUpmut = 0.;
	TUdanti = 0.1f;
	TUdadju = 1.f;
	damageTresh = 5;
	// END INITIALIZE TUMOR CELLS ---------------------------------------------


	// START INITIALIZE LYMPHOCYTES ------------------------------------------
	IMkmax = 5;
	IMpmax = 10;
	IMpmig = 0.8;
	IMpkill = 0.1;
	IMrwalk = 0.8f;
	IMspeed = 30;
	IMpprol = 0.049 / (double)IMspeed;
	IMpdeath = 0.0147;
	engagementDuration = 48;
	adjuTresh = 0.0f;
	antiTresh = 0.3f;
	adjuDecay = 0.95f;
	IMinfluxRate = 1;
	// END INITIALIZE LYMPHOCYTES --------------------------------------------



	// START INITIALIZE MACROPHAGES ------------------------------------------
	MPpmax = 10;
	MPpprol = 0.049 / (double)MPspeed;
	MPpmig = 0.8;
	MPpdeath = 0.0147;
	MPrwalk = 0.8f;
	MPspeed = 30;
	MPppola = 0.1;
	MPprepola = 0.1;
	MPeffect = -0.1f;
	adjuRange = 7;
	MPinfluxRate = 1;
	// END INITIALIZE MACROPHAGES --------------------------------------------

	// START INITIALIZE CHEMOTAXIS MAP------------------------------------------
	DCchemo = 1440.;
	SCchemo = 10.;
	// END INITIALIZE CHEMOTAXIS MAP------------------------------------------
	

	// START INITIALIZE NECROSIS------------------------------------------
	DCnecro = 1.;
	TCnecro = 1.;
	necThresh = 0.5;
	// END INITIALIZE NECROSIS------------------------------------------
	
	// START INITIALIZE FIBROSIS  ---------------------------------
	smoothRadius = 3;
	probSeedFibr = 0.008;
	fibrFrac = 0.3;
	stromaPerm = 0.1;
	// END INITIALIZE FIBROSIS  ---------------------------------
	
	defaultAntigenicity = 0.05f;
	maxAntigenicity = 1.;
	

	N1 = 351; N2 = 351; N3 = 1;
}


SIMparameters::~SIMparameters(){};


void SIMparameters::listSIMParameters() {
	printf("Current simulation parameters in SIMengine: \n");
	printf("_____________________________________________ \n");
	printf("TUpmax: %d \n", TUpmax);
	printf("TUpprol: %f \n", TUpprol);
	printf("TUpdeath: %f \n", TUpdeath);
	printf("TUps: %f \n", TUps);
	printf("TUpmig: %f \n", TUpmig);
	printf("TUpmut: %f \n", TUpmut);
	printf("TUdanti: %f \n", TUdanti);
	printf("TUdadju: %f \n", TUdadju);
	printf("adjuRange: %d \n", adjuRange);
	printf("adjuDecay: %f \n", adjuDecay);
	printf("DCnecro: %f \n", DCnecro);
	printf("necThresh: %f \n", necThresh);
	printf("DCchemo: %f \n", DCchemo);
	printf("SCchemo: %f \n", SCchemo);
	printf("defaultAntigenicity: %f \n", defaultAntigenicity);
	printf("maxAntigenicity: %f \n", maxAntigenicity);
	printf("IMpmax: %d \n", IMpmax);
	printf("IMkmax: %d \n", IMkmax);
	printf("IMpmig: %f \n", IMpmig);
	printf("IMpprol: %f \n", IMpprol);
	printf("IMpdeath: %f \n", IMpdeath);
	printf("IMpkill: %f \n", IMpkill);
	printf("IMrwalk: %f \n", IMrwalk);
	printf("IMinfluxRate: %d \n", IMinfluxRate);
	printf("IMspeed: %d \n", IMspeed);
	printf("probSeedFibr: %f \n", probSeedFibr);
	printf("fibrFrac: %f \n", fibrFrac);
	printf("damageTresh: %d \n", damageTresh);
	printf("engagementDuration: %d \n", engagementDuration);
	printf("antiTresh: %f \n", antiTresh);
	printf("adjuTresh: %f \n", adjuTresh);
	printf("MPpmax: %d \n", MPpmax);
	printf("MPinfluxRate: %d \n", MPinfluxRate);
	printf("stromaPerm: %f \n", stromaPerm);
	printf("initSeed: %f \n", initSeed);
	printf("smoothRadius: %d \n", smoothRadius);
	printf("N1xN2xN3: %dx%dx%d \n", N1, N2, N3);
}

void SIMparameters::initialize(const mxArray* sysTempl, const mxArray* cnst) {
	//parsing input parameters
	mxArray *params = mxGetField(sysTempl, 0, "params");
    CellDieAtProl = mxGetPr(mxGetField(params, 0, "CellDieAtProl"))[0];
	TUpmax = (unsigned char)mxGetPr(mxGetField(params, 0, "TUpmax"))[0]; //divisions before proliferative capacity exhaustion
	TUpprol = mxGetPr(mxGetField(params, 0, "TUpprol"))[0]; //division probability
	TUpdeath = mxGetPr(mxGetField(params, 0, "TUpdeath"))[0]; //spontaneous death probability
	TUps = mxGetPr(mxGetField(params, 0, "TUps"))[0]; //probability of symmetric division
	TUpmig = mxGetPr(mxGetField(params, 0, "TUpmig"))[0]; //probability of migration
	TUpmut = mxGetPr(mxGetField(params, 0, "TUpmut"))[0]; //probability of mutation
	TUdanti = (float)mxGetPr(mxGetField(params, 0, "TUdanti"))[0];
	TUdadju = (float)mxGetPr(mxGetField(params, 0, "TUdadju"))[0];

	initSeed = mxGetPr(mxGetField(params, 0, "initialSeed"))[0];
	stromaPerm = (double)mxGetPr(mxGetField(params, 0, "stromaPerm"))[0];

	adjuRange = (int)mxGetPr(mxGetField(params, 0, "adjuRange"))[0];
	adjuDecay = (float)mxGetPr(mxGetField(params, 0, "adjuDecay"))[0];

	smoothRadius = (int)mxGetPr(mxGetField(params, 0, "smoothRadius"))[0];
    seedUnderneath = mxGetPr(mxGetField(params, 0, "seedUnderneath"))[0] > 0;

	DCnecro = mxGetPr(mxGetField(params, 0, "DCnecro"))[0];
	TCnecro = mxGetPr(mxGetField(params, 0, "TCnecro"))[0];
	necThresh = mxGetPr(mxGetField(params, 0, "necThresh"))[0];

	DCchemo = mxGetPr(mxGetField(params, 0, "DCchemo"))[0];
	SCchemo = mxGetPr(mxGetField(params, 0, "SCchemo"))[0];

	IMpmax = (unsigned char)mxGetPr(mxGetField(params, 0, "IMpmax"))[0];
	IMkmax = (unsigned char)mxGetPr(mxGetField(params, 0, "IMkmax"))[0];
	IMpmig = mxGetPr(mxGetField(params, 0, "IMpmig"))[0];
	IMpprol = mxGetPr(mxGetField(params, 0, "IMpprol"))[0];
	IMpdeath = mxGetPr(mxGetField(params, 0, "IMpdeath"))[0];
	IMpkill = mxGetPr(mxGetField(params, 0, "IMpkill"))[0];
	IMrwalk = (float)mxGetPr(mxGetField(params, 0, "IMrwalk"))[0];
	probSeedFibr = mxGetPr(mxGetField(params, 0, "probSeedFibr"))[0];
	fibrFrac = mxGetPr(mxGetField(params, 0, "fibrFrac"))[0];

	IMinfluxRate = (int)mxGetPr(mxGetField(params, 0, "IMinfluxRate"))[0];
	IMspeed = (int)mxGetPr(mxGetField(params, 0, "IMspeed"))[0];

	damageTresh = (unsigned char)mxGetPr(mxGetField(params, 0, "TUdamageThresh"))[0];
	engagementDuration = (unsigned char)mxGetPr(mxGetField(params, 0, "engagementDuration"))[0];

	antiTresh = (float)mxGetPr(mxGetField(params, 0, "antiThresh"))[0];
	adjuTresh = (float)mxGetPr(mxGetField(params, 0, "adjuThresh"))[0];

	MPpmax = (unsigned char)mxGetPr(mxGetField(params, 0, "MPpmax"))[0];
	MPinfluxRate = (int)mxGetPr(mxGetField(params, 0, "MPinfluxRate"))[0];

	MPspeed = (int)mxGetPr(mxGetField(params, 0, "MPspeed"))[0];
	MPpmig = mxGetPr(mxGetField(params, 0, "MPpmig"))[0];
	MPrwalk = (float)mxGetPr(mxGetField(params, 0, "MPrwalk"))[0];
	MPpdeath = mxGetPr(mxGetField(params, 0, "MPpdeath"))[0];
	MPpprol = mxGetPr(mxGetField(params, 0, "MPpprol"))[0];
	MPppola = mxGetPr(mxGetField(params, 0, "MPppola"))[0];
	MPprepola = mxGetPr(mxGetField(params, 0, "MPprepola"))[0];
	MPeffect = (float)mxGetPr(mxGetField(params, 0, "MPeffect"))[0];


	mxArray *grid = mxGetField(sysTempl, 0, "grid");
	N1 = (unsigned int)mxGetPr(mxGetField(grid, 0, "N"))[0];
	N2 = (unsigned int)mxGetPr(mxGetField(grid, 0, "M"))[0];
	if (mxGetField(grid, 0, "P") != NULL) {
		N3 = (unsigned int)mxGetPr(mxGetField(grid, 0, "P"))[0];
	}
	else {
		N3 = 1;
	}

	defaultAntigenicity = (float)mxGetPr(mxGetField(cnst, 0, "defaultAntigenicity"))[0];
	maxAntigenicity = (float)mxGetPr(mxGetField(cnst, 0, "maxAntigenicity"))[0];
}
