/*
 *  Created by Jan Poleszczuk
 *  Last modified September, 2017 by Jan Poleszczuk
 */


#pragma once
#include "stdafx.h"

class SIMparameters
{
public:

	//tumor cell related parameters
	unsigned char TUpmax; //divisions before proliferative capacity exhaustion
    
    double CellDieAtProl; // cell dies at proliferation attempt, default 0, increase by chemo
	double TUpprol; //division probability
	double TUpdeath; //spontaneous death probability
	double TUps; //probability of symmetric division
	double TUpmig; //probability of migration
	double TUpmut;
	float TUdanti;

	float defaultAntigenicity;
	double stromaPerm;
	float maxAntigenicity;

	float TUdadju;
	int adjuRange;
	float adjuDecay;

	int smoothRadius;
    bool seedUnderneath;
	
	//chemotaxis map associated parameters
	double DCchemo;
	double SCchemo;

	//necrosis map associated parameters
	double DCnecro;
	double TCnecro;
	double necThresh;


	//lymphocytes assiociated parameters
	unsigned char IMpmax;
	unsigned char IMkmax;
	double IMpmig;
	double IMpprol;
	double IMpdeath;
	float IMrwalk;
	double IMpkill;
	int IMinfluxRate;
	int IMspeed;
	double probSeedFibr;
	double fibrFrac;

	float antiTresh;
	float adjuTresh;

	unsigned char damageTresh;
	unsigned char engagementDuration;

	//macrophages associated parameters
	unsigned char MPpmax;
	int MPinfluxRate;
	int MPspeed;
	double MPpmig;
	float MPrwalk;
	double MPpdeath;
	double MPpprol;
	double MPppola;
	double MPprepola;
	float MPeffect;

	unsigned int N1, N2, N3;; //environment dimensions

	double initSeed;

	SIMparameters();
	~SIMparameters();

	void listSIMParameters();

	void initialize(const mxArray*, const mxArray*);

};

