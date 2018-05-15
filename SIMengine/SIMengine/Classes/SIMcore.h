/*
 *  Created by Jan Poleszczuk
 *  Last modified September, 2017 by Jan Poleszczuk
 */


#pragma once

#include "stdafx.h"
#include "SIMparameters.h"
#include "Environment.h"
#include "TumorCells.h"
#include "Lymphocytes.h"
#include "Macrophages.h"

class SIMcore
{
private:
	SIMparameters params;
	Environment env;
	TumorCells TUcells;
	Lymphocytes lymphocytes;
	Macrophages macrophages;

	CRandomMersenne randGen;
	
public:
	SIMcore();
	~SIMcore();

	void TU_go_grow_die(){ TUcells.go_grow_die(); };
	void modulateAdjuvanticity() { TUcells.modulateAdjuvanticity(); };
	void decayAdjuvanicityMap() { env.decayAdjuvanicityMap(); };
	void IMinflux() { lymphocytes.influx(); };
	void MPinflux() { macrophages.influx(); };
	void IMinfluxInput(const mxArray* numNewCellsIn) { lymphocytes.influxInput(numNewCellsIn); };
	void MPinfluxInput(const mxArray* numNewCellsIn) { macrophages.influxInput(numNewCellsIn); };

	void lymphocytesAct() { lymphocytes.action(); };
	void macrophagesAct() { macrophages.action(); };
	void updateNecroMap() { env.updateNecroMap(); };
	void updateChemoMap() { env.updateChemoMap(); };
	void seedFibrosis() { lymphocytes.fibrosify(); };

	void initialize(const mxArray*, const mxArray*);
	void initializeFromState(const mxArray*, const mxArray*);
	void readState(const mxArray*);
	void initializeMexFree();
	void getState(mxArray***);
	void TUcellsNum(mxArray***);
	
	double generateRandom(){ return randGen.Random(); };
	
	size_t NumTUcells();
};

