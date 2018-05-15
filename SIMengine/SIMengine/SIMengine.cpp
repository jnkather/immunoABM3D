/*
*  Created by Jan Poleszczuk
*  Last modified September, 2017 by Jan Poleszczuk
*/


// SIMengine.cpp : Defines the entry point for the console application.
//

#include "Classes/stdafx.h"

#include "Classes/SIMcore.h"

void runSimulation() {
	//initialize simulation
	SIMcore simCore;
	simCore.initializeMexFree();


	for (int i = 1; i < 140; ++i) {
		//simCore.TU_go_grow_die();
		if (simCore.NumTUcells() > 0) {
			std::cout << simCore.NumTUcells() << std::endl;
		}
		else {
			std::cout << "All TU cells died." << std::endl;
			break;
		}

		simCore.modulateAdjuvanticity();
		
		simCore.updateNecroMap();
		simCore.updateChemoMap();
		simCore.decayAdjuvanicityMap();

		if (simCore.generateRandom() < 0)
			simCore.IMinflux();

		if (simCore.generateRandom() < 0)
			simCore.MPinflux();

		simCore.lymphocytesAct();

		simCore.macrophagesAct();

		simCore.seedFibrosis();
	}
	
	system("pause");
}

int main()
{
	runSimulation();
	
	_CrtDumpMemoryLeaks();//this is for detecting memory leaks if in DEBUG mode
    return 0;
}

