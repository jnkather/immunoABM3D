/*
*  Created by Jan Poleszczuk
*  Last modified November, 2017 by Jan Poleszczuk
*/

#include "stdafx.h"
#include "Macrophages.h"

size_t Macrophages::NumMacrophages() {
	return cells.size();
}

void Macrophages::initialize(SIMparameters *parSet, Environment *envSet, CRandomMersenne* rG) {
	env = envSet;
	randGen = rG;
	params = parSet;

	cells.clear(); //just in case if simulation is reinitialized
}

void Macrophages::initializeFromState(const mxArray* sysTempl, SIMparameters *parSet, Environment *envSet, CRandomMersenne* rG) {
	env = envSet;
	randGen = rG;
	params = parSet;

	cells.clear(); //just in case if simulation is reinitialized

	//reading from state
	mxArray *IM = mxGetField(sysTempl, 0, "MP");

	mxArray *val = mxGetField(IM, 0, "MPcells");

	const mwSize *dimensionPtr = mxGetDimensions(val);
	int numIMCells = (int)dimensionPtr[1]; //assuming it is column vector

	unsigned int* pos = (unsigned int*)mxGetPr(val);

	mxArray *IMprops = mxGetField(IM, 0, "MPprop");
	val = mxGetField(IMprops, 0, "Pcap");
	unsigned char* Pc = (unsigned char*)mxGetPr(val);
	val = mxGetField(IMprops, 0, "State");
	unsigned char* state = (unsigned char*)mxGetPr(val);

	for (int i = 0; i < numIMCells; ++i) {
		Macrophage cellIMAdd = { pos[i]-1, Pc[i], state[i] };
		env->setOccupancy(pos[i]-1, true);
		cells.push_back(cellIMAdd);
	}
}

void Macrophages::action() {

	//MOVEMENT PART
	unsigned int newSite;
	for (int iter = 0; iter < params->MPspeed; ++iter) {//number of time action to be taken

		//unsigned int newSite;
		std::random_shuffle(cells.begin(), cells.end(), std::bind(&CRandomMersenne::Rand, *randGen, std::placeholders::_1)); //shuffling cells

		for (size_t i = 0; i < cells.size(); ++i){//go through each macrophage
			if (randGen->Random() < params->MPpmig) { //if to migrate
				newSite = env->returnEmptyPlaceImmune(cells.at(i).place, params->MPrwalk);
				if (newSite) {
					env->setOccupancy(cells.at(i).place, false);
					cells.at(i).place = newSite;
					env->setOccupancy(newSite, true);
				}
			}
		}
	}//MP speed end

	//GO GROW DIE DIFF PART
	std::vector<Macrophage> cellsTmp;
	cellsTmp.reserve(cells.size());
	Macrophage currCell, newCell;

		std::random_shuffle(cells.begin(), cells.end(), std::bind(&CRandomMersenne::Rand, *randGen, std::placeholders::_1)); //shuffling cells

		while (!cells.empty()) {//go through each lymhpcyte
			currCell = cells.back(); //pick the cell
			cells.pop_back();
			

				if (randGen->Random() < params->MPpdeath) {
					env->setOccupancy(currCell.place, false);
				}
				else {
					newSite = env->returnEmptyPlaceImmune(currCell.place, params->MPrwalk);
					if (newSite) {
						if (randGen->Random() < params->MPpprol) { // if to proliferate
                            // cell wants to divide
                            if (params->CellDieAtProl && (randGen->Random() < params->CellDieAtProl)) {
                                env->setOccupancy(currCell.place, false);
                            }
                            else
                            {
							if (currCell.p > 0) {//there is proliferation capacity
								currCell.p--;
								newCell = currCell;
								newCell.place = newSite;
								env->setOccupancy(newSite, true);
								cellsTmp.push_back(currCell);
								cellsTmp.push_back(newCell);
							}
							else {//proliferation capacity exhausted
								env->setOccupancy(currCell.place, false);
							}
                            }
						}
						else if (randGen->Random() < params->MPpmig) { //if to migrate
							env->setOccupancy(currCell.place, false);
							currCell.place = newSite;
							env->setOccupancy(newSite, true);
							cellsTmp.push_back(currCell);
						}
						else {//do nothing
							cellsTmp.push_back(currCell);
						}
					}
					else {//no free site
						cellsTmp.push_back(currCell);
					}
				}//end if not dying

		}//for each cell
		cells.swap(cellsTmp);


		//POLARIZATION REPOLARIZATION
		std::vector<unsigned int> stateTwo;
		stateTwo.reserve(cells.size());
		for (size_t i = 0; i < cells.size(); ++i){//go through each macrophage
			if (cells.at(i).state == 1 && randGen->Random() < params->MPppola) {
				cells.at(i).state = 2;
				stateTwo.push_back(cells.at(i).place);
			}
			else if (cells.at(i).state == 2){
				if (randGen->Random() < params->MPprepola) {
					cells.at(i).state = 1;
				}
				else {
					stateTwo.push_back(cells.at(i).place);
				}
			}
		}

	//MODULATE ADJUVANTICITY
	env->modulateAdjuvanticity(stateTwo.begin(), stateTwo.end(), params->MPeffect);
	stateTwo.clear();

}

void Macrophages::influx() {

	int numNewCells = params->MPinfluxRate;
	
	unsigned int* positions = env->generatePostions(&numNewCells);

	if (positions != NULL) {//there are some free spots

		for (int i = 0; i < numNewCells; ++i) {
			env->setOccupancy(positions[i], true);
			Macrophage macrophage = { positions[i], params->MPpmax, 1};
			cells.push_back(macrophage);
		}

		delete[] positions;
	}
}


void Macrophages::influxInput(const mxArray* numNewCellIn) {

	int numNewCells = (int)mxGetPr(numNewCellIn)[0];
	//printf("Adding %d new Macrophages \n", numNewCells);

	unsigned int* positions = env->generatePostions(&numNewCells);

	if (positions != NULL) {//there are some free spots

		for (int i = 0; i < numNewCells; ++i) {
			env->setOccupancy(positions[i], true);
			Macrophage macrophage = { positions[i], params->MPpmax, 1 };
			cells.push_back(macrophage);
		}

		delete[] positions;
	}
}

mxArray* Macrophages::getState() {

	const char *field_names[] = { "position", "Pcap", "state" };
	mwSize dimsStruct[2] = { 1, 1 };
	mxArray* IMcells = mxCreateStructArray(1, dimsStruct, 3, field_names);

	if (cells.size() > 0) {//if there is any lymphocyte
						   //2 - copying Lymphocytes positions
		mwSize dims[2];
		dims[0] = 1; dims[1] = cells.size();
		mxSetFieldByNumber(IMcells, 0, mxGetFieldNumber(IMcells, "position"),
			mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL));
		int *IMcellsPos = (int*)mxGetPr(mxGetField(IMcells, 0, "position"));
		for (size_t i = 0; i < cells.size(); ++i)
			IMcellsPos[i] = cells.at(i).place + 1;

		//3 - copying properties
		mxSetFieldByNumber(IMcells, 0, mxGetFieldNumber(IMcells, "Pcap"),
			mxCreateNumericArray(2, dims, mxUINT8_CLASS, mxREAL));
		unsigned char *Pcap = (unsigned char*)mxGetPr(mxGetField(IMcells, 0, "Pcap"));
		mxSetFieldByNumber(IMcells, 0, mxGetFieldNumber(IMcells, "state"),
			mxCreateNumericArray(2, dims, mxUINT8_CLASS, mxREAL));
		unsigned char *st = (unsigned char*)mxGetPr(mxGetField(IMcells, 0, "state"));

		for (size_t i = 0; i < cells.size(); ++i) {
			Pcap[i] = cells.at(i).p;
			st[i] = cells.at(i).state;
		}
	}

	return IMcells;
}
