/*
*  Created by Jan Poleszczuk
*  Last modified November, 2017 by Jan Poleszczuk
*/

#include "stdafx.h"
#include "Lymphocytes.h"

size_t Lymphocytes::NumLymphocytes() {
	return cells.size();
}

void Lymphocytes::initialize(SIMparameters *parSet, Environment *envSet, CRandomMersenne* rG, TumorCells* TUc) {
	env = envSet;
	randGen = rG;
	params = parSet;
	TUcells = TUc;
	
	cells.clear(); //just in case if simulation is reinitialized
}


void Lymphocytes::initializeFromState(const mxArray* sysTempl, SIMparameters *parSet, Environment *envSet, CRandomMersenne* rG, TumorCells* TUc) {
	env = envSet;
	randGen = rG;
	params = parSet;
	TUcells = TUc;

	cells.clear(); //just in case if simulation is reinitialized

	//reading from state
	mxArray *IM = mxGetField(sysTempl, 0, "IM");

	mxArray *val = mxGetField(IM, 0, "IMcells");

	const mwSize *dimensionPtr = mxGetDimensions(val);
	int numIMCells = (int)dimensionPtr[1]; //assuming it is column vector

	unsigned int* pos = (unsigned int*)mxGetPr(val);

	mxArray *IMprops = mxGetField(IM, 0, "IMprop");
	val = mxGetField(IMprops, 0, "Kcap");
	unsigned char* Kcap = (unsigned char*)mxGetPr(val);
	val = mxGetField(IMprops, 0, "Pcap");
	unsigned char* Pc = (unsigned char*)mxGetPr(val);
	val = mxGetField(IMprops, 0, "engaged");
	unsigned char* eng = (unsigned char*)mxGetPr(val);

	for (int i = 0; i < numIMCells; ++i) {
		Lymphocyte cellIMAdd = { pos[i]-1, Pc[i], Kcap[i], eng[i] };
		env->setOccupancy(pos[i]-1, true);
		cells.push_back(cellIMAdd);
	}
}

void Lymphocytes::fibrosify() {

	Lymphocyte currCell;
	std::vector<Lymphocyte> cellsTmp;
	cellsTmp.reserve(cells.size());

	std::vector<unsigned int> seeds;
	seeds.reserve(cells.size());

	while (!cells.empty()) {//go through each lymhpcyte
		currCell = cells.back(); //pick the cell
		cells.pop_back();
		if (currCell.kcap == 0 && randGen->Random() < params->probSeedFibr) {
			env->setOccupancy(currCell.place, false);
			seeds.push_back(currCell.place);
		}
		else {
			cellsTmp.push_back(currCell);
		}
	}

	cells.swap(cellsTmp);//copy remaining cells

	if (!seeds.empty())
		env->modulateFibrosis(seeds.begin(), seeds.end());
}

void Lymphocytes::action() {
	
	//create new lattice with tumor cell position
	//std::map<unsigned int, std::pair<unsigned int, float>> indices;
	std::map<unsigned int, std::vector<TumorCells::TumorCell>::iterator> TUcellsMap;
	std::map<unsigned int, std::vector<TumorCells::TumorCell>::iterator>::iterator it;

	TUcells->getStateForImmune(&TUcellsMap);

	std::vector<Lymphocyte> cellsTmp;
	cellsTmp.reserve(cells.size());


	for (int iter = 0; iter < params->IMspeed; ++iter) {//number of time action to be taken
		
		cellsTmp.clear();
		Lymphocyte currCell, newCell;
		unsigned int newSite;

		//unsigned int newSite;
		std::random_shuffle(cells.begin(), cells.end(), std::bind(&CRandomMersenne::Rand, *randGen, std::placeholders::_1)); //shuffling cells

		while (!cells.empty()) {//go through each lymhpcyte
			currCell = cells.back(); //pick the cell
			cells.pop_back();
			bool skipAttack = false;

			if (currCell.engaged > 0) {//check if engaged
				currCell.engaged--;
				cellsTmp.push_back(currCell); //just push back and do nothing
			}
			else {

				if (randGen->Random() < params->IMpdeath) {
					env->setOccupancy(currCell.place, false);
					skipAttack = true;
				}
				else {
					newSite = env->returnEmptyPlaceImmune(currCell.place, params->IMrwalk);
					if (newSite) {
						if (randGen->Random() < params->IMpprol) { // if to proliferate
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
								newCell.kcap = params->IMkmax;
								env->setOccupancy(newSite, true);
								cellsTmp.push_back(currCell);
								cellsTmp.push_back(newCell);
							}
							else {//proliferation capacity exhausted
								env->setOccupancy(currCell.place, false);
							}
                            }
							skipAttack = true;
						}
						else if (randGen->Random() < params->IMpmig) { //if to migrate
							env->setOccupancy(currCell.place, false);
							currCell.place = newSite;
							env->setOccupancy(newSite, true);
							//cellsTmp.push_back(currCell);
						}
						else {//do nothing
							//cellsTmp.push_back(currCell);
						}
					}
					else {//no free site
						//cellsTmp.push_back(currCell);
					}
				}//end if not dying

				//handling the possible attack event
				if (!skipAttack) {
					
					if (currCell.kcap > 0 && randGen->Random() < params->IMpkill) { //this cell is ready to attack
						//find the target
						unsigned int target = env->findTarget(currCell.place);
						if (target) {
							it = TUcellsMap.find(target);
							if (it != TUcellsMap.end()) {
								if (it->second->Antigen >= params->antiTresh && env->getAdjuValue(currCell.place) >= params->adjuTresh) {//successful kill
									//deposit information about attack in the tumor cell
									currCell.kcap--;
									currCell.engaged = params->engagementDuration;
									if (it->second->damage < params->damageTresh)
										it->second->damage++;
								}
							}
							else {
								printf("Warning: Target not found!!\n");
							}
						}
					}

					cellsTmp.push_back(currCell);
				}//end if not to skip attack
			}//end if not engaged

		}//for each cell
		cells.swap(cellsTmp);
	}//IM speed

}

void Lymphocytes::influx() {
	
	int numNewCells = params->IMinfluxRate;
	
	unsigned int* positions = env->generatePostions(&numNewCells);

	if (positions != NULL) {//there are some free spots

		for (int i = 0; i < numNewCells; ++i) {
			env->setOccupancy(positions[i], true);
			Lymphocyte lymphocyte = {positions[i], params->IMpmax, params->IMkmax, 0 };
			cells.push_back(lymphocyte);
		}

		delete[] positions;
	}
}

void Lymphocytes::influxInput(const mxArray* numNewCellIn) {

	int numNewCells = (int)mxGetPr(numNewCellIn)[0];
	//printf("Adding %d new T cells \n", numNewCells);

	unsigned int* positions = env->generatePostions(&numNewCells);

	if (positions != NULL) {//there are some free spots

		for (int i = 0; i < numNewCells; ++i) {
			env->setOccupancy(positions[i], true);
			Lymphocyte lymphocyte = { positions[i], params->IMpmax, params->IMkmax, 0 };
			cells.push_back(lymphocyte);
		}

		delete[] positions;
	}
}

mxArray* Lymphocytes::getState() {

	const char *field_names[] = { "position", "Pcap", "Kcap","engaged" };
	mwSize dimsStruct[2] = { 1, 1 };
	mxArray* IMcells = mxCreateStructArray(1, dimsStruct, 4, field_names);

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
		mxSetFieldByNumber(IMcells, 0, mxGetFieldNumber(IMcells, "Kcap"),
			mxCreateNumericArray(2, dims, mxUINT8_CLASS, mxREAL));
		unsigned char *Kcap = (unsigned char*)mxGetPr(mxGetField(IMcells, 0, "Kcap"));
		mxSetFieldByNumber(IMcells, 0, mxGetFieldNumber(IMcells, "engaged"),
			mxCreateNumericArray(2, dims, mxUINT8_CLASS, mxREAL));
		unsigned char *eng = (unsigned char*)mxGetPr(mxGetField(IMcells, 0, "engaged"));


		for (size_t i = 0; i < cells.size(); ++i) {
			Pcap[i] = cells.at(i).p;
			Kcap[i] = cells.at(i).kcap;
			eng[i] = cells.at(i).engaged;
		}
	}

	return IMcells;
}
