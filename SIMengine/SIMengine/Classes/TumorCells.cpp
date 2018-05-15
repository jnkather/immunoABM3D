/*
*  Created by Jan Poleszczuk
*  Last modified November, 2017 by Jan Poleszczuk
*/

#include "stdafx.h"
#include "TumorCells.h"

size_t TumorCells::NumTUcells() {
	return cells.size();
}

void TumorCells::initialize(SIMparameters *parSet, Environment *envSet, CRandomMersenne* rG) {
	env = envSet;
	randGen = rG;
	params = parSet;

	//intializing with single cancer stem cell in the lattice center
	unsigned int center = env->getCenter();
	env->newTUcell(center);
	TumorCell initialCell = { center, params->TUpmax, true, params->defaultAntigenicity, 0 };
	cells.clear(); //just in case if simulation is reinitialized
	cells.push_back(initialCell);

	dyingCells.clear();
}

void TumorCells::initializeFromState(const mxArray* sysTempl, SIMparameters *parSet, Environment *envSet, CRandomMersenne* rG) {
	env = envSet;
	randGen = rG;
	params = parSet;

	cells.clear(); //just in case
	dyingCells.clear();

	//reading from state
	mxArray *TU = mxGetField(sysTempl, 0, "TU");

	mxArray *val = mxGetField(TU, 0, "TUcells");

	const mwSize *dimensionPtr = mxGetDimensions(val);
	int numTUCells = (int)dimensionPtr[1]; //assuming it is column vector

	unsigned int* pos = (unsigned int*)mxGetPr(val);
	
	mxArray *TUprops = mxGetField(TU, 0, "TUprop");
	val = mxGetField(TUprops, 0, "isStem");
	bool* iS = (bool*)mxGetPr(val);
	val = mxGetField(TUprops, 0, "Pcap");
	unsigned char* Pc = (unsigned char*)mxGetPr(val);
	val = mxGetField(TUprops, 0, "Antigen");
	float* An = (float*)mxGetPr(val);
	val = mxGetField(TUprops, 0, "damage");
	unsigned char* dem = (unsigned char*)mxGetPr(val);

	for (int i = 0; i < numTUCells; ++i) {
		TumorCell cellTuAdd = { pos[i]-1, Pc[i], iS[i], An[i], dem[i] };
		env->newTUcell(pos[i]-1);
		cells.push_back(cellTuAdd);
	}
}

void TumorCells::readState(mxArray* init) {
	

}

void TumorCells::getStateForImmune(std::map<unsigned int, std::vector<TumorCells::TumorCell>::iterator>* TUmap) {

	std::vector<TumorCell>::iterator it;
	for (it = cells.begin(); it != cells.end(); it++) {
		(*TUmap)[it->place] = it;
	}

};

void TumorCells::modulateAdjuvanticity() {
	if (!dyingCells.empty()) {
		env->modulateAdjuvanticity(dyingCells.begin(), dyingCells.end(),params->TUdadju);
	}
}

void TumorCells::go_grow_die() {

	dyingCells.clear();

	unsigned int newSite;
	TumorCell currCell, newCell;
	std::vector<TumorCell> cellsTmp;
	cellsTmp.reserve(cells.size());

	std::random_shuffle(cells.begin(), cells.end(), std::bind(&CRandomMersenne::Rand, *randGen, std::placeholders::_1)); //shuffling cells
	while (!cells.empty()) {
		currCell = cells.back(); //pick the cell
		cells.pop_back();

		if (currCell.damage == params->damageTresh || (!currCell.is_stem && randGen->Random() < params->TUpdeath)) {//cell is killed by the immune system, dies spontaneously
			env->deleteTUcell(currCell.place);
			dyingCells.push_back(currCell.place);
		}
		else if (env->getNecrosisValue(currCell.place) < params->necThresh && randGen->Random() > env->getNecrosisValue(currCell.place)){//necrtoic death
			env->deleteTUcell(currCell.place);
			env->markNecrosis(currCell.place);
			dyingCells.push_back(currCell.place);
		}
		else {

			newSite = env->returnEmptyPlace(currCell.place);

			if (newSite) {//if there is a new spot
				newCell = currCell;
				newCell.place = newSite;
				if (randGen->Random() < params->TUpprol) {
                    // START Chemotherapy
                    // tumor cell wants to divide
                    if (params->CellDieAtProl && (randGen->Random() < params->CellDieAtProl)) {
                        // tumor cell will die
                        env->deleteTUcell(currCell.place);
						dyingCells.push_back(currCell.place);
                    }
                    // END Chemotherapy
                    else
                    {
					if (currCell.is_stem) {
						env->newTUcell(newSite);
						if (randGen->Random() > params->TUps) {//asymmetric division
							newCell.is_stem = false;
						}
						if (currCell.Antigen < params->maxAntigenicity && randGen->Random() < params->TUpmut) {
							currCell.Antigen += params->TUdanti;
							if (currCell.Antigen > params->maxAntigenicity)
								currCell.Antigen = params->maxAntigenicity;
							newCell.Antigen = currCell.Antigen;
						}
						newCell.damage = 0;
						cellsTmp.push_back(currCell);
						cellsTmp.push_back(newCell);
					}
					else if (currCell.p > 0) {
						currCell.p--;
						newCell.p--;
						env->newTUcell(newSite);
						if (currCell.Antigen < params->maxAntigenicity && randGen->Random() < params->TUpmut) {
							currCell.Antigen += params->TUdanti;
							if (currCell.Antigen > params->maxAntigenicity)
								currCell.Antigen = params->maxAntigenicity;
							newCell.Antigen = currCell.Antigen;
						}
						newCell.damage = 0;
						cellsTmp.push_back(currCell);
						cellsTmp.push_back(newCell);
					}
					else {
						env->deleteTUcell(currCell.place);
						dyingCells.push_back(currCell.place);
					}
                    }
				}
				else if (randGen->Random() < params->TUpmig) {
					env->deleteTUcell(currCell.place);
					env->newTUcell(newSite);
					cellsTmp.push_back(newCell);
				}
				else {//doing nothing
					cellsTmp.push_back(currCell);
				}
			}
			else {//no free spot
				cellsTmp.push_back(currCell);
			}
		}
	}
	cells.swap(cellsTmp);
}

mxArray* TumorCells::getState() {

	const char *field_names[] = { "position", "isStem", "Pcap", "Antigen","damage"};
	mwSize dimsStruct[2] = { 1, 1 };
	mxArray* TUcells = mxCreateStructArray(1, dimsStruct, 5, field_names);

	//2 - copying TUcells positions
	mwSize dims[2];
	dims[0] = 1; dims[1] = cells.size();
	mxSetFieldByNumber(TUcells, 0, mxGetFieldNumber(TUcells, "position"),
		mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL));
	int *TUcellsPos = (int*)mxGetPr(mxGetField(TUcells, 0, "position"));
	for (size_t i = 0; i < cells.size(); ++i)
		TUcellsPos[i] = cells.at(i).place + 1;

	//3 - copying TUprop
	mxSetFieldByNumber(TUcells, 0, mxGetFieldNumber(TUcells, "isStem"),
		mxCreateLogicalArray(2, dims));
	bool *isStem = (bool*)mxGetPr(mxGetField(TUcells, 0, "isStem"));

	mxSetFieldByNumber(TUcells, 0, mxGetFieldNumber(TUcells, "Pcap"),
		mxCreateNumericArray(2, dims, mxUINT8_CLASS, mxREAL));
	unsigned char *Pcap = (unsigned char*)mxGetPr(mxGetField(TUcells, 0, "Pcap"));

	mxSetFieldByNumber(TUcells, 0, mxGetFieldNumber(TUcells, "Antigen"),
		mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL));
	float *Antigen = (float*)mxGetPr(mxGetField(TUcells, 0, "Antigen"));

	mxSetFieldByNumber(TUcells, 0, mxGetFieldNumber(TUcells, "damage"),
		mxCreateNumericArray(2, dims, mxUINT8_CLASS, mxREAL));
	unsigned char *dam = (unsigned char*)mxGetPr(mxGetField(TUcells, 0, "damage"));

	for (size_t i = 0; i < cells.size(); ++i) {
		isStem[i] = cells.at(i).is_stem;
		Pcap[i] = cells.at(i).p;
		Antigen[i] = cells.at(i).Antigen;
		dam[i] = cells.at(i).damage;
	}

	return TUcells;
}
