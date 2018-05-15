/*
*  Created by Jan Poleszczuk
*  Last modified November, 2017 by Jan Poleszczuk
*/

#include "stdafx.h"
#include "Environment.h"

Environment::Environment() {
	lattice = NULL;
	indcNeigh = NULL;
	neigh = NULL;
	necrosis = NULL;
	fibrosis = NULL;
	AdjuMap = NULL;
	latticeOnlyTU = NULL;
	chemoGrad = NULL;
	//setting default lattice sizes
	N1 = 100;
	N2 = 100;
	N3 = 1;
	numSpots = 10000;
	is3D = false;
	numNeigh = 8;
	numOccupiedSpots = 0;
	numNecroticSpots = 0;
	numFibroticSpots = 0;
};

Environment::~Environment() {
	if (lattice != NULL)
		delete[] lattice;
	if (indcNeigh != NULL)
		delete[] indcNeigh;
	if (neigh != NULL)
		delete[] neigh;
	if (necrosis != NULL)
		delete[] necrosis;
	if (fibrosis != NULL)
		delete[] fibrosis;
	if (AdjuMap != NULL)
		delete[] AdjuMap;
	if (latticeOnlyTU != NULL)
		delete[] latticeOnlyTU;
	if (chemoGrad != NULL)
		delete[] chemoGrad;
}

unsigned int Environment::getCenter() {
	return N1 / 2 + (N2 / 2)*N1 + (N3 / 2)*N1*N2;
}

void Environment::setOccupancy(unsigned int pos, bool state) {
	lattice[pos] = state;
	if (state)
		numOccupiedSpots++;
	else
		numOccupiedSpots--;

}

void Environment::newTUcell(unsigned int pos) {
	lattice[pos] = true;
	latticeOnlyTU[pos] = true;
	necrosis[pos] = false;
	numOccupiedSpots++;
	ChtaxMap.addSource(pos);
	NecroMap.addSink(pos);
}

void Environment::deleteTUcell(unsigned int pos) {
	lattice[pos] = false;
	latticeOnlyTU[pos] = false;
	numOccupiedSpots--;
	ChtaxMap.removeSource(pos);
	NecroMap.removeSink(pos);
}

unsigned int* Environment::generatePostions(int* N) {

	unsigned int* positions = NULL;

	std::vector<unsigned int> availableSpots;
	availableSpots.reserve(numSpots - numOccupiedSpots - numNecroticSpots);

	for (unsigned int i = 0; i < numSpots; ++i)
		if (!lattice[i] && !necrosis[i] && (!fibrosis[i] || randGen->Random() < params->stromaPerm))
			availableSpots.push_back(i);

	if (availableSpots.size() < *N) {
		*N = (int)availableSpots.size();
		if (*N > 0) {
			positions = DBG_NEW unsigned int[*N];
			for (int i = 0; i < *N; ++i)
				positions[i] = availableSpots.at(i);
		}
	}
	else {
		positions = DBG_NEW unsigned int[*N];
		int gP;
		for (int i = 0; i < *N; ++i) {
			gP = randGen->Rand((int)availableSpots.size());
			positions[i] = availableSpots.at(gP);
			availableSpots.erase(availableSpots.begin() + gP);
		}
	}
	availableSpots.clear();
	return positions;
}

void Environment::decayAdjuvanicityMap() {
	for (unsigned int i = 0; i < numSpots; ++i)
		if (AdjuMap[i] > 0)
			AdjuMap[i] *= params->adjuDecay;
}

void Environment::modulateFibrosis(std::vector<unsigned int>::iterator bg, std::vector<unsigned int>::iterator ed) {
	for (std::vector<unsigned int>::iterator it = bg; it != ed; it++) {//iterate through seeds
		fibrosis[(*it)] = true;
		if (is3D) {

			//getting coordinates
			int aux = (int)((*it) % (N1*N2));
			int cx = aux % (int)N1;
			int cy = (int)((double)aux / (double)N1);
			int cz = (int)((double)(*it) / (double)(N1*N2));

			//going through neighborhood
			unsigned int pos;
			for (int i = 0; i < smoothMask.size(); ++i) {
				int cxN = cx + smoothMask.at(i).x;
				int cyN = cy + smoothMask.at(i).y;
				int czN = cz + smoothMask.at(i).z;
				if (cxN >= 0 && cxN < (int)N1 && cyN >= 0 && cyN < (int)N2 && czN >= 0 && czN < (int)N3 && randGen->Random() < params->fibrFrac) {
					pos = (unsigned int)czN*N1*N2 + (unsigned int)cyN*N1 + (unsigned int)cxN;
                    if (params->seedUnderneath) {
                        fibrosis[pos] = true;
                    } else if(!lattice[pos]) {//if the lattice spot is empty
                        fibrosis[pos] = true;
                    }
				}

			}


		}
		else {
			//getting coordinates
			int cx = (int)((*it) % N1);
			int cy = (int)((double)((*it)) / (double)N1);

			//going through neighborhood
			unsigned int pos;
			for (int i = 0; i < smoothMask.size(); ++i) {
				int cxN = cx + smoothMask.at(i).x;
				int cyN = cy + smoothMask.at(i).y;
				if (cxN >= 0 && cxN < (int)N1 && cyN >= 0 && cyN < (int)N2 && randGen->Random() < params->fibrFrac) {
					pos = (unsigned int)cyN*N1 + (unsigned int)cxN;
                    if (params->seedUnderneath) {
                        fibrosis[pos] = true;
                    } else if (!lattice[pos]){
                        fibrosis[pos] = true;
                    }
				}

			}
		}
	}
}

void Environment::modulateAdjuvanticity(std::vector<unsigned int>::iterator bg, std::vector<unsigned int>::iterator ed, float addAdju) {
	
	for (std::vector<unsigned int>::iterator it = bg; it != ed; it++) {//iterate through seeds
		if (is3D) {
			//getting coordinates
			int aux = (int)((*it) % (N1*N2));
			int cx = aux % (int)N1;
			int cy = (int)((double)aux / (double)N1);
			int cz = (int)((double)(*it) / (double)(N1*N2));

			//going through neighborhood
			unsigned int pos;
			for (int i = 0; i < adjuMask.size(); ++i) {
				int cxN = cx + adjuMask.at(i).x;
				int cyN = cy + adjuMask.at(i).y;
				int czN = cz + adjuMask.at(i).z;
				if (cxN >= 0 && cxN < (int)N1 && cyN >= 0 && cyN < (int)N2 && czN >= 0 && czN < (int)N3) {
					pos = (unsigned int)czN*N1*N2 + (unsigned int)cyN*N1 + (unsigned int)cxN;
					AdjuMap[pos] += addAdju;
				}

			}

		}
		else {
			//getting coordinates
			int cx = (int)((*it) % N1);
			int cy = (int)((double)((*it)) / (double)N1);

			//going through neighborhood
			unsigned int pos;
			for (int i = 0; i < adjuMask.size(); ++i) {
				int cxN = cx + adjuMask.at(i).x;
				int cyN = cy + adjuMask.at(i).y;
				if (cxN >= 0 && cxN < (int)N1 && cyN >= 0 && cyN < (int)N2) {
					pos = (unsigned int)cyN*N1 + (unsigned int)cxN;
						AdjuMap[pos] += addAdju;
				}

			}

		}
	}
}


unsigned int Environment::returnEmptyPlaceImmune(unsigned int indx, float CErwalk) {
	int nF = 0;
	float maxChemo = 0.0f;

	for (int j = 0; j < numNeigh; j++) {//searching through neighborhood
		if (!lattice[indx - indcNeigh[j]] && (!fibrosis[indx - indcNeigh[j]] || randGen->Random() < params->stromaPerm)) {
			neigh[nF] = indx - indcNeigh[j];
			chemoGrad[nF] = ChtaxMap.getValue(neigh[nF]);
			if (chemoGrad[nF] > maxChemo)
				maxChemo = chemoGrad[nF];
			nF++;
		}
	}
	for (int j = 0; j < numNeigh; j++) {//searching through neighborhood
		if (!lattice[indx + indcNeigh[j]] && (!fibrosis[indx + indcNeigh[j]] || randGen->Random() < params->stromaPerm)) {
			neigh[nF] = indx + indcNeigh[j];
			chemoGrad[nF] = ChtaxMap.getValue(neigh[nF]);
			if (chemoGrad[nF] > maxChemo)
				maxChemo = chemoGrad[nF];
			nF++;
		}
	}

	if (nF) {//selecting free spot taking into account chemotaxis

		if (maxChemo > 0){
			float minChemo = (1.0f - CErwalk)*chemoGrad[0] / maxChemo + CErwalk*(float)randGen->Random();
			int which = 0;
			for (int j = 1; j < nF; j++) {
				float propMinChemo = (1.0f - CErwalk)*chemoGrad[j] / maxChemo + CErwalk*(float)randGen->Random();
				if (propMinChemo > minChemo) {
					which = j;
					minChemo = propMinChemo;
				}
			}

			return neigh[which];
		}
		else {
			return neigh[randGen->Rand(nF)];
		}
	}
	else {//no free spot
		return 0;
	}
}

unsigned int Environment::returnEmptyPlace(unsigned int indx) {
	int nF = 0;
	for (int j = 0; j < numNeigh; j++) {//searching through neighborhood
		if (!lattice[indx - indcNeigh[j]] && (!fibrosis[indx - indcNeigh[j]] || randGen->Random() < params->stromaPerm) ) {
			neigh[nF] = indx - indcNeigh[j];
			nF++;
		}
	}
	for (int j = 0; j < numNeigh; j++) {//searching through neighborhood
		if (!lattice[indx + indcNeigh[j]] && (!fibrosis[indx + indcNeigh[j]] || randGen->Random() < params->stromaPerm)) {
			neigh[nF] = indx + indcNeigh[j];
			nF++;
		}
	}
	if (nF) {//selecting free spot at random
		return neigh[randGen->Rand(nF)];
	}
	else {//no free spot
		return 0;
	}
}

unsigned int Environment::findTarget(unsigned int indx) {
	int nF = 0;
	for (int j = 0; j < numNeigh; j++) {//searching through neighborhood
		if (latticeOnlyTU[indx - indcNeigh[j]]) {
			neigh[nF] = indx - indcNeigh[j];
			nF++;
		}
	}
	for (int j = 0; j < numNeigh; j++) {//searching through neighborhood
		if (latticeOnlyTU[indx + indcNeigh[j]]) {
			neigh[nF] = indx + indcNeigh[j];
			nF++;
		}
	}
	if (nF) {//selecting free spot at random
		return neigh[randGen->Rand(nF)];
	}
	else {//no free spot
		return 0;
	}
}


void Environment::initialize(SIMparameters* parIn, CRandomMersenne *rG) {
	//setting enviroment dimension
	params = parIn;
	N1 = params->N1; N2 = params->N2; N3 = params->N3;
	numSpots = N1*N2*N3;
	is3D = (N3 > 1);

	numOccupiedSpots = 0;
	numNecroticSpots = 0;
	numFibroticSpots = 0;

	randGen = rG;

	if (lattice != NULL)
		delete[] lattice;
	if (necrosis != NULL)
		delete[] necrosis;
	if (fibrosis != NULL)
		delete[] fibrosis;
	if (AdjuMap != NULL)
		delete[] AdjuMap;
	if (latticeOnlyTU != NULL)
		delete[] latticeOnlyTU;

	adjuMask.clear();
	smoothMask.clear();

	lattice = DBG_NEW bool[numSpots];
	memset(lattice, 0, numSpots * sizeof(bool));
	necrosis = DBG_NEW bool[numSpots];
	memset(necrosis, 0, numSpots * sizeof(bool));
	fibrosis = DBG_NEW bool[numSpots];
	memset(fibrosis, 0, numSpots * sizeof(bool));
	latticeOnlyTU = DBG_NEW bool[numSpots];
	memset(latticeOnlyTU, 0, numSpots * sizeof(bool));

	AdjuMap = DBG_NEW float[numSpots];
	memset(AdjuMap, 0, numSpots * sizeof(float));

	ChtaxMap.initialize(params);
	NecroMap.initialize(params);

	if (is3D) {//3D implementation
		//setting lattice boundary
		for (unsigned int i = 0; i < N1; ++i)
			for (unsigned int j = 0; j < N2; ++j)
				for (unsigned int k = 0; k < N3; ++k)
					if (i == 0 || i == N1 - 1 || j == 0 || j == N2 - 1 || k == 0 || k == N3 - 1) {
						lattice[k*N1*N2 + j*N1 + i] = true;
						numOccupiedSpots++;
					}

		//defining the neighborhood
		if (indcNeigh != NULL)
			delete[] indcNeigh;
		indcNeigh = DBG_NEW unsigned int[13];
		int indx, aux = 0;
		for (int i = -1; i < 2; ++i)
			for (int j = -1; j < 2; ++j)
				for (int k = -1; k < 2; ++k) {
					indx = k*N1*N2+j*N1+i;
					if (indx > 0) {
						indcNeigh[aux] = (unsigned int)indx;
						aux++;
					}
						
				}

		
		if (neigh != NULL)
			delete[] neigh;
		neigh = DBG_NEW unsigned int[26];
		numNeigh = 13;

		if (chemoGrad != NULL)
			delete[] chemoGrad;
		chemoGrad = DBG_NEW float[26];
		//defining adju mask
		int r2 = params->adjuRange*params->adjuRange;
		for (int i = -params->adjuRange; i <= params->adjuRange; ++i)
			for (int j = -params->adjuRange; j <= params->adjuRange; ++j)
				for (int k = -params->adjuRange; k <= params->adjuRange; ++k)
				if (i*i + j*j + k*k <= r2) {//within the sphere
					point point = { i, j, k };
					adjuMask.push_back(point);
				}
		//defining smoothMask
		r2 = params->smoothRadius*params->smoothRadius;
		for (int i = -params->smoothRadius; i <= params->smoothRadius; ++i)
			for (int j = -params->smoothRadius; j <= params->smoothRadius; ++j)
				for (int k = -params->smoothRadius; k <= params->smoothRadius; ++k)
					if (i*i + j*j + k*k <= r2) {//within the sphere
						point point = { i, j, k };
						smoothMask.push_back(point);
					}


	}
	else {//2D case

		//setting lattice boundary
		unsigned int i;
		for (i = 0; i < N1; i++) { lattice[i] = true; numOccupiedSpots++; };//left
		for (i = 0; i < N1*N2; i = i + N1) { lattice[i] = true; numOccupiedSpots++;};//top
		for (i = N1 - 1; i < N1*N2; i = i + N1) { lattice[i] = true; numOccupiedSpots++;};//bottom
		for (i = N1*(N2 - 1); i < N1*N2; i++) { lattice[i] = true;  numOccupiedSpots++;};//right

		//defining the neighborhood
		if (indcNeigh != NULL)
			delete[] indcNeigh;

		indcNeigh = DBG_NEW unsigned int[4];
		indcNeigh[0] = N1 + 1; //right-bottom
		indcNeigh[1] = N1; //right
		indcNeigh[2] = N1 - 1; //right-top
		indcNeigh[3] = 1; //bottom

		if (neigh != NULL)
			delete[] neigh;
		neigh = DBG_NEW unsigned int[8];
		numNeigh = 4; //half of the actual number

		if (chemoGrad != NULL)
			delete[] chemoGrad;
		chemoGrad = DBG_NEW float[8];

		//defining adju mask
		int r2 = params->adjuRange*params->adjuRange;
		for (int i = -params->adjuRange; i <= params->adjuRange; ++i)
			for (int j = -params->adjuRange; j<=params->adjuRange; ++j)
				if (i*i + j*j <= r2) {//within the circle
					point point = { i, j, 0 };
					adjuMask.push_back(point);
				}

		//defining smoothMask
		r2 = params->smoothRadius*params->smoothRadius;
		for (int i = -params->smoothRadius; i <= params->smoothRadius; ++i)
			for (int j = -params->smoothRadius; j <= params->smoothRadius; ++j)
				if (i*i + j*j <= r2) {//within the circle
					point point = { i, j, 0 };
					smoothMask.push_back(point);
				}
					
	}

}

void Environment::readState(const mxArray* sysTempl) {
	mxArray *grid = mxGetField(sysTempl, 0, "grid");
	float* AdjuMapS = (float*)mxGetPr(mxGetField(grid, 0, "AdjuMap"));
	memcpy(AdjuMap, AdjuMapS, numSpots * sizeof(float));
	bool* Lf = (bool*)mxGetPr(mxGetField(grid, 0, "Lf"));
	memcpy(fibrosis, Lf, numSpots * sizeof(bool));

}

mxArray* Environment::getState() {

	const char *field_names[] = { "L", "Ln", "Lf", "AdjuMap", "ChtaxMap","NecroMap"};
	mwSize dimsStruct[2] = { 1, 1 };
	mxArray* envState = mxCreateStructArray(1, dimsStruct, 6, field_names);

	
	//2 - copying the L matrix
	mwSize dims[3] = { (mwSize)N1, (mwSize)N2, (mwSize)N3 };
	mxSetField(envState, 0, "L",mxCreateLogicalArray(3, dims));
	bool *L = (bool*)mxGetPr(mxGetField(envState, 0, "L"));
	memcpy(L, lattice, numSpots* sizeof(bool));

	mxSetField(envState, 0, "Ln", mxCreateLogicalArray(3, dims));
	bool *Ln = (bool*)mxGetPr(mxGetField(envState, 0, "Ln"));
	memcpy(Ln, necrosis, numSpots * sizeof(bool));

	mxSetField(envState, 0, "Lf", mxCreateLogicalArray(3, dims));
	bool *Lf = (bool*)mxGetPr(mxGetField(envState, 0, "Lf"));
	memcpy(Lf, fibrosis, numSpots * sizeof(bool));

	mxSetField(envState, 0, "AdjuMap", mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL));
	float *Amap = (float*)mxGetPr(mxGetField(envState, 0, "AdjuMap"));
	memcpy(Amap, AdjuMap, numSpots * sizeof(float));

	mxSetField(envState, 0, "ChtaxMap", ChtaxMap.getState());
	mxSetField(envState, 0, "NecroMap", NecroMap.getState());
	return envState;
	
}
