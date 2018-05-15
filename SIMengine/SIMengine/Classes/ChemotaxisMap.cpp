/*
*  Created by Jan Poleszczuk
*  Last modified November, 2017 by Jan Poleszczuk
*/

#include "stdafx.h"
#include "ChemotaxisMap.h"

ChemotaxisMap::ChemotaxisMap() {
	map = NULL;
	kernel = NULL;
	changes = NULL;

	//setting default lattice sizes
	N1 = 100;
	N2 = 100;
	N3 = 1;
	is3D = false;
	numSpots = 10000;
};

ChemotaxisMap::~ChemotaxisMap() {
	if (map != NULL)
		delete[] map;
	if (kernel != NULL)
		delete[] kernel;
	if (changes != NULL)
		delete[] changes;
};

void ChemotaxisMap::initialize(SIMparameters* parIn) {
	params = parIn;

	N1 = params->N1; N2 = params->N2; N3 = params->N3;
	numSpots = N1*N2*N3;
	is3D = (N3 > 1);

	if (map != NULL)
		delete[] map;

	map = DBG_NEW float[numSpots];
	memset(map, 0, numSpots * sizeof(float));

	//initializing diffusion kernel
	if (kernel != NULL)
		delete[] kernel;

	kernelDim = (int)params->radCutoff;
	kernelWidth = (2 * kernelDim + 1);
	
	if (is3D)
		kernelNum = kernelWidth*kernelWidth*kernelWidth;
	else
		kernelNum = kernelWidth*kernelWidth;

	kernel = DBG_NEW float[kernelNum];
	memset(kernel, 0, kernelNum*sizeof(float));

	//evaluating shift
	double r = params->radCutoff + 1.;
	float shift = (float)(exp(-params->effDiff *r) / r);
	int pos = 0;
	
	if (is3D) {
		for (int i = -kernelDim; i <= kernelDim; ++i)
			for (int j = -kernelDim; j <= kernelDim; ++j)
				for (int k = -kernelDim; k <= kernelDim; ++k) {
				r = sqrt((double)i*(double)i + (double)j*(double)j + (double)k*(double)k);
				if (r == 0.0)
					r = 1.0;
				if (r <= params->radCutoff)
					kernel[pos] = (float)(exp(-params->effDiff *r) / r) - shift;
				pos++;
			}
	}
	else {
		for (int i = -kernelDim; i <= kernelDim; ++i)
			for (int j = -kernelDim; j <= kernelDim; ++j) {
				r = sqrt((double)i*(double)i + (double)j*(double)j);
				if (r == 0.0)
					r = 1.0;
				if (r <= params->radCutoff)
					kernel[pos] = (float)(exp(-params->effDiff *r) / r) - shift;
				pos++;
			}
	}


	if (changes != NULL)
		delete[] changes;

	changes = DBG_NEW char[numSpots];
	memset(changes, 0, numSpots * sizeof(char));

};

void ChemotaxisMap::addSource(unsigned int src) {
	changes[src] += 1;
}

void ChemotaxisMap::removeSource(unsigned int src) {
	changes[src] -= 1;
}

void ChemotaxisMap::updateMap() {
	int aux, cx, cy, cz;
	int auxX, auxY, auxZ;
	int begMapX, begMapY, begMapZ;
	int endMapX, endMapY, endMapZ;
	int begKernelX, begKernelY, begKernelZ;

	int widthX, widthY, widthZ;
	unsigned int shiftM, shiftK, posK, posM;
	unsigned int shiftM2, shiftK2;
	bool add;

	for (unsigned int i = 0; i < numSpots; ++i) {
		
		if (((int)changes[i]) != 0) {
			add = ((int)changes[i]) > 0;

			aux = (int)(i % (N1*N2));
			cx = aux % (int)N1;
			cy = (int)((double)aux / (double)N1);
			cz = (int)((double)i / (double)(N1*N2));
			
			if (is3D) {
				auxX = cx - kernelDim;
				auxY = cy - kernelDim;
				auxZ = cz - kernelDim;

				begMapX = std::max(auxX, 0);
				begMapY = std::max(auxY, 0);
				begMapZ = std::max(auxZ, 0);
				endMapX = std::min(cx + kernelDim, (int)N1 - 1);
				endMapY = std::min(cy + kernelDim, (int)N2 - 1);
				endMapZ = std::min(cz + kernelDim, (int)N3 - 1);

				begKernelX = std::max(-auxX, 0);
				begKernelY = std::max(-auxY, 0);
				begKernelZ = std::max(-auxZ, 0);


				widthZ = endMapZ - begMapZ;
				widthY = endMapY - begMapY;
				widthX = endMapX - begMapX;

				shiftM = N1 - (unsigned int)widthX - 1;
				shiftK = (unsigned int)(kernelWidth - widthX - 1);
				shiftM2 = N1*N2 - (widthY + 1)*shiftM - (widthX + 1)*(widthY + 1);
				shiftK2 = (unsigned int)kernelWidth*(unsigned int)kernelWidth - (widthY+1)*shiftK - (widthX+1)*(widthY+1);
				

				posK = (unsigned int)(begKernelZ)*(unsigned int)kernelWidth*(unsigned int)kernelWidth + (unsigned int)begKernelY*(unsigned int)kernelWidth + (unsigned int)begKernelX;
				posM = (unsigned int)(begMapZ)*N1*N2 + (unsigned int)begMapY*N1 + (unsigned int)begMapX;

				for (int k = 0; k <= widthZ; ++k) {

					for (int j = 0; j <= widthY; ++j) {
						for (int i = 0; i <= widthX; ++i) {
							if (add) {
								map[posM] += kernel[posK];
							}
							else {
								map[posM] -= kernel[posK];
							}
							posK++;
							posM++;
						}
						posK += shiftK;
						posM += shiftM;
					}

					posK += shiftK2;//(unsigned int)(begKernelZ + k)*(unsigned int)kernelWidth*(unsigned int)kernelWidth + (unsigned int)begKernelY*(unsigned int)kernelWidth + (unsigned int)begKernelX;
					posM += shiftM2;//(unsigned int)(begMapZ + k)*N1*N2 + (unsigned int)begMapY*N1 + (unsigned int)begMapX;
				}
			}
			else {//is 2D
				auxX = cx - kernelDim;
				auxY = cy - kernelDim;
				begMapX = std::max(auxX, 0);
				begMapY = std::max(auxY, 0);
				endMapX = std::min(cx + kernelDim, (int)N1 - 1);
				endMapY = std::min(cy + kernelDim, (int)N2 - 1);

				begKernelX = std::max(-auxX, 0);
				begKernelY = std::max(-auxY, 0);

				posK = (unsigned int)begKernelY*(unsigned int)kernelWidth + (unsigned int)begKernelX;
				posM = (unsigned int)begMapY*N1 + (unsigned int)begMapX;

				widthY = endMapY - begMapY;
				widthX = endMapX - begMapX;

				shiftM = N1 - (unsigned int)widthX - 1;
				shiftK = (unsigned int)(kernelWidth - widthX - 1);

				
				for (int j = 0; j <= widthY; ++j) {
					for (int i = 0; i <= widthX; ++i) {
						if (add) {
							map[posM] += kernel[posK];
						}
						else {
							map[posM] -= kernel[posK];
						}
						posK++;
						posM++;
					}
					posK += shiftK;
					posM += shiftM;
				}
			}//end is 3D

		}//end changes[i]!=0
	}
	
	memset(changes, 0, numSpots * sizeof(char));
}

mxArray* ChemotaxisMap::getState() {

	mwSize dims[3] = { (mwSize)N1, (mwSize)N2, (mwSize)N3 };
	mxArray *outMap = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
	float *ptr = (float*)mxGetPr(outMap);
	memcpy(ptr, map, numSpots*sizeof(float));

	return outMap;
}