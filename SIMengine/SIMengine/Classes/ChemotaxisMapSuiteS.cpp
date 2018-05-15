/*
*  Created by Jan Poleszczuk
*  Last modified November, 2017 by Jan Poleszczuk
*/

#include "stdafx.h"
#include "ChemotaxisMapSuiteS.h"

ChemotaxisMap::ChemotaxisMap() {
	//setting default lattice sizes
	N1 = 100;
	N2 = 100;
	N3 = 1;
	is3D = false;
	numSpots = N1*N2;
	afterInit = false;
};

ChemotaxisMap::~ChemotaxisMap() {
	if (afterInit) {
		if (!is3D) {
			cholmod_free_sparse(&A, &c);
			cholmod_free_sparse(&Interp, &c);
			cholmod_free_dense(&state, &c);
			cholmod_free_dense(&b, &c);
			cholmod_free_dense(&bInterp, &c);
			cholmod_free_factor(&L, &c);  /* free matrices */
			cholmod_finish(&c); /* finish CHOLMOD */
		}
		afterInit = false;
		delete[] changes;
	}
};

void ChemotaxisMap::initialize(SIMparameters* parIn) {
	params = parIn;

	N1 = params->N1; N2 = params->N2; N3 = params->N3;
	numSpots = N1*N2*N3;
	is3D = (N3 > 1);

	N1D = N1 / 3; N2D = N2 / 3;

	N3D = is3D ? N3 / 3 : 1;

	numSpotsD = N1D*N2D*N3D;

	if (is3D) {
		if (afterInit) {
			delete[] changes;
		}
		afterInit = true;

		changes = DBG_NEW char[numSpots];
		memset(changes, 0, numSpots * sizeof(char));

		SpMat Snew(numSpotsD, numSpotsD);
		std::vector<T> coefficients;
		coefficients.reserve(7 * numSpotsD);
		for (unsigned int i = 0; i < numSpotsD; ++i) {
			coefficients.push_back(T(i, i, 1.+params->DCchemo*6.));//0 diagonal
			if ((i + 1) % N1D)
				coefficients.push_back(T(i + 1, i, -params->DCchemo));//-1 diagonal
			if (i % N1D)
				coefficients.push_back(T(i - 1, i, -params->DCchemo));//1 diagonal
			if ((i % (N1D*N2D)) + N1D < N1D*N2D)
				coefficients.push_back(T(i + N1D, i, -params->DCchemo));//-N1 diagonal
			if (((int)(i % (N1D*N2D)) - (int)N1D) >= 0)
				coefficients.push_back(T(i - N1D, i, -params->DCchemo));//N1 diagonal
			if (i + N1D*N2D < numSpotsD)
				coefficients.push_back(T(i + N1D*N2D, i, -params->DCchemo));//-N1*N2 diagonal
			if (((int)i - (int)(N1D*N2D)) >= 0)
				coefficients.push_back(T(i - N1D*N2D, i, -params->DCchemo));//N1*N2 diagonal
		}


		Snew.setFromTriplets(coefficients.begin(), coefficients.end());
		Snew.makeCompressed();
		RmS3D = Snew;

		b3D.resize(numSpotsD);
		b3D.fill(0.);
		

		//initializing Interpolation matrix
		coefficients.clear();
		SpMat InterpTmp(numSpots, numSpotsD);
		int aux, cx, cy, cz, whX, whY, whZ, rPx, rPy, rPz, whX2, whY2, whZ2;
		for (unsigned int i = 0; i < numSpots; ++i) {
			//casting to super node
			aux = (int)(i % (N1*N2));
			cx = aux % (int)N1;
			cy = (int)((double)aux / (double)N1);
			cz = (int)((double)i / (double)(N1*N2));

			whX = cx / 3;
			whY = cy / 3;
			whZ = cz / 3;

			rPx = (cx % 3) - 1;
			rPy = (cy % 3) - 1;
			rPz = (cz % 3) - 1;

			if (rPx != 0 || rPy != 0 || rPz != 0) {
				whX2 = whX + rPx;
				whY2 = whY + rPy;
				whZ2 = whZ + rPz;

				coefficients.push_back(T(i, whX + whY*N1D + whZ*N1D*N2D, 2. / 3.));//closer super node

				if (whX2 >= 0 && whX2 < (int)N1D && whY2 >= 0 && whY2 < (int)N2D && whZ2 >= 0 && whZ2 < (int)N3D)
					coefficients.push_back(T(i, whX2 + whY2*N1D + whZ2*N1D*N2D, 1. / 3.));//other node
			}
			else {
				coefficients.push_back(T(i, whX + whY*N1D + whZ*N1D*N2D, 1.));//center super node
			}
		}

		InterpTmp.setFromTriplets(coefficients.begin(), coefficients.end());
		InterpTmp.makeCompressed();
		Interp3D = InterpTmp;


		bInterp3D.resize(numSpots);
		bInterp3D.fill(0.);

		for (unsigned int i = 0; i < N1; ++i)
			for (unsigned int j = 0; j < N2; ++j)
				for (unsigned int k = 0; k < N3; ++k)
					if (i == 0 || i == N1 - 1 || j == 0 || j == N2 - 1 || k == 0 || k == N3 - 1) {
						bInterp3D(k*N1*N2 + j*N1 + i) = 1. / 3.;
					}
		state3D.resize(numSpots);
		state3D.fill(0.); //initial state

	}
	else {
		cholmod_start(&c); /* start CHOLMOD */

		if (afterInit) {//already initialized -> clear the workspace
			cholmod_free_sparse(&A, &c);
			cholmod_free_sparse(&Interp, &c);
			cholmod_free_dense(&state, &c);
			cholmod_free_dense(&b, &c);
			cholmod_free_dense(&bInterp, &c);
			cholmod_free_factor(&L, &c);  /* free matrices */
			cholmod_finish(&c); /* finish CHOLMOD */
			delete[] changes;
		}
		afterInit = true;

		changes = DBG_NEW char[numSpots];
		memset(changes, 0, numSpots * sizeof(char));

		
		cholmod_triplet *T = cholmod_allocate_triplet(numSpotsD, numSpotsD, 3 * numSpotsD, 1, CHOLMOD_REAL, &c); //symmetric upper

		int *Ti = (int*)T->i, *Tj = (int*)T->j;
		double *Tx = (double *)T->x;
		int k = 0;
		for (int i = 0; i < (int)numSpotsD; ++i) {
			Ti[k] = i; Tj[k] = i; Tx[k] = 1. + params->DCchemo*4.; //0 diagonal
			k++;
			if (i % N1D) {
				Ti[k] = i - 1; Tj[k] = i; Tx[k] = -params->DCchemo;
				k++;
			}
			if (((int)i - (int)N1D) >= 0){
				Ti[k] = i - N1D; Tj[k] = i; Tx[k] = -params->DCchemo;
				k++;
			}
		}
		T->nnz = k;

		A = cholmod_triplet_to_sparse(T, 0, &c);

		cholmod_free_triplet(&T, &c);

		
		L = cholmod_analyze(A, &c);
		cholmod_factorize(A, L, &c);


		b = cholmod_zeros(numSpotsD, 1, A->xtype, &c);

		state = cholmod_zeros(numSpots, 1, A->xtype, &c);

		
		T = cholmod_allocate_triplet(numSpots, numSpotsD, 2 * numSpots, 0, CHOLMOD_REAL, &c); //unsymmetric
		Ti = (int*)T->i;
		Tj = (int*)T->j;
		Tx = (double *)T->x;

		int cx, cy, whX, whY, rPx, rPy, whX2, whY2;
		k = 0;
		for (unsigned int i = 0; i < numSpots; ++i) {
			//casting to super node
			cx = i % (int)N1;
			cy = (int)((double)i / (double)N1);

			whX = cx / 3;
			whY = cy / 3;

			rPx = (cx % 3) - 1;
			rPy = (cy % 3) - 1;


			if (rPx != 0 || rPy != 0) {
				whX2 = whX + rPx;
				whY2 = whY + rPy;

				Ti[k] = i; Tj[k] = whX + whY*N1D; Tx[k] = 2. / 3.;//center super node
				k++;
				if (whX2 >= 0 && whX2 < (int)N1D && whY2 >= 0 && whY2 < (int)N2D) {
					Ti[k] = i; Tj[k] = whX2 + whY2*N1D; Tx[k] = 1. / 3.;//other node
					k++;
				}
			}
			else {
				Ti[k] = i; Tj[k] = whX + whY*N1D; Tx[k] = 1.;//center super node
				k++;
			}
		}
		T->nnz = k;
		Interp = cholmod_triplet_to_sparse(T, 0, &c);
		cholmod_free_triplet(&T, &c);

		//cholmod_print_sparse(Interp, "Interp", &c);

		bInterp = cholmod_zeros(numSpots, 1, A->xtype, &c);
		double *bx = (double*)bInterp->x;
		for (unsigned int i = 0; i < N1; ++i)
			for (unsigned int j = 0; j < N2; ++j)
				if (i == 0 || i == N1 - 1 || j == 0 || j == N2 - 1) {
					bx[j*N1 + i] = 1. / 3.;
				}
				
	}//end is 3D
};

void ChemotaxisMap::addSource(unsigned int src) {
	changes[src] += 1;
}

void ChemotaxisMap::removeSource(unsigned int src) {
	changes[src] -= 1;
}

void ChemotaxisMap::updateMap() {
	int *iCC = DBG_NEW int[numSpotsD];
	memset(iCC, 0, numSpotsD*sizeof(int));
	int iC;
	int aux, cx, cy, cz;

	for (unsigned int i = 0; i < numSpots; ++i)
		if (changes[i] == 1) {
			aux = (int)(i % (N1*N2));
			cx = aux % (int)N1;
			cy = (int)((double)aux / (double)N1);
			cz = (int)((double)i / (double)(N1*N2));
			iC = (int)((cx / 3) + (cy / 3)*N1D + (cz / 3)*(N1D*N2D));
			iCC[iC]++;
		}
		else if (changes[i] == -1) {
			aux = (int)(i % (N1*N2));
			cx = aux % (int)N1;
			cy = (int)((double)aux / (double)N1);
			cz = (int)((double)i / (double)(N1*N2));
			iC = (cx / 3) + (cy / 3)*N1D + (cz / 3)*(N1D*N2D);
			iCC[iC]--;
		}

		
		
		if (is3D) {
			int nnz = 0;
			for (unsigned int i = 0; i < numSpotsD; i++)
				if (iCC[i] != 0){
					b3D(i) += ((double)iCC[i])*params->SCchemo;
					nnz++;
				}

			if (nnz) {

				solver3D.compute(RmS3D);
				if (solver3D.info() != Eigen::Success) {
					printf("NecrosisMap -> Compute step failed! \n");
					return;
				}
				Eigen::VectorXd x = solver3D.solve(b3D);
				if (solver3D.info() != Eigen::Success) {
					printf("Error: NecrosisMap -> Solving failed! \n");
					return;
				}
				//interpolating the state
				state3D = Interp3D*x + bInterp3D;
			}
		}
		else {
			int nnz = 0;
			double *bx = (double*)b->x;
			for (unsigned int i = 0; i < numSpotsD; i++)
				if (iCC[i] != 0){
					bx[i] += ((double)iCC[i])*params->SCchemo;
					nnz++;
				}

			if (nnz) {
				cholmod_dense *x = cholmod_solve(CHOLMOD_A, L, b, &c);


				double one[2] = { 1, 0 }, zero[2] = { 0, 0 };
				cholmod_sdmult(Interp, 0, one, zero, x, state, &c);
				cholmod_free_dense(&x, &c);

				for (int i = 0; i < state->nrow; ++i)
					((double*)state->x)[i] += ((double*)bInterp->x)[i];
			}
		}

				delete[]  iCC;
				//clearing changes vector
				memset(changes, 0, numSpots*sizeof(char));
}

mxArray* ChemotaxisMap::getState() {

	mwSize dims[3] = { (mwSize)N1, (mwSize)N2, (mwSize)N3 };
	mxArray *outMap = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
	double *ptr = mxGetPr(outMap);
	if (is3D){
		memcpy(ptr, state3D.data(), numSpots*sizeof(double));
	}
	else {
		memcpy(ptr, state->x, numSpots*sizeof(double));
	}

	return outMap;
}