/*
*  Created by Jan Poleszczuk
*  Last modified November, 2017 by Jan Poleszczuk
*/

#include "stdafx.h"
#include "NecrosisMapSuiteS.h"


NecrosisMap::NecrosisMap() {
	//setting default lattice sizes
	N1 = 100;
	N2 = 100;
	N3 = 1;
	is3D = false;
	numSpots = N1*N2;
	afterInit = false;
};

NecrosisMap::~NecrosisMap() {
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
		delete[] changes;
		afterInit = false;
	}
};



void NecrosisMap::initialize(SIMparameters* parIn) {
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
			coefficients.push_back(T(i, i, params->DCnecro*6.));//0 diagonal
			if ((i + 1) % N1D)
				coefficients.push_back(T(i + 1, i, -params->DCnecro));//-1 diagonal
			if (i % N1D)
				coefficients.push_back(T(i - 1, i, -params->DCnecro));//1 diagonal
			if ((i % (N1D*N2D)) + N1D < N1D*N2D)
				coefficients.push_back(T(i + N1D, i, -params->DCnecro));//-N1 diagonal
			if (((int)(i % (N1D*N2D)) - (int)N1D) >= 0)
				coefficients.push_back(T(i - N1D, i, -params->DCnecro));//N1 diagonal
			if (i + N1D*N2D < numSpotsD)
				coefficients.push_back(T(i + N1D*N2D, i, -params->DCnecro));//-N1*N2 diagonal
			if (((int)i - (int)(N1D*N2D)) >= 0)
				coefficients.push_back(T(i - N1D*N2D, i, -params->DCnecro));//N1*N2 diagonal
		}


		Snew.setFromTriplets(coefficients.begin(), coefficients.end());
		Snew.makeCompressed();
		RmS3D = Snew;

		b3D.resize(numSpotsD);
		b3D.fill(0.);
		//vector b
		for (unsigned int i = 0; i < N1D; ++i)
			for (unsigned int j = 0; j < N2D; ++j)
				for (unsigned int k = 0; k < N3D; ++k)
					if (i == 0 || i == N1D - 1 || j == 0 || j == N2D - 1 || k == 0 || k == N3D - 1) {
						char numFS = (char)(i == 0) + (char)(i == N1D - 1) + (char)(j == 0) + (char)(j == N2D - 1) + (char)(k == 0) + (char)(k == N3D - 1);
						b3D(k*N1D*N2D + j*N1D + i) = (double)numFS*params->DCnecro;
					}


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
		state3D.fill(1.); //initial state

	} else {//is 2D, using SuiteSparse

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
			Ti[k] = i; Tj[k] = i; Tx[k] = params->DCnecro*4.; //0 diagonal
			k++;
			if (i % N1D) {
				Ti[k] = i - 1; Tj[k] = i; Tx[k] = -params->DCnecro;
				k++;
			}
			if (((int)i - (int)N1D) >= 0){
				Ti[k] = i - N1D; Tj[k] = i; Tx[k] = -params->DCnecro;
				k++;
			}
		}
		T->nnz = k;

		A = cholmod_triplet_to_sparse(T, 0, &c);

		cholmod_free_triplet(&T, &c);

		//cholmod_print_sparse(A, "A", &c);

		L = cholmod_analyze(A, &c);
		cholmod_factorize(A, L, &c);


		b = cholmod_zeros(numSpotsD, 1, A->xtype, &c);
		double *bx = (double*)b->x;
		for (unsigned int i = 0; i < N1D; ++i)
			bx[i] += params->DCnecro;
		for (unsigned int i = 0; i < numSpotsD; i = i + N1D)
			bx[i] += params->DCnecro;
		for (unsigned int i = (N2D - 1)*N1D; i < numSpotsD; ++i)
			bx[i] += params->DCnecro;
		for (unsigned int i = N1D - 1; i < numSpotsD; i = i + N1D)
			bx[i] += params->DCnecro;

		//cholmod_print_dense(b, "b", &c);

		state = cholmod_ones(numSpots, 1, A->xtype, &c);

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


		bInterp = cholmod_zeros(numSpots, 1, A->xtype, &c);
		bx = (double*)bInterp->x;
		for (unsigned int i = 0; i < N1; ++i)
			for (unsigned int j = 0; j < N2; ++j)
				if (i == 0 || i == N1 - 1 || j == 0 || j == N2 - 1) {
					bx[j*N1 + i] = 1. / 3.;
				}
	}//end is 3D
};

void NecrosisMap::addSink(unsigned int src) {
	changes[src] += 1;
}

void NecrosisMap::removeSink(unsigned int src) {
	changes[src] -= 1;
}

void NecrosisMap::updateMap() {

	char *iCC = DBG_NEW char[numSpotsD];
	memset(iCC, 0, numSpotsD*sizeof(char));
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

	int nnzU = 0, nnzD = 0;

	for (unsigned int i = 0; i < numSpotsD; i++)
		if (iCC[i]>0){
			nnzU++;
		}
		else if (iCC[i] < 0) {
			nnzD++;
		}
		
		if (is3D) {
			if (nnzU || nnzD) {

				for (unsigned int i = 0; i < numSpotsD; ++i)
					if (iCC[i]!= 0) { //update
						RmS3D.coeffRef(i, i) += iCC[i]*params->TCnecro;
					}

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
		else {//is 2D
			if (nnzU || nnzD) {//if any change to be made
				cholmod_triplet *Tu = cholmod_allocate_triplet(numSpotsD, nnzU, nnzU, 0, CHOLMOD_REAL, &c); //unsymmetric, permuttation does not support symmetric upper/lower
				Tu->nnz = nnzU;
				int *Tiu = (int*)Tu->i, *Tju = (int*)Tu->j;
				double *Txu = (double *)Tu->x;

				cholmod_triplet *Td = cholmod_allocate_triplet(numSpotsD, nnzD, nnzD, 0, CHOLMOD_REAL, &c); //unsymmetric, permuttation does not support symmetric upper/lower
				Td->nnz = nnzD;
				int *Tid = (int*)Td->i, *Tjd = (int*)Td->j;
				double *Txd = (double *)Td->x;


				int ku = 0, kd = 0;
				for (unsigned int i = 0; i < numSpotsD; ++i)
					if (iCC[i] > 0) { //update
						Tiu[ku] = i; Tju[ku] = ku; Txu[ku] = sqrt((double)iCC[i]*params->TCnecro);
						ku++;
					}
					else if (iCC[i] < 0) {//downdate
						Tid[kd] = i; Tjd[kd] = kd; Txd[kd] = sqrt(-(double)iCC[i]*params->TCnecro);
						kd++;
					}

					cholmod_sparse *C, *CP;
					if (nnzU) {
						C = cholmod_triplet_to_sparse(Tu, 0, &c);
						CP = cholmod_submatrix(C, (int*)L->Perm, L->n, NULL, -1, TRUE, TRUE, &c);
						cholmod_updown(1, CP, L, &c);
						cholmod_free_sparse(&C, &c);
						cholmod_free_sparse(&CP, &c);
					}

					if (nnzD) {
						C = cholmod_triplet_to_sparse(Td, 0, &c);
						CP = cholmod_submatrix(C, (int*)L->Perm, L->n, NULL, -1, TRUE, TRUE, &c);
						cholmod_updown(0, CP, L, &c);
						cholmod_free_sparse(&C, &c);
						cholmod_free_sparse(&CP, &c);
					}

					cholmod_free_triplet(&Td, &c);
					cholmod_free_triplet(&Tu, &c);

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


mxArray* NecrosisMap::getState() {

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
