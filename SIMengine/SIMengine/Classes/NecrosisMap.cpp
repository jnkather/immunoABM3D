/*
*  Created by Jan Poleszczuk
*  Last modified November, 2017 by Jan Poleszczuk
*/

#include "stdafx.h"
#include "NecrosisMap.h"

NecrosisMap::NecrosisMap() {
	//setting default lattice sizes
	N1 = 100;
	N2 = 100;
	N3 = 1;
	is3D = false;
	numSpots = N1*N2;
};

NecrosisMap::~NecrosisMap() {

};



void NecrosisMap::initialize(SIMparameters* parIn) {
	params = parIn;

	N1 = params->N1; N2 = params->N2; N3 = params->N3;
	numSpots = N1*N2*N3;
	is3D = (N3 > 1);
	
	N1D = N1 / 3; N2D = N2 / 3; 
	
	N3D = is3D ? N3 / 3 : 1;
	
	numSpotsD = N1D*N2D*N3D;

	//building the problem
	SpMat Snew(numSpotsD, numSpotsD);

	if (is3D){
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
		RmS = Snew;
				
		b.resize(numSpotsD);
		b.fill(0.);
		//vector b
		for (unsigned int i = 0; i < N1D; ++i)
			for (unsigned int j = 0; j < N2D; ++j)
				for (unsigned int k = 0; k < N3D; ++k)
					if (i == 0 || i == N1D - 1 || j == 0 || j == N2D - 1 || k == 0 || k == N3D - 1) {
						char numFS = (char)(i == 0) + (char)(i == N1D - 1) + (char)(j == 0) + (char)(j == N2D - 1) + (char)(k == 0) + (char)(k == N3D - 1);
						b(k*N1D*N2D + j*N1D + i) = (double)numFS*params->DCnecro;
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
		Interp = InterpTmp;

		bInterp.resize(numSpots);
		bInterp.fill(0.);

		for (unsigned int i = 0; i < N1; ++i)
			for (unsigned int j = 0; j < N2; ++j)
				for (unsigned int k = 0; k < N3; ++k)
					if (i == 0 || i == N1 - 1 || j == 0 || j == N2 - 1 || k == 0 || k == N3 - 1) {
						bInterp(k*N1*N2 + j*N1 + i) = 1. / 3.;
					}

	}
	else {
		std::vector<T> coefficients;
		coefficients.reserve(5 * numSpotsD);
		for (unsigned int i = 0; i < numSpotsD; ++i) {
			coefficients.push_back(T(i, i, params->DCnecro*4.));//0 diagonal
			if ((i + 1) % N1D)
				coefficients.push_back(T(i + 1, i, -params->DCnecro));//-1 diagonal
			if (i % N1D)
				coefficients.push_back(T(i - 1, i, -params->DCnecro));//1 diagonal
			if (i + N1D < numSpotsD)
				coefficients.push_back(T(i + N1D, i, -params->DCnecro));//-N1 diagonal
			if (((int)i - (int)N1D) >= 0)
				coefficients.push_back(T(i - N1D, i, -params->DCnecro));//N1 diagonal
		}

		Snew.setFromTriplets(coefficients.begin(), coefficients.end());
		Snew.makeCompressed();
		RmS = Snew;

		solver.analyzePattern(RmS);
		if (solver.info() != Eigen::Success) {
			printf("NecrosisMap -> Analyze pattern failed! \n");
			return;
		}

		b.resize(numSpotsD);
		b.fill(0.);

		for (unsigned int i = 0; i < N1D; ++i)
			b(i) += params->DCnecro;
		for (unsigned int i = 0; i < numSpotsD; i = i + N1D)
			b(i) += params->DCnecro;
		for (unsigned int i = (N2D - 1)*N1D; i < numSpotsD; ++i)
			b(i) += params->DCnecro;
		for (unsigned int i = N1D - 1; i < numSpotsD; i = i + N1D)
			b(i) += params->DCnecro;

		//initializing Interpolation matrix
		coefficients.clear();
		SpMat InterpTmp(numSpots, numSpotsD);
		int cx, cy, whX, whY, rPx, rPy, whX2, whY2;
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

				coefficients.push_back(T(i, whX + whY*N1D, 2. / 3.));//center super node

				if (whX2 >= 0 && whX2 < (int)N1D && whY2 >= 0 && whY2 < (int)N2D)
					coefficients.push_back(T(i, whX2 + whY2*N1D, 1. / 3.));//other node

				

			}
			else {
				coefficients.push_back(T(i, whX + whY*N1D, 1.));//center super node
			}
		}

		InterpTmp.setFromTriplets(coefficients.begin(), coefficients.end());
		InterpTmp.makeCompressed();
		Interp = InterpTmp;

		bInterp.resize(numSpots);
		bInterp.fill(0.);

		for (unsigned int i = 0; i < N1; ++i)
			for (unsigned int j = 0; j < N2; ++j)
				if (i == 0 || i == N1 - 1 || j == 0 || j == N2 - 1) {
						bInterp(j*N1 + i) = 1./3.;
					}
	}

	state.resize(numSpots);
	state.fill(1.); //initial state

	changes.resize(numSpots);
	changes.fill(0.);
};

void NecrosisMap::addSink(unsigned int src) {
	changes(src) += 1.;
}

void NecrosisMap::removeSink(unsigned int src) {
	changes(src) -= 1;
}

void NecrosisMap::updateMap() {
	
	unsigned int iC;
	int aux, cx, cy, cz;
	for (unsigned int i = 0; i < numSpots; ++i)
		if (changes(i) == 1.) {
			//casting to super node
			aux = (int)(i % (N1*N2));
			cx = aux % (int)N1;
			cy = (int)((double)aux / (double)N1);
			cz = (int)((double)i / (double)(N1*N2));
			iC = (cx / 3) + (cy / 3)*N1D + (cz / 3)*(N1D*N2D);
			RmS.coeffRef(iC, iC) += 1.;
		}
		else if (changes(i) == -1) {
			aux = (int)(i % (N1*N2));
			cx = aux % (int)N1;
			cy = (int)((double)aux / (double)N1);
			cz = (int)((double)i / (double)(N1*N2));
			iC = (cx / 3) + (cy / 3)*N1D + (cz / 3)*(N1D*N2D);
			RmS.coeffRef(iC, iC) -= 1.;
		}

	RmS.makeCompressed();
	
	if (is3D) {
		solver3D.compute(RmS);
		if (solver3D.info() != Eigen::Success) {
			printf("NecrosisMap -> Compute step failed! \n");
			return;
		}
		Eigen::VectorXd x = solver3D.solve(b);
		if (solver3D.info() != Eigen::Success) {
			printf("Error: NecrosisMap -> Solving failed! \n");
			return;
		}
		//interpolating the state
		state = Interp*x + bInterp;
	}
	else {
		solver.factorize(RmS);
		if (solver.info() != Eigen::Success) {
			printf("NecrosisMap -> Factorization failed! \n");
			return;
		}

		Eigen::VectorXd x = solver.solve(b);
		if (solver.info() != Eigen::Success) {
			printf("Error: NecrosisMap -> Solving failed! \n");
			return;
		}
		//interpolating the state
		state = Interp*x + bInterp;
	}

	changes.fill(0.);
}


mxArray* NecrosisMap::getState() {

	mwSize dims[3] = { (mwSize)N1, (mwSize)N2, (mwSize)N3 };
	mxArray *outMap = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
	double *ptr = mxGetPr(outMap);
	memcpy(ptr, state.data(), numSpots*sizeof(double));

	return outMap;
}
