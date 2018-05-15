/*
*  Created by Jan Poleszczuk
*  Last modified November, 2017 by Jan Poleszczuk
*/

#pragma once

#include "stdafx.h"
#include "SIMparameters.h"

#include "suitesparse/cholmod.h"
//#include "suitesparse/cholmod_internal.h"


typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

class NecrosisMap
{

private:

	unsigned int N1, N2, N3;
	unsigned int numSpots;

	unsigned int N1D, N2D, N3D;
	unsigned int numSpotsD;

	bool is3D;

	SIMparameters* params;

	cholmod_sparse *A, *Interp;
	cholmod_dense *state, *b, *bInterp;
	cholmod_factor *L;
	cholmod_common c;

	SpMat RmS3D, Interp3D;
	Eigen::VectorXd b3D, bInterp3D, state3D;
	Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper> solver3D;

	char* changes;
	bool afterInit;

public:

	NecrosisMap();
	~NecrosisMap();

	void initialize(SIMparameters*);
	void updateMap();

	double getValue(unsigned int pos){ 
		if (is3D) {
			return state3D(pos);
		}
		else {
			return ((double*)state->x)[pos];
		}
	};

	void addSink(unsigned int);
	void removeSink(unsigned int);

	mxArray* getState();
};
