/*
*  Created by Jan Poleszczuk
*  Last modified November, 2017 by Jan Poleszczuk
*/

#pragma once

#include "stdafx.h"
#include "SIMparameters.h"
#include "Eigen/Sparse"
#include "Eigen/SparseCholesky"
#include "Eigen/IterativeLinearSolvers"

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

	SpMat RmS, Interp;
	Eigen::VectorXd b, bInterp, state;
	Eigen::SimplicialLDLT<SpMat> solver;
	Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper> solver3D;
	//Eigen::BiCGSTAB<SpMat> solver3D;

	Eigen::VectorXd changes;


public:

	NecrosisMap();
	~NecrosisMap();

	void initialize(SIMparameters*);
	void updateMap();

	double getValue(unsigned int pos){ return state(pos); };

	void addSink(unsigned int);
	void removeSink(unsigned int);

	mxArray* getState();
};
