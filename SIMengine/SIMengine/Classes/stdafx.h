/*
*  Created by Jan Poleszczuk
*  Last modified January, 2018 by Jan Poleszczuk
*/


// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#include <stdio.h>
//#include <tchar.h>


#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <functional>
#include "mex.h"

#include "randomc.h"
#include "time.h"

#include "Eigen/Core"
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"

#define _CRTDBG_MAP_ALLOC  

#include <stdlib.h>  

#ifdef _MSC_VER
	#include <crtdbg.h>
#endif

#ifdef _DEBUG
#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
// Replace _NORMAL_BLOCK with _CLIENT_BLOCK if you want the
// allocations to be of _CLIENT_BLOCK type
#else
#define DBG_NEW new
#endif

#ifndef TRUE
#define TRUE true
#endif



// TODO: reference additional headers your program requires here
