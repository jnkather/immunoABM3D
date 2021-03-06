%libeigen3-dev to be installed
%libsuitesparse
%libmetis

mex -v '-I/usr/include/eigen3/' ...
     SIMengine/Classes/stdafx.cpp SIMengine/Classes/mersenne.cpp ...
     SIMengine/Classes/SIMcore.cpp SIMengine/Classes/SIMparameters.cpp ...
     SIMengine/Classes/ChemotaxisMapSuiteS.cpp SIMengine/Classes/NecrosisMapSuiteS.cpp...
     SIMengine/Classes/Environment.cpp SIMengine/Classes/TumorCells.cpp ...
     SIMengine/Classes/Lymphocytes.cpp SIMengine/Classes/Macrophages.cpp ...
     SIMengine_interface_mex.cpp ...
 -lamd -lcamd -lccolamd -lcholmod -lcolamd -lmetis -lsuitesparseconfig -lcxsparse -lblas -llapack

 copyfile('stdafx.mexa64','SIMengine_interface_mex.mexa64')