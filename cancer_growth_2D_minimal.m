% 2016-2018, created by JN Kather and J Poleszczuk
% license: see separate license file

% What is this? this is an agent-based model of tumor cell - immune cells 
% interactions
% How does it work? open this file in Matlab and run it. 
% 
% the code relies on a simulation engine that it written in C++. It has
% been compiled on a Windows computer. Before you run the model on a 
% Linux or MacOS computer, you might have to re-compile it. 
% See instructions in ./SIMengine
% 
% the C++ code relies on the open source software SuiteSparse (DLLs are
% included in the ./SIMengine folder. For license and disclaimer please see
% https://github.com/jlblancoc/suitesparse-metis-for-windows/

close all; clear variables; format compact; clc % avoid spillover 
addpath('./SIMengine/'); % include SIMengine (MEX-based simulation engine)
addpath('./subroutines/'); % include generic subroutines for 2D and 3D
addpath('./subroutines_2D/'); % include generic subroutines for 2D modeling
addpath('./subroutines_plot/'); % include advanced subroutines for plotting

% all parameters for the model are stored in the structure "sysTempl".
% Hyperparameters are stored in the structure "cnst". If you want to 
% manually change parameters, you need to overwrite the respective value 
% in sysTempl.params or in cnst, for example by adding
% "sysTempl.params.TUps = 0.65" after the call to "getSystemParams"

% get system parameters. 2nd argument is domain size which must be
% multiplication of 3. [270 270] is a 4 mm square domain. minimum domain
% size is 150x150 (because otherwise the scale bar does not fit in).
%
[sysTempl, cnst] = getSystemParams('2D',[300 300]);
sysTempl.params.mycellsize = 4; % cell size for plotting

% override some system parameters (just for this example)
numExp = 1;                 % number of simulation runs, default 1
cnst.verbose = true;        % show simulation output on screen
saveImage = false;           % save simulation output image. requires verbose true
saveVideo = false;          % save simulation output video. requires verbose true
sysTempl.params.IMrateDynamic = sysTempl.params.IMrateDynamic*0.4; % override for this example
sysTempl.params.MPrateDynamic = sysTempl.params.MPrateDynamic*0.4; % override for this example
cnst.nSteps   = 400;        % iterations (1 iteration = 12 hours)
cnst.drawWhen = 10;         % draw after ... iterations
cnst.printImages = false;   % print high resolution image

for i=1:numExp
sysTempl.params.initialSeed = i;

expname = ['minimal_2D_',num2str(i)];    % experiment name for saving

[sysOut, lastFrame, summary, imWin, fcount, masterID] = ...
    runSystem(sysTempl,cnst,expname,saveImage,saveVideo); % run the system
end
