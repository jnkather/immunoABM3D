% 2016-2018, created by JN Kather and J Poleszczuk
% license: see separate license file

% What is this? this is an agent-based model of tumor cell - immune cells interactions
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
addpath('./subroutines_3D/'); % include generic subroutines for 3D modeling
addpath('./subroutines_plot/'); % include advanced subroutines for plotting

% all parameters for the model are stored in the structure "sysTempl".
% Hyperparameters are stored in the structure "cnst". If you want to 
% manually change parameters, you need to overwrite the respective value 
% in sysTempl.params or in cnst, for example by adding
% "sysTempl.params.TUps = 0.65" after the call to "getSystemParams"

% get system parameters. 2nd argument is domain size which must be
% multiplication of 3. [135 135 135] is a 2 mm cube domain.
[sysTempl, cnst] = getSystemParams('3D',[90 90 90]);  

cnst.nSteps   = 100; % how many iterations. 1 iteration = 12 hours
cnst.drawWhen = 20;  % update plot after ... iterations

% define vector of six free parameters and then copy this vector to the
% respetive fields in the parameter structure
pVec = [0.3933    0.3490    0.0067    0.0045    0.0120    0.0741]; % sample param vector
sysTempl.params.TUps = pVec(1);             % stem cell symmetric division
sysTempl.params.TUpmut = pVec(2);           % mutation probability
sysTempl.params.IMpkill =  pVec(3);         % killing probability
sysTempl.params.IMrateDynamic = pVec(4);    % lymphocyte influx
sysTempl.params.MPrateDynamic = pVec(5);    % macrophage influx
sysTempl.params.probSeedFibr =  pVec(6);    % fibrosis (stroma) generation

% add/override some global variables after loading the system
numExp = 1;                 % number of simulation runs, default 1
cnst.VideoFrameRep = 12;    % frame repetition if video is recorded, default 12
cnst.verbose = true;        % show simulation output on screen
saveImage = true;           % save simulation output image. requires verbose true
saveVideo = false;          % save simulation output video. requires verbose true
cnst.printImages = false;   % save high resolution simulation output image
cnst.lossFunction = 'default_stem'; % specify loss function for simulation

for i=1:numExp % for each simulation run
    
sysTempl.params.initialSeed = i; % reset random seed for reproducibility
expname = ['minimal_3D_',num2str(i)]; % experiment name for saving

[sysOut, lastFrame, summary, imWin, fcount, masterID] = ...
    runSystem(sysTempl,cnst,expname,saveImage,saveVideo); % run the system
end
