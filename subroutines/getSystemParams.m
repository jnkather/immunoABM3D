% JN Kather, NCT Heidelberg, 2017
%
% IMinfluxRate will only be used if IMrateDynamic is not supplied, same for
% MP

function [mySystem, cnst] = getSystemParams(dims,spaceSetting)
disp('requested system parameters');

% START GENERAL SYSTEM PROPERTIES -------------------------------------------
mySystem.params.initialSeed = 1;        % initial random seed, default 1
mySystem.params.seedUnderneath = false; % if fibrosis can be seeded below cell
mySystem.params.CellDieAtProl = 0; % cell dies upon proliferation attempt, default 0 (increase by chemo)
% END SYSTEM PROPERTIES -------------------------------------------

% START INITIALIZE TUMOR CELLS -------------------------------------------
mySystem.params.TUpprol = 0.5055;   % HISTO GENERATED - probability of proliferation
mySystem.params.TUpmig = 0.35;      % probability of migrating, default 0.35
mySystem.params.TUpdeath = 1-(1-0.0319)^4;  % HISTO GENERATED - probability of spontaneous death
mySystem.params.TUpmax = 10;        % FIXED max. proliferation capacity, default 10
mySystem.params.TUdanti = 0.1;      % FIXED antigenicity strength of mutating tumor cell, default 1(arbitrary)
mySystem.params.TUdadju = 1;        % FIXED adjuvanticity strength of dying tumor cell, default 1 (arbitrary)
mySystem.params.TUdamageThresh = 2; % T cell inflicted damage threshold for tumor cell, default 2
if strcmp(dims,'2D') % parameters that depend on the dimensionality
    mySystem.params.TUps = 0.42;    % 2D probability of symmetric division, default 0.4
	mySystem.params.TUpmut = 0.075;     % mutation probability (increases antigenicity)
	mySystem.params.bounds.TUpmut = [0.02 0.2];
else % 3D 
    mySystem.params.TUps = 0.42;    % 3D probability of symmetric division, default 0.42
    mySystem.params.bounds.TUps = [0.35 0.50]; % in 3D, this is a viable free parameter for model fitting
	mySystem.params.TUpmut = 0.075;     % mutation probability (increases antigenicity)
	mySystem.params.bounds.TUpmut = [0.010 0.400];
end
% END INITIALIZE TUMOR CELLS ---------------------------------------------

% START INITIALIZE LYMPHOCYTES ------------------------------------------
mySystem.params.IMkmax = 5;          % FIXED killing capacity of immune cells, default 5
mySystem.params.IMpmax = 5;          % FIXED proliferation capacity of immune cells, default 10
mySystem.params.IMpmig = 0.8;        % probability of lymphocyte migration, default 0.8
mySystem.params.IMpkill = 0.02;      % probability of killing
mySystem.params.bounds.IMpkill = [0.005 0.1];
mySystem.params.IMrwalk = 0.5;       % random influence on movement, default 0.75
mySystem.params.IMspeed = 97;        % speed of immune cell movement, default 97
mySystem.params.IMpprol = 0.0449/mySystem.params.IMspeed;   % HISTO GENERATED - probability of proliferation
mySystem.params.IMpdeath = (1-(1-0.0037)^4)/mySystem.params.IMspeed;  % HISTO GENERATED - probability of spontaneous death
mySystem.params.engagementDuration = 48; % how many intermed. steps is a killing cell engaged? default 48 (=6 hrs)
mySystem.params.adjuThresh = 0;      % FIXED adjuvanticity threshold for lymphocyte activation, default 0 (arbitrary)
mySystem.params.antiThresh = 0.3;    % FIXED antigenicity threshold for lymphocyte activation, default 0.3 (arbitrary)
mySystem.params.adjuDecay = 0.2;     % adjuvanticity decay in each iteration
mySystem.params.IMinfluxProb = 0.33; % probability of immune cell influx, default 0.33
if strcmp(dims,'2D') % parameters that depend on the dimensionality
    mySystem.params.IMinfluxRate = 1;      % 2D how many lymphocytes appear simultaneously, always fixed at 1
    mySystem.params.IMrateDynamic = 0.0075;   % 2D how does lymphocyte influx scale with increasing tumor size
    mySystem.params.bounds.IMrateDynamic = [0.0002 0.02];
else % 3D
    mySystem.params.IMinfluxRate = 1; 	   % 3D how many lymphocytes appear simultaneously, always fixed at 1
    mySystem.params.IMrateDynamic = 0.0075;   % 3D how does lymphocyte influx scale with increasing tumor size
    mySystem.params.bounds.IMrateDynamic = [0.0002 0.0250];
end
% END INITIALIZE LYMPHOCYTES --------------------------------------------

% START INITIALIZE MACROPHAGES ------------------------------------------
mySystem.params.MPpmax = 5;        % proliferation capacity of macrophages, default 10
mySystem.params.MPpprol = mySystem.params.IMpprol;     % macrophage probability of proliferation
mySystem.params.MPpmig = mySystem.params.IMpmig;       % macrophage probability of migrating
mySystem.params.MPpdeath = mySystem.params.IMpdeath;   % macrophage probability of spontaneous death
mySystem.params.MPrwalk = 0.5;       	% FIXED random walking of macrophages 0...1, default 0.5
mySystem.params.MPspeed = mySystem.params.IMspeed;       % speed of macrophage movement, default 97
mySystem.params.MPppola = 0.1;     		% FIXED polarization: transition probability from inactive to pro-tumor macrophages
mySystem.params.MPprepola = 0.16666;	% FIXED repol. based on M2/(M1+M2) ratio from PMC4839334, N=205 pts. 0.4 = 1-(0.1/0.16666) 
mySystem.params.MPeffect = -1;    		% FIXED effect strength of tumor promoting macrophages, don't forget the minus
mySystem.params.adjuRange = 7;      	% effect range of macrophages on AdjuMap (radius of circle), default 7
mySystem.params.MPinfluxProb = 0.33; 	% probability of macrophage influx, fixed at 0.33
if strcmp(dims,'2D') % parameters that depend on the dimensionality
    mySystem.params.MPinfluxRate = 1;    % 2D how many macrophages appear simultaneously, always fixed at 1
    mySystem.params.MPrateDynamic = 0.0075;   % 2D how does macrophage influx scale with increasing tumor size
    mySystem.params.bounds.MPrateDynamic = [0.0020 0.1000];
else % 3D
    mySystem.params.MPinfluxRate = 1;    % 3D how many macrophages appear simultaneously, always fixed at 1
    mySystem.params.MPrateDynamic = 0.0075;   % 3D how does macrophage influx scale with increasing tumor size
    mySystem.params.bounds.MPrateDynamic = [0.0020 0.1000];
end
% END INITIALIZE MACROPHAGES --------------------------------------------

% START INITIALIZE CHEMOTAXIS MAP ------------------------------------------
mySystem.params.DCchemo = 100; %diffusion/consumption in the stationary diffusion-consumption equation
mySystem.params.SCchemo = 1; %secretion/consumption in the stationary diffusion-consumption equation
% END INITIALIZE CHEMOTAXIS MAP  ---------------------------------

% START INITIALIZE NECROSIS
oxygenDiffusion = 2.5*1e-5; %ref (Powathil, et al., Comput Math Met Med, 2012), in cm^2/s
oxygenPointConsumption = 3.8*1e-13; %ref (Powathil, et al., Comput Math Met Med, 2012), in cm^2*O2/(s*cell)
carryingCapacity = 2.1*10^11/2; %ref (Powathil, et al., Comput Math Met Med, 2012); scaled by 2, cells
dx = 14.9*1e-4;%grid spot width, in cm
mySystem.params.DCnecro = oxygenDiffusion/oxygenPointConsumption/carryingCapacity/dx^2; %diffucion/(TU cell consumption) in the stationary diffusion-consumption equation + numerical scheme correction
mySystem.params.TCnecro = 1; % FIXED total consumption by a single TU cell (because of normalization D/c)
mySystem.params.necThresh = 0.04; % threshold of nutrients value below which the TU cell might die
% END INITIALIZE NECROSIS

% START INITIALIZE FIBROSIS  ---------------------------------
mySystem.params.smoothRadius = 3; 		% for smoothing fibrotic maps, default 3
mySystem.params.probSeedFibr = 0.06;   	% probability of turning into fibrosis, model fitting: 0.06
mySystem.params.bounds.probSeedFibr = [0.002 0.2];   % probability of turning into fibrosis, model fitting: 0.10
mySystem.params.fibrFrac = 0.3;    		% FIXED size of fibrotic seed, 0...1, default 0.3
mySystem.params.stromaPerm = 0.0005;	% 0 = stroma not permeable, 1 = fully permeable, default very small (0.0005)
% END INITIALIZE FIBROSIS  ---------------------------------

% START DEFINING ADDITIONAL CONSTANTS -----------------------------------
cnst.verbose = true;            % draw intermediary steps? default true
cnst.createNewSystem = true;    % create new system at the start, default true
cnst.maxCells = Inf; %maximal number of tumor cells
cnst.averageOut = true;

% cnst.saveImage = true;          % save image directly after each iteration, default true
% cnst.doImage = false;           % plot result again afterwards, default false
% cnst.doVideo = false;           % create a video afterwards, default false
% cnst.doSummary = true;          % summarize the result, default true
cnst.inTumor = 1;               % defines "in tumor" ROI, default 1
cnst.marginSize = round(67/2);  % default "invasive margin" ROI, default 67/2 = 0.5 mm
cnst.around = round(67*2);      % defines "adjacent tissue" ROI, default 67*2 = 2 mm
cnst.requireAlive = 150;        % require tumor to be alive for some time
cnst.maxAntigenicity = 1;       % maximum antigenicity for tumor cells
cnst.tumorColorLevels = 100;    % how many tumor color shades (must be multiple of 4!)
cnst.antigenKernelWidth = 0.05; % width of ksdensity kernel for antigenicity plot
cnst.defaultAntigenicity = 0.05; % antigenicity of first tumor cell
cnst.smoRegion = strel('disk',5,0); % region smoother for topography statistics
cnst.topoColorNorm = 0.03;      % max value for topography maps (red/blu)
cnst.VideoFrameRep = 5;         % repetition of video frames
cnst.lossFunction = 'default';  % define the loss function for system summary, 'default' or 'none'
cnst.penalty = 50;              % define penalty for immune cell win for loss function, default 50
% END DEFINING ADDITIONAL CONSTANTS -----------------------------------

% START DEFINING DIMENSION VARIABLES  -----------------------------------
if strcmp(dims,'2D') % parameters that are specific for 2D
mySystem.is3D = false;
mySystem.grid.N = spaceSetting(1);  % domain dimension vertical, default 380
mySystem.grid.M = spaceSetting(2);  % domain dimension horizontal, default 380
mySystem.params.mycellsize = 6;     % circle size in scatter plot, default 6
elseif strcmp(dims,'3D') % parameters that are specific for 3D
mySystem.is3D = true;
mySystem.grid.N = spaceSetting(1);  % domain dimension 1, default 90
mySystem.grid.M = spaceSetting(2);  % domain dimension 2, default 90
mySystem.grid.P = spaceSetting(3);  % domain dimension 3, default 90
else
    error('dimension must be 2D or 3D')
end

% START PLAUSIBILITY CHECK ----------------------------------------------
checkPlausibility(mySystem);
% END PLAUSIBILITY CHECK ------------------------------------------------

end
