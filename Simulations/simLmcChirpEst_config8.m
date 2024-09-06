clear all
clc

% LMC with Single Loop Gaussian Smoothing
% This is the beginning of trying out amplitude polynomial function
% rng(1);

%% 

% Runs
numRuns = 3;

% Setup Chirp parameters
fs = 1000;
Td = [0.05, 0.06, 0.7];

phi{1,1} = [0, 60,  90, 20].'; 
phi{1,2} = [0, 90, 100, 15].';  
rho{1,1} = 1;%flipud([-1.79867186746021	4.12672434196200	-3.29928369720428	0.969161440768821].'); % Gaussian window (0.9 variance)
rho{1,2} = 1;%flipud([-1.79867186746021	4.12672434196200	-3.29928369720428	0.969161440768821].');
snr = 15;
tol = 1e-8;

% Visualization
bDisplayPlots = true;

% Tuning for LMC
numParticles  = [50, 10, 5];
stepSizePhi   = [1e-2, 1e-2, 1e-2, 1e-2].'; % We do not care about phi0. 
stepSizeConst = 0.1; 
stepSizeMax   = 1;
stepSizeMin   = 5e-6;
stepNoiseVar  = [5e-6, 7e-5, 7e-5];
avgConst      = 0.95; 
tempConst     = [1, 1, 1];
numIterLmc    = [200, 200, 500];   % This is T. Having more than one
noiseVarInit  = [0.5, 0.5, 1];
noiseVarMin   = 0.001;
numIterSmooth = [20, 50, 50];
bGaussSmooth  = [true, true, true];
bEnableLangevin = true;

gamma = [0.1, 0.05, 0];
initValMinMax = [0, 100];

%% Build setup & tuning structures

% How many params?
Nc = length(phi); % Number of chirps
numParams = 0;
for c = 1:Nc
    numParams = numParams + numel(phi{1,c}(2:end));
end
stepSize = repmat(stepSizePhi, Nc, 1);

cpeSetting.fs        = fs;
cpeSetting.Nc        = Nc;  
cpeSetting.phi       = phi;
cpeSetting.rho       = rho;
cpeSetting.snr       = snr;
cpeSetting.numParams = numParams;
cpeSetting.minObjTol = tol;

lmcTuning.stepSize      = stepSize;
lmcTuning.stepSizeConst = stepSizeConst;
lmcTuning.stepSizeMax   = stepSizeMax;
lmcTuning.stepSizeMin   = stepSizeMin;
lmcTuning.avgConst      = avgConst;
lmcTuning.noiseVarMin   = noiseVarMin;
lmcTuning.bDisplayPlots = bDisplayPlots;
lmcTuning.bEnableLangevin = bEnableLangevin;
lmcTuning.initValMinMax = initValMinMax;
lmcTuning.initParams    = [];

%% Simulation

if numRuns ~= numel(Td) || numRuns ~= numel(numIterLmc)
    error('Number of runs does not match setup for Signal Length and LMC Iterations');
end

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('Starting simulation ... ');
fprintf('\n');
tic;
for runInd = 1:numRuns

    disp('--------------------------------------------------------');
    disp(['Run Number ', num2str(runInd), ' of ', num2str(numRuns)]);
    disp('--------------------------------------------------------');

    % Some more settings which are run dependent
    cpeSetting.Td    = Td(1,runInd);
    cpeSetting.gamma = gamma(1,runInd);
    lmcTuning.noiseVarInit  = noiseVarInit(1,runInd);
    lmcTuning.numIterLmc    = numIterLmc(1,runInd);
    lmcTuning.numIterSmooth = numIterSmooth(1,runInd);
    lmcTuning.bGaussSmooth  = bGaussSmooth(1,runInd);
    lmcTuning.tempConst     = tempConst(1,runInd);
    lmcTuning.numParticles  = numParticles(1, runInd);
    lmcTuning.stepNoiseVar  = stepNoiseVar(1,runInd);

    if runInd == numRuns
        % bEnableLangevin = false;
    end

    % Call the constructor to setup
    lmc = classLangevinMonteCarlo_SL(lmcTuning, cpeSetting);

    % Start the simulation
    tic;
    lmc = lmc.runLmcCore_SL();

    if lmc.bStopSim
        return;
    end
    toc;

    % Init next run using params from previous run
    bestParticleInd      = lmc.bestParticleInd ;
    lmcTuning.initParams = repmat(lmc.param(:,bestParticleInd), 1, numParticles(1, runInd));
    % lmcTuning.stepSize   = lmc.stepSize(:,bestParticleInd);

    if runInd ~= numRuns
        delete(lmc)
    end
end
toc;

%% Generate Plots

% Setup class for plottingf
plt = classGenPlots_SL(lmc);

% Call function to produce plots
plt.genAllPlots();
