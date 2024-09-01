clear all
clc

% LMC with Single Loop Gaussian Smoothing
% This is the beginning of trying out amplitude polynomial function
% rng(1);

%% 

% Runs
numRuns = 5;

% Setup Chirp parameters
fs = 1000;
Td = [0.02, 0.075, 0.1, 0.5, 1];

phi{1,1} = [0, 100, 50].'; 
phi{1,2} = [0, 130, 40].'; 
rho{1,1} = 1; %flipud([-56.7741802753438, 170.236120002308, -175.321926949479, 67.0086426225657, -6.58823998781464, 1.43827467211612, 0.0733688653715676].'); % Gaussian window (0.9 variance)
rho{1,2} = 1; %flipud([-56.7741802753438, 170.236120002308,	-175.321926949479,	67.0086426225657, -6.58823998781464, 1.43827467211612, 0.0733688653715676].');
snr = 60;
tol = 1e-8;

% Visualization
bDisplayPlots = true;

% Tuning for LMC
numParticles  = 20;
stepSizePhi   = [1e-5, 1e-5, 5e-5, 5e-5].'; % We do not care about phi0. 
stepSizeConst = 0.1; 
stepSizeMax   = 2;
stepSizeMin   = 5e-6;
stepNoiseVar  = 1e-5;
avgConst      = 0.8; 
tempConst     = 1;
numIterLmc    = [500, 500, 500, 600, 1000];   % This is T. Having more than one
noiseVarInit  = [0.1, 0.1, 1, 1, 1];
noiseVarMin   = 0.001;
numIterSmooth = [5, 5, 10, 10, 50];
bGaussSmooth  = [true, true, true, true, true];
bEnableLangevin = true;

gamma = [0.8, 0.5, 0, 0, 0];
initValMinMax = [0, 200];

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
lmcTuning.tempConst     = tempConst;
lmcTuning.stepSizeConst = stepSizeConst;
lmcTuning.stepSizeMax   = stepSizeMax;
lmcTuning.stepSizeMin   = stepSizeMin;
lmcTuning.stepNoiseVar  = stepNoiseVar;
lmcTuning.avgConst      = avgConst;
lmcTuning.numParticles  = numParticles;
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

for runInd = 1:numRuns

    disp('------------------------------');
    disp(['Run Number ', num2str(runInd)]);
    disp('------------------------------');

    % Some more settings which are run dependent
    cpeSetting.Td     = Td(1,runInd);
    cpeSetting.gamma  = gamma(1,runInd);
    lmcTuning.noiseVarInit  = noiseVarInit(1,runInd);
    lmcTuning.numIterLmc    = numIterLmc(1,runInd);
    lmcTuning.numIterSmooth = numIterSmooth(1,runInd);
    lmcTuning.bGaussSmooth  = bGaussSmooth(1,runInd);

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
    lmcTuning.initParams = repmat(lmc.param(:,bestParticleInd), 1, numParticles);
    % lmcTuning.stepSize   = lmc.stepSize(:,bestParticleInd);

    if runInd ~= numRuns
        delete(lmc)
    end
end

%% Generate Plots

% Setup class for plottingf
plt = classGenPlots_SL(lmc);

% Call function to produce plots
plt.genAllPlots();
