clear all
clc

% LMC with Annealed Iterative Gaussian Smoothing
% This is the beginning of trying out amplitude polynomial function
rng(1);

%% 

% Setup Chirp parameters
fs = 300;
Td = 1;

phi{1,1} = [0, 60, 40, 10].'; 
phi{1,2} = [0, 40, 35, 12].'; 
rho{1,1} = [-1.83106519981518, 4.17637124212063, -3.32116203820533, 0.971682696168882].';
rho{1,2} = [-0.523573178531541, 1.57521448968648, -1.91623270407228, 0.996272121844967].';
snr = 20;
tol = 1e-8;

% Visualization
bDisplayPlots = true;

% Tuning for LMC
numParticles  = 10;
stepSizePhi   = 0.01 * [1, 1, 1, 1].'; % We do not care about phi0. 
stepSizeConst = 0.01; % Turn ON/OFF Gradient based stepsize normalization
stepSizeMax   = 1;
stepNoiseVar  = 0.001;
avgConst      = 0.9;  % Coefficient for exponential smoothing of Gradient and trace of Hessian
tempConst     = 1;
numIterLmc    = 300;  % This is T
numIterNoise  = 10;   % This is L
noiseVarInit  = 1;
noiseVarFinal = 0.00001;
bGaussSmooth  = true;
bEnableLangevin = true;
numIterSmooth = 50;

%% Build setup & tuning structures

% How many params?
Nc = length(phi); % Number of chirps
numParams = 0;
for c = 1:Nc
    numParams = numParams + numel(phi{1,c}(2:end));
end
stepSize = repmat(stepSizePhi, Nc, 1);

cpeSetting.fs        = fs;
cpeSetting.Td        = Td;
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
lmcTuning.stepNoiseVar  = stepNoiseVar;
lmcTuning.avgConst      = avgConst;
lmcTuning.numParticles  = numParticles;
lmcTuning.numIterLmc    = numIterLmc;
lmcTuning.noiseVarInit  = noiseVarInit;
lmcTuning.noiseVarFinal = noiseVarFinal;
lmcTuning.numIterNoise  = numIterNoise;
lmcTuning.bDisplayPlots = bDisplayPlots;
lmcTuning.numIterSmooth = numIterSmooth;
lmcTuning.bGaussSmooth  = bGaussSmooth;
lmcTuning.bEnableLangevin = bEnableLangevin;

%% Simulation

% Call the constructor to setup
lmc = classLangevinMonteCarlo(lmcTuning, cpeSetting);

% Start the simulation
tic;
lmc = lmc.runLmcCore();
toc;

%% Generate Plots

% Setup class for plottingf
plt = classGenPlots(lmc);

% Call function to produce plots
plt.genAllPlots();

