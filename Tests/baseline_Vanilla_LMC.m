clear all
clc

% LMC with Single Loop Gaussian Smoothing
% This is the beginning of trying out amplitude polynomial function
% rng(1);
% Vanilla Langevin Monte Carlo

%% 

% Runs
numOptRuns = 2;

% Setup Chirp parameters
fs = 500;
Td = [0.25, 0.75];

phi{1,1} = [0, 10, 60].'; 
phi{1,2} = [0, 20, 85].';  
rho{1,1} = 1;
rho{1,2} = 1;
snr = [12];
tol = 1e-8;

% Visualization
bDisplayPlots = false;

% Tuning for LMC
numParticles  = [50, 20];
stepSizePhi{1,1} = 0.01  * [1, 1, 1, 1].'; % We do not care about phi0.
stepSizePhi{1,2} = 0.001 * [1, 1, 1, 1].'; % We do not care about phi0.
stepSizeConst = 0.1; 
stepSizeMax   = 0.5;
stepSizeMin   = 5e-6;
stepNoiseVar  = [5e-6, 5e-6];
avgConst      = [1, 1]; 
tempConst     = [1, 1];
numIterLmc    = [100, 100];   % This is T. Having more than one
noiseVarInit  = [1, 1];
noiseVarMin   = 0.001;
numIterSmooth = [10, 30];
bGaussSmooth  = [true, true];
bEnableLangevin = true; 
bMetropolisOn   = [true, true];
bApplyWin       = [false, false];
gamma           = [0.5, 0];
initValMinMax   = [0, 100;
                   0, 100];
perturbVar = 0;

% How many params?
Nc = length(phi); % Number of chirps
numParams = 0;
for c = 1:Nc
    numParams = numParams + numel(phi{1,c}(2:end));
end

%% Statistics Tests

numSnrRuns  = numel(snr);
numStatRuns = 1;

cellSnrRunStats = cell(1, numSnrRuns);

for snrRunInd = 1:numSnrRuns
    cellSnrRunStats{1, snrRunInd} = zeros(numParams, numStatRuns);
end

%% Build setup & tuning structures

cpeSetting.fs        = fs;
cpeSetting.Nc        = Nc;  
cpeSetting.phi       = phi;
cpeSetting.rho       = rho;
cpeSetting.numParams = numParams;
cpeSetting.minObjTol = tol;

lmcTuning.stepSizeConst = stepSizeConst;
lmcTuning.stepSizeMax   = stepSizeMax;
lmcTuning.stepSizeMin   = stepSizeMin;
lmcTuning.noiseVarMin   = noiseVarMin;
lmcTuning.bDisplayPlots = bDisplayPlots;
lmcTuning.bEnableLangevin = bEnableLangevin;
lmcTuning.initValMinMax = initValMinMax;
lmcTuning.initParams    = [];

%% Simulation
% 
% if numRuns ~= numel(Td) || numRuns ~= numel(numIterLmc)
%     error('Number of runs does not match setup for Signal Length and LMC Iterations');
% end

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('Starting simulation ... ');
fprintf('\n');
tic
for snrRunInd = 1:numSnrRuns
    tic;

    disp('--------------------------------------------------------');
    disp(['Stats Run Number ', num2str(snrRunInd), ' of ', num2str(numSnrRuns)]);
    disp('--------------------------------------------------------');

    % Choose snr for test
    cpeSetting.snr = snr(1,snrRunInd);

    % Run loops for statistics (VoE)
    for statRunInd = 1:numStatRuns

        disp('--------------------------------------------------------');
        disp(['Stats Run Number ', num2str(statRunInd), ' of ', num2str(numStatRuns)]);
        disp('--------------------------------------------------------');

        % Run optimization loops
        for optRunInd = 1:numOptRuns
        
            disp('--------------------------------------------------------');
            disp(['Optimization Run Number ', num2str(optRunInd), ' of ', num2str(numOptRuns)]);
            disp('--------------------------------------------------------');
        
            % Some more settings which are run dependent
            cpeSetting.Td           = Td(1,optRunInd);
            cpeSetting.gamma        = gamma(1,optRunInd);
            cpeSetting.bApplyWin    = bApplyWin(1,optRunInd);
            lmcTuning.noiseVarInit  = noiseVarInit(1,optRunInd);
            lmcTuning.numIterLmc    = numIterLmc(1,optRunInd);
            lmcTuning.numIterSmooth = numIterSmooth(1,optRunInd);
            lmcTuning.bGaussSmooth  = bGaussSmooth(1,optRunInd);
            lmcTuning.tempConst     = tempConst(1,optRunInd);
            lmcTuning.numParticles  = numParticles(1, optRunInd);
            lmcTuning.stepNoiseVar  = stepNoiseVar(1, optRunInd);
            lmcTuning.bMetropolisOn = bMetropolisOn(1, optRunInd);
            lmcTuning.stepSize      = stepSizePhi{1,optRunInd}; %repmat(stepSizePhi{1,runInd}, Nc, 1);
            lmcTuning.avgConst      = avgConst(1,optRunInd);
        
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
            bestParticleInd = lmc.bestParticleInd ;
            bestParamsSoFar = lmc.param(:,bestParticleInd);
            noisyMat        = randn(numParams, numParticles(1,optRunInd));
            lmcTuning.initParams = repmat(bestParamsSoFar, 1, numParticles(1, optRunInd)) + perturbVar * noisyMat;
        
            if optRunInd ~= numOptRuns
                delete(lmc)
            end
        end % opt run end

        cellSnrRunStats{1, snrRunInd}(:, statRunInd) = lmc.sqrError; 

        if statRunInd ~= numStatRuns
            delete(lmc)
        end

        lmcTuning.initParams = [];
    end % stat run end
    
end % snr run end
toc;

%% Calculate Variance of Error 

VoE = zeros(numParams, numSnrRuns);
for snrRunInd = 1:numSnrRuns
    VoE(:, snrRunInd) = mean(cellSnrRunStats{1, snrRunInd}(:, snrRunInd), 2);

    disp(['Variance of Error at ', num2str(snr(snrRunInd)), ' dB = ' , num2str(log10(VoE(:, snrRunInd).'))]);
end

%% Generate Plots

% Setup class for plottingf
plt = classGenPlots_SL(lmc);

% Call function to produce plots
plt.genAllPlots();
