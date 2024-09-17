clear all
clc

% NA-LMC Tests

%% 
exptName = 'NA_LMC_Test';
snr = 3;

% Runs #
numSnrRuns  = numel(snr);
numStatRuns = 1;
numOptRuns  = 2;

% Setup Chirp parameter
fs = 1000;
Td = [0.5, 1];

phi{1,1} = [0, 10, 40, -70, 110].'; 
phi{1,2} = [0, 50, 60, -90, 105].';  
rho{1,1} = flipud([-1.79867186746021, 4.12672434196200, -3.29928369720428, 0.969161440768821].');
rho{1,2} = flipud([-1.79867186746021, 4.12672434196200, -3.29928369720428, 0.969161440768821].');
tol = 1e-8;

% Visualization
bDisplayPlots = true;

% Tuning for LMC
numParticles  = [50, 10]; % Number of sample runs in parallel
stepSizePhi{1,1} = 0.01      * [1, 1, 1, 1, 1, 1, 1, 1].';
stepSizePhi{1,2} = 0.0000001 * [1, 1, 1, 1, 1, 1, 1, 1].'; 
stepSizeConst = 0.1; 
stepSizeMax   = 0.5;
stepSizeMin   = 5e-6;
stepNoiseVar  = [5e-6, 5e-6];
avgConst      = [1, 1]; 
tempConst     = [1, 2];
numIterNoise  = 10;
numIterLmc    = [100, 1000]; 
noiseVarInit  = [1, 1];
noiseVarFinal   = 0.001;
numIterSmooth = [1, 1];
bGaussSmooth  = [false, false];
bEnableLangevin = true; 
bMetropolisOn   = [true, true];
bApplyWin       = [false, false];
gamma           = [0.6, 0];
initValMinMax   = [-500, 500;
                   -500, 500;
                   -500, 500;
                   -500, 500];
perturbVar = 0;

% How many params?
Nc = length(phi); % Number of chirps
numParams    = 0;
numAmpParams = 0;
for c = 1:Nc
    numParams    = numParams    + numel(phi{1,c}(2:end));
    numAmpParams = numAmpParams + numel(rho{1,c});
end

%% Statistics Tests

cellSnrRunStats = cell(2, numSnrRuns);

for snrRunInd = 1:numSnrRuns
    cellSnrRunStats{1, snrRunInd} = zeros(numParams, numStatRuns);
    cellSnrRunStats{2, snrRunInd} = zeros(numAmpParams, numStatRuns);
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
lmcTuning.noiseVarFinal   = noiseVarFinal;
lmcTuning.bDisplayPlots = bDisplayPlots;
lmcTuning.bEnableLangevin = bEnableLangevin;
lmcTuning.numIterNoise  = numIterNoise;
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
% bStopSim = false;

tic;
for snrRunInd = 1:numSnrRuns

    disp('--------------------------------------------------------');
    disp(['SNR Run Number ', num2str(snrRunInd), ' of ', num2str(numSnrRuns)]);
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
            lmc = classLangevinMonteCarlo(lmcTuning, cpeSetting);
        
            % Start the simulation
            lmc = lmc.runLmcCore();
        
            if lmc.bStopSim
                bStopSim = true;
                return;
            else
                bStopSim = false;
            end
        
            % Init next run using params from previous run
            bestParticleInd = lmc.bestParticleInd ;
            bestParamsSoFar = lmc.param(:,bestParticleInd);
            noisyMat        = randn(numParams, numParticles(1,optRunInd));
            lmcTuning.initParams = repmat(bestParamsSoFar, 1, numParticles(1, optRunInd)) + perturbVar * noisyMat;
        
            if bStopSim
                return;
            end
            if optRunInd ~= numOptRuns
                delete(lmc)
            end

        end % opt run END

        if bStopSim
            return;
        end

        % Save errors
        cellSnrRunStats{1, snrRunInd}(:, statRunInd) = lmc.cpe{1, lmc.bestParticleInd}.sqrPhiError; 
        cellSnrRunStats{2, snrRunInd}(:, statRunInd) = lmc.cpe{1, lmc.bestParticleInd}.sqrRhoError; 

        % Save data to MAT files
        saveSimData.classData = lmc;
        saveSimData.runStats  = cellSnrRunStats;      
        saveFileFormat = '.mat';
        saveFilePath   = 'Data/Output/ForPaper/MATFiles';
        saveFileName   = [fullfile(saveFilePath,exptName), '_', num2str(snr(snrRunInd)), '_dB_StRun_', num2str(statRunInd), saveFileFormat];
        save(saveFileName, "saveSimData", '-v7.3');

        if statRunInd ~= numStatRuns
            delete(lmc)
        end
        lmcTuning.initParams = [];

    end % stat run END    

end % snr run END
toc;

%% Calculate Variance of Error 

VoE_phi = zeros(numParams,    numSnrRuns);
VoE_rho = zeros(numAmpParams, numSnrRuns);
for snrRunInd = 1:numSnrRuns
    VoE_phi(:, snrRunInd) = mean(cellSnrRunStats{1, snrRunInd}, 2);
    VoE_rho(:, snrRunInd) = mean(cellSnrRunStats{2, snrRunInd}, 2);

    disp(['Phase VoE (log10) at ', num2str(snr(snrRunInd)), ' dB SNR     = ' , num2str(log10(VoE_phi(:, snrRunInd).'))]);
    disp(['Amplitude VoE (log10) at ', num2str(snr(snrRunInd)), ' dB SNR = ' , num2str(log10(VoE_rho(:, snrRunInd).'))]);
    fprintf('\n');
end

saveVoeData.phase     = VoE_phi;
saveVoeData.amplitude = VoE_rho;
saveFileFormat = '.mat';
saveFilePath   = 'Data/Output/ForPaper/MATFiles';
saveFileName   = [fullfile(saveFilePath,exptName), '_VoE', saveFileFormat];
save(saveFileName, "saveVoeData");

%% Generate Plots

% Setup class for plottingf
plt = classGenPlots_SL(lmc);

% Call function to produce plots
plt.genAllPlots();
