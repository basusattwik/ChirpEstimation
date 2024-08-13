clear all
clc

%% 

% Setup Chirp parameters
fs = 200;
Td = 1;
N  = Td * fs;

alpha = [1,1].';
phi{1,1} = [0, 30.6, 50].'; 
phi{1,2} = [0, 20.5, 32].'; 
snr = 20;

% Tuning for LMC
numParticles  = 20;
stepSizePhi   = [0.00001, 0.00001, 0.00001].';
stepSizeConst = 10;
stepSizeMax   = 1;
avgConst      = 0.97;
tempConst     = 2;
numIterNoise  = 20;  % This is L
numIterLmc    = 1000; % This is T
noiseVarInit  = 2;
noiseVarFinal = 0.00001;

%% Build setup & tuning structures

% How many params?
Nc = length(phi); % Number of chirps
numParams = 0;
for c = 1:Nc
    numParams = numParams + numel(phi{1,c});
end

stepSize = repmat(stepSizePhi, Nc, 1);

cpeSetting.fs        = fs;
cpeSetting.Td        = Td;
cpeSetting.Nc        = Nc;  
cpeSetting.alpha     = alpha;
cpeSetting.phi       = phi;
cpeSetting.snr       = snr;
cpeSetting.numParams = numParams;

lmcTuning.stepSize      = stepSize;
lmcTuning.tempConst     = tempConst;
lmcTuning.stepSizeConst = stepSizeConst;
lmcTuning.stepSizeMax   = stepSizeMax;
lmcTuning.avgConst      = avgConst;
lmcTuning.numParticles  = numParticles;
lmcTuning.numIterNoise  = numIterNoise;
lmcTuning.numIterLmc    = numIterLmc;
lmcTuning.noiseVarInit  = noiseVarInit;
lmcTuning.noiseVarFinal = noiseVarFinal;

%% Simulation

% Call the constructor to setup
lmc = classLangevinMonteCarlo(lmcTuning, cpeSetting);

% Start the simulation
tic;
lmc = lmc.runLmcCore();
toc;

%% Generate Plots

% figure('windowstyle','docked')
% tiledlayout flow
% nexttile
%     plot(1:N, real(lmc.cpe.ym));
%     xlabel('samples'); 
%     ylabel('Amplitude');
%     grid on; grid minor;
%     title('Mixture of Chirps i.e., measured signal');
% nexttile
%     plot(1:numIter, lmc.saveObjectiveFunc);
%     xlabel('Iterations');
%     ylabel('Objective Function');
%     grid on; grid minor;
%     title('Objective Function vs Iterations');
% nexttile
%     stairs(1:numIter, lmc.savebAccept);
%     xlabel('Iterations');
%     yticks([0 1]);
%     ylim([-0.5 1.5]);
%     yticklabels({'Reject', 'Accept'});
%     grid on; grid minor;
%     title('Accept or Reject');
% nexttile
%     plot(1:numIter, lmc.saveGradNorm);
%     xlabel('Iterations');
%     ylabel('Norm of the Gradient');
%     grid on; grid minor;
%     title('Gradient Norm');
% nexttile
%     plot(1:numIter, lmc.saveTemp);
%     xlabel('Iterations');
%     ylabel('Inverse Temp');
%     grid on; grid minor;
%     title('Inverse Temperature');
% nexttile
%     plot(1:numIter, lmc.saveStepSize);
%     xlabel('Iterations');
%     ylabel('Stepsize');
%     grid on; grid minor;
%     title('Stepsize');
