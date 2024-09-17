clear all
clc

%% 

% Setup Chirp parameters
fs = 500;
Nc = 2;
Td = 1;
N  = Td * fs;

alpha = [1, 1];
beta  = [5, 7];
gamma = [5, 4];

phi   = cell(1,Nc);
phi{1,1} = [0, 50, 2].'; 
phi{1,2} = [0, 70, -2].'; 

% How many params?
numParams = 0;
for c = 1:Nc
    numParams = numParams + numel(phi{1,c});
end
numParams = numParams + numel(beta) + numel(gamma);
snr = 120;

% Tuning for LMC
stepSizeBeta  = [0.1, 0.1].'; % should be different for different params
stepSizeGamma = [0.1, 0.1].';
stepSizePhi   = [0.00, 0.09, 0.01, 0.0, 0.09, 0.01].';
stepSize = [stepSizeBeta; stepSizeGamma; stepSizePhi];
numIter  = 10000;
tempSwap = 100;
initInvTemp = 2;
numSamplesToUse = 20;
avgConst = 0.9;
stepSizeConst = 0.02;
stepSizeMax  = 0.9;
numParticles = 20;

%% Build setup structure

cpeSetting.fs    = fs;
cpeSetting.Td    = Td;
cpeSetting.Nc    = Nc;  
cpeSetting.alpha = alpha;
cpeSetting.beta  = beta;
cpeSetting.gamma = gamma;
cpeSetting.phi   = phi;
cpeSetting.snr   = snr;
cpeSetting.numParams = numParams;

tune.stepSize      = stepSize;
tune.numIter       = numIter;
tune.initInvTemp   = initInvTemp;
tune.tempSwap      = tempSwap;
tune.stepSizeConst = stepSizeConst;
tune.stepSizeMax   = stepSizeMax;
tune.avgConst      = avgConst;
tune.numParticles  = numParticles;
tune.numSamplesToUse = numSamplesToUse;

%% Simulation

lmc = classLangevinMonteCarlo(tune, cpeSetting);
lmc = lmc.runLmcCore();

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
%     plot(1:numIter, lmc.saveInvTemp);
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