clear all
clc

%% 

% Setup Chirp parameters
fs = 500;
Nc = 1;
Td = 1;
N  = Td * fs;

alpha = [1];

phi = cell(1,Nc);
phi{1,1} = [0, 30].'; 

% How many params?
numParams = 0;
for c = 1:Nc
    numParams = numParams + numel(phi{1,c});
end
snr = 6;

% Tuning for LMC

stepSizePhi     = [0.0, 0.6].';
stepSize        = stepSizePhi;
numIter         = 10000;
tempSwap        = 70;
initTemp        = 10; % T in exp(-1/T * F(theta)), start high and then drop
numSamplesToUse = 20;
avgConst        = 0.99;
stepSizeConst   = 0.1;
stepSizeMax     = 2;
numParticles    = 1;

%% Build setup structure

cpeSetting.fs    = fs;
cpeSetting.Td    = Td;
cpeSetting.Nc    = Nc;  
cpeSetting.alpha = alpha;
cpeSetting.phi   = phi;
cpeSetting.snr   = snr;
cpeSetting.numParams = numParams;

tune.stepSize      = stepSize;
tune.numIter       = numIter;
tune.initTemp      = initTemp;
tune.tempSwap      = tempSwap;
tune.stepSizeConst = stepSizeConst;
tune.stepSizeMax   = stepSizeMax;
tune.avgConst      = avgConst;
tune.numParticles  = numParticles;
tune.numSamplesToUse = numSamplesToUse;

%% Simulation


lmc = classLangevinMonteCarlo(tune, cpeSetting);

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