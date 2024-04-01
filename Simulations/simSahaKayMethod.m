clearvars
close all
clc

% Random seed
seed = 12345;
rng(seed);

%% Generate chirps

% Setup Chirp parameters
Fs = 50;
Td = 0.5; % sec

f0  = 5; 
f1  = 20;
k   = (f1 - f0) / Td;
phi = [f0, k; ...
       0, 0];
amp = [1, 0];

nPhi = size(phi, 1);

% Generate the chirp
[x, t] = genExpPolyChirp3(Fs, Td, amp, phi);

% Add White Gaussian Noise to x
snr = 120;
y   = addWhiteGaussianNoise(x, snr); 

%% Maximum Likelikhood Estimation 

% 1. Make the 2D Joint PDF g
% Create a 2D rectangular grid of M x M points alpha and beta
M   = 2000;
dx  = 1/M;
rho = 0.8;
lam = 4;
P   = numel(amp);
jointPdf = genPseudoPdf(y, nPhi, rho, M);

% % load joint pdf
% filePath = '/Users/sattwikbasu/Documents/MATLAB/ProjectData/ChirpParameterEstimation/';
% fileName   = ['jointPdf_M_', num2str(M), '_Fs_', num2str(Fs), '.mat'];
% dataStruct = load([filePath, fileName]);
% jointPdf   = real((dataStruct.jointPdf).^rho);

% 2. Obtain the Marginal PDF of alpha
marAlphaPdf = mean(jointPdf, 2);

% 3. Obtain the Marginal CDF of alpha
marAlphaCdf = cumtrapz(marAlphaPdf) * dx; % OR cumsum?

% 4. Obtain Conditional PDF of beta when given alpha
conBetaGivAlphaPdf = zeros(size(jointPdf));
for i = 1:M
    for j = 1:M
        conBetaGivAlphaPdf(i,j) = jointPdf(i,j) / marAlphaPdf(i);
    end
end

% 5. Obtain Conditional CDF of beta when given alpha
conBetaGivAlphaCdf = zeros(size(jointPdf)); %cumtrapz(conBetaGivAlphaPdf); % OR cumsum?
for i = 1:M
    conBetaGivAlphaCdf(i,:) = cumtrapz(conBetaGivAlphaPdf(i,:)) * dx;
end

%% 6. Inverse CDF transform
R = 4000;

alphaGrid = linspace(0,1,M);
betaGrid  = linspace(0,2,M);

alphaIs = zeros(P,R);
betaIs  = zeros(P,R);
for r = 1:R

    % Sample uniform distribution
    u1 = rand(P,1);
    u2 = rand(P,1);

    % Get a sample of alpha
    [uMarAlphaCdf, uidx] = unique(marAlphaCdf);
    uAlphaGrid = alphaGrid(uidx);
    temp = interp1(uMarAlphaCdf, uAlphaGrid, u1, 'spline');
    if temp < 0
        temp = 0;
    elseif temp > 1
        temp = 1;
    end
    alphaIs(:,r) = temp;

    % Get a sample of beta
    for p = 1:P
        alphaInd = find(alphaGrid <= alphaIs(p,r), 1, "last");
        [uConBetaGivAlphaCdf, uidx] = unique(conBetaGivAlphaCdf(alphaInd, :));
        uBetaGrid = betaGrid(uidx);

        temp = interp1(uConBetaGivAlphaCdf,  uBetaGrid, u2(p,1), 'spline');
        if temp < 0
            temp = 0;
        elseif temp > 1
            temp = 1;
        end
        betaIs(p,r) = temp;
    end
end

% Calculate circular means
L = zeros(R,1);
g = zeros(R,1);
for r = 1:R
    L(r,1) = genLiklihoodFunc(y, P, alphaIs(:,r), betaIs(:,r), lam);
    g(r,1) = genProposalFunc(y, P, alphaIs(:,r), betaIs(:,r), rho);
end

arg = complex(0,0); % Init

% Frequencies
alphaEst = zeros(P,1);
for p = 1:P   
    for r = 1:R
        arg = arg + exp(2*pi*1i * alphaIs(p,r)) * exp(L(r,1) - g(r,1));       
    end
    alphaEst(p,1) = (1/(2*pi)) * angle((1/R) * arg) * Fs;
    arg(:) = 0 + 0*1j;
end

% Chirp rate
betaEst = zeros(P,1);
for p = 1:P   
    for r = 1:R
        arg = arg + exp(pi*1i * betaIs(p,r)) * exp(L(r,1) - g(r,1));       
    end
    betaEst(p,1) = (1/pi) * angle((1/R) * arg) * Fs;
    arg(:) = 0 + 0*1j;
end


