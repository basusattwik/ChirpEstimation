clearvars
close all
clc

% Random seed
seed = 12345;
rng(seed);

%% Generate chirps

% Setup Chirp parameters
Fs = 1000;
Td = 0.5; % sec
L  = Td * Fs;

phi = [1, 0.5; ...
       0.2, 0.3]; % frequency and chirp rate
amp = [1; 1];

% Generate the chirp
[x, t] = genExpPolyChirp2(Fs, Td, amp, phi);

% Add White Gaussian Noise to x
snr = 3;
y   = addWhiteGaussianNoise(x, snr); 

%% Maximum Likelikhood Estimation 

% 1. Make the 2D Joint PDF g
% Create a 2D rectangular grid of M x M points alpha and beta
M   = 100;
dx  = 1/M;
rho = 0.4;
lam = 0.1;
P   = 2;
jointPdf = genPseudoPdf(x, phi, rho, M, Fs);

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
conBetaGivAlphaCdf = cumtrapz(conBetaGivAlphaPdf); % OR cumsum?
for i = 1:M
    conBetaGivAlphaCdf(i,:) = cumtrapz(conBetaGivAlphaPdf(i,:)) * dx;
end

% 6. Inverse CDF transform
R = 3000;

alphaGrid = linspace(0,1,M);
betaGrid  = linspace(0,1,M);

alphaIs = zeros(P,R);
betaIs  = zeros(P,R);
for r = 1:R

    % Sample uniform distribution
    u1 = rand(2,1);
    u2 = rand(2,1);

    % Get a sample of alpha
    temp = interp1(marAlphaCdf, alphaGrid, u1, 'spline');
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
        uBetaGrid   = betaGrid(uidx);

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
N = length(x) / Fs;
L = zeros(R, 1);
I = zeros(R, 1);
for r = 1:R
    L(r,1) = genLiklihoodFunc(N, Fs, P, alphaIs(:,r), beta(:,r), lam);
    c = genExpPolyChirp2(Fs, N, 1, -[alphaIs(:,r), beta(:,r)]);
    I(r,1) = exp(rho * (1/N) * abs(x.' * c)^2);
end




