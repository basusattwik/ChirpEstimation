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
M   = 10;
dx  = 1/M;
rho = 0.4;
P   = 2;
jointPdf = genPseudoPdf(x, phi, rho, M, Fs);

% 2. Obtain the Marginal PDF of alpha
marAlphaPdf = mean(jointPdf, 2);

% 3. Obtain the Marginal CDF of alpha
marAlphaCdf = cumtrapz(marAlphaPdf) * dx; % OR cumsum?

% 4. Obtain Conditional PDF of beta when given alpha
conBetaPdf = zeros(size(jointPdf));
for i = 1:M
    for j = 1:M
        conBetaPdf(i,j) = jointPdf(i,j) / marAlphaPdf(i);
    end
end

% 5. Obtain Conditional CDF of beta when given alpha
conBetaCdf = cumtrapz(conBetaPdf); % OR cumsum?

% 6. Inverse CDF transform
u1 = unifrnd(0, 1, [P,1]);
u2 = unifrnd(0, 1, [P,1]);

invMarAlphaCdf = interp1(marAlphaCdf, linspace(0,1,M) , u1, 'spline','extrap');
% invConAlphaCdf = interp1(conBetaCdf,  linspace(0,1,M) , u2, 'spline','extrap');


