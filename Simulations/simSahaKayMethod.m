clearvars
close all
clc

% Random seed
seed = 12345;
rng(seed);

%% Generate chirps

% Setup Chirp parameters
Fs = 50;
Td = 1; % sec

f0_1 = 5; 
f1_1 = 15;
k_1  = (f1_1 - f0_1) / Td;

f0_2 = 7; 
f1_2 = 21;
k_2  = (f1_2 - f0_2) / Td;

phi = [f0_1, k_1; ...
       f0_2, k_2];
amp = [1, 1];

% phi  = [f0_1, k_1];
% amp  = 1;
nPhi = size(phi, 1);

M    = 2000;
rho  = 7;
rho1 = 0.4;

%% Generate the chirp

[x, t] = genExpPolyChirp3(Fs, Td, amp, phi);

% Add White Gaussian Noise to x
snr = 80;
y   = addWhiteGaussianNoise(x, snr); 

%% Maximum Likelikhood Estimation 

% 1. Make the 2D Joint PDF g
% Create a 2D rectangular grid of M x M points alpha and beta
P = numel(amp);

alphaGrid = linspace(0,1/2,M);
betaGrid  = linspace(0,2/2,M);

[jointPdf,~] = genPseudoPdf(y, rho1, M);

% check that total area under pdf = 1
areaJointPdf = trapz(alphaGrid, trapz(betaGrid, jointPdf, 2));
disp(['total area of Joint PDF = ', num2str(areaJointPdf)]);

%%
% 2. Obtain the Marginal PDF of alpha

marAlphaPdf = trapz(betaGrid, jointPdf, 2);

% 3. Obtain the Marginal CDF of alpha
marAlphaCdf = cumtrapz(alphaGrid, marAlphaPdf) + 10e-6*rand(M,1); % adding random noise makes the CDF values "unique". Creates fewer issues during inv. CDF transformation.

% 4. Obtain Conditional PDF of beta when given alpha
conBetaGivAlphaPdf = zeros(size(jointPdf));
for i = 1:M
    for j = 1:M
        conBetaGivAlphaPdf(i,j) = jointPdf(i,j) / marAlphaPdf(i);
    end
end

% 5. Obtain Conditional CDF of beta when given alpha
conBetaGivAlphaCdf = zeros(size(jointPdf));
for i = 1:M
    conBetaGivAlphaCdf(i,:) = cumtrapz(betaGrid, conBetaGivAlphaPdf(i,:), 2);
end
conBetaGivAlphaCdf = conBetaGivAlphaCdf + 10e-6*rand(M,M);

%% 6. Inverse CDF transform

R = 1000;

alphaSampled  = zeros(P,R);
betaSampled   = zeros(P,R);
tempAlphaGrid = alphaGrid;
for r = 1:R
    % Sample uniform distribution
    u1 = rand(P,1);
    u2 = rand(P,1);
    for p = 1:P
    
        % Get one sample estimate of alpha
        tempAlphaSample   = interp1(marAlphaCdf, alphaGrid, u1(p,1), 'linear');
        [~, alphaInd]     = min(abs(tempAlphaGrid - tempAlphaSample));
        alphaSampled(p,r) = alphaGrid(alphaInd);

        % Get one sample estimate of beta
        tempBetaSample   = interp1(conBetaGivAlphaCdf(alphaInd, :), betaGrid, u2(p,1), 'linear'); 
        [~, betaInd]     = min(abs(betaGrid - tempBetaSample));
        betaSampled(p,r) = betaGrid(betaInd);

        % Save previous alphaInd
        tempAlphaGrid(alphaInd) = [];
    end
    tempAlphaGrid = alphaGrid;
end

% Calculate circular means
L = zeros(R,1);
g = zeros(R,1);
for r = 1:R
    L(r,1) = genLiklihoodFunc(y, P, alphaSampled(:,r), betaSampled(:,r), rho);
    g(r,1) = genProposalFunc(y,  P, alphaSampled(:,r), betaSampled(:,r), rho1);
end

% Frequencies & Chirp rate Circular Means
alphaEstArr = complex(zeros(P,R), zeros(P,R));
betaEstArr  = complex(zeros(P,R), zeros(P,R));
for p = 1:P   
    for r = 1:R
        expLMinusG       = exp(L(r,1) - g(r,1));
        alphaEstArr(p,r) = exp(2*pi*1j * alphaSampled(p,r)) * expLMinusG;
        betaEstArr(p,r)  = exp(pi*1j   * betaSampled(p,r))  * expLMinusG;     
    end
end

alphaEst = (1/(2*pi)) * angle(mean(alphaEstArr, 2));
betaEst  = (2/(2*pi)) * angle(mean(betaEstArr,  2));

phiHat = [alphaEst * Fs, betaEst * Fs^2];

%% Plots
close all

% Reconstruct signal
yHat = genExpPolyChirp3(Fs, Td, amp, phiHat);

x1 = genExpPolyChirp3(Fs, Td, amp(1), phi(1,:));
x2 = genExpPolyChirp3(Fs, Td, amp(2), phi(2,:));

x1Hat = genExpPolyChirp3(Fs, Td, amp(1), phiHat(1,:));
x2Hat = genExpPolyChirp3(Fs, Td, amp(2), phiHat(2,:));

% Time
tx = 0:1/Fs:Td-1/Fs;
figure(1)
subplot(3,1,1)
    plot(tx, imag(y)); hold on;
    plot(tx, imag(yHat));
    grid on; grid minor;
    xlabel('Time (s)');
    ylabel('Amp.');
    legend('Original', 'Recon.');
    title('$x_1 + x_2$', 'interpreter', 'latex');
subplot(3,1,2)
    plot(tx, imag(x1)); hold on;
    plot(tx, imag(x2Hat));
    grid on; grid minor;
    xlabel('Time (s)');
    ylabel('Amp.');
    legend('Original', 'Recon.');
    title('$x_1$', 'interpreter', 'latex');
subplot(3,1,3)
    plot(tx, imag(x2)); hold on;
    plot(tx, imag(x1Hat));
    grid on; grid minor;
    xlabel('Time (s)');
    ylabel('Amp.');
    legend('Original', 'Recon.');
    title('$x_2$', 'interpreter', 'latex');
sgtitle('Original vs Reconstructed superimposed chirps');

figure(2)
subplot(2,1,1)
    plot(alphaGrid * Fs, marAlphaPdf);
    xlabel('$\alpha$','interpreter','latex');
    ylabel('$g(\alpha)$', 'interpreter', 'latex');
    grid on; grid minor;
    title('PDF of frequency parameter');
subplot(2,1,2)
    plot(alphaGrid * Fs, marAlphaCdf);
    xlabel('$\alpha$','interpreter','latex');
    ylabel('$G(\alpha)$', 'interpreter', 'latex');
    grid on; grid minor;
    title('CDF of frequency parameter');

figure(3)
surf(betaGrid * Fs^2, alphaGrid * Fs, jointPdf, 'EdgeColor', 'none');
ylabel('$\alpha$','interpreter','latex');
xlabel('$\beta$','interpreter','latex');
zlabel('$g(\alpha, \beta)$', 'interpreter', 'latex');
title('Joint PDF');

figure(4)
surf(betaGrid * Fs^2, alphaGrid * Fs, conBetaGivAlphaPdf, 'EdgeColor', 'none');
ylabel('$\alpha$','interpreter','latex');
xlabel('$\beta$','interpreter','latex');
zlabel('$g(\beta | \alpha)$', 'interpreter', 'latex');
title('Conditional PDF');
