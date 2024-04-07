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

f0_1 = 6; 
f1_1 = 6;
k_1  = (f1_1 - f0_1) / Td;

f0_2 = 17; 
f1_2 = 17;
k_2  = (f1_2 - f0_2) / Td;

phi = [f0_1, k_1; ...
       f0_2, k_2];
amp = [1, 1];

nPhi = size(phi, 1);

% Generate the chirp
[x, t] = genExpPolyChirp3(Fs, Td, amp, phi);

% Add White Gaussian Noise to x
snr = 30;
y   = addWhiteGaussianNoise(x, snr); 

%% Maximum Likelikhood Estimation 

% 1. Make the 2D Joint PDF g
% Create a 2D rectangular grid of M x M points alpha and beta
M    = 2000;
rho  = 4;
rho1 = 0.4;
P    = numel(amp);

alphaGrid = linspace(0,1,M);
betaGrid  = linspace(0,2,M);

[jointPdf,~] = genPseudoPdf(y, rho1, M);

% check that total area under pdf = 1
areaJointPdf = trapz(alphaGrid, trapz(betaGrid, jointPdf, 2));
disp(['total area of Joint PDF = ', num2str(areaJointPdf)]);

%%
% 2. Obtain the Marginal PDF of alpha

marAlphaPdf = trapz(betaGrid, jointPdf, 2);

% 3. Obtain the Marginal CDF of alpha
marAlphaCdf = cumtrapz(alphaGrid, marAlphaPdf) + 10e-5*rand(M,1); % adding random noise makes the CDF values "unique". Creates fewer issues during inv. CDF transformation.

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
conBetaGivAlphaCdf = conBetaGivAlphaCdf + 10e-5*rand(M,M);

%% 6. Inverse CDF transform

R = 10000;

alphaSampled = zeros(P,R);
betaSampled  = zeros(P,R);
for r = 1:R
    for p = 1:P
        % while(true)
            % Sample uniform distribution
            u1 = rand(1);
        
            % Get one sample estimate of alpha
            tempAlphaSample   = interp1(marAlphaCdf, alphaGrid, u1, 'linear');
            [~, alphaInd]     = min(abs(alphaGrid - tempAlphaSample));
            alphaSampled(p,r) = alphaGrid(alphaInd);
    
            % Choose alpha index to condition beta's inv. cdf generation on
            % alphaInd = find(alphaGrid >= alphaSampled(p,r), 1, "first");

            % if ~isempty(alphaInd)
            %     break % break the while loop
            % end
        % end
        
        % Sample uniform distribution 
        u2 = rand(1);

        % Get one sample estimate of beta
        tempBetaSample   = interp1(conBetaGivAlphaCdf(alphaInd, :), betaGrid, u2, 'linear'); 
        [~, betaInd]     = min(abs(betaGrid - tempBetaSample));
        betaSampled(p,r) = betaGrid(betaInd);
    end
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
tempAlphaSample = complex(0,0);
for p = 1:P   
    for r = 1:R

        expLMinusG = exp(L(r,1) - g(r,1));

        tempAlphaSample = exp(2*pi*1j * alphaSampled(p,r)) * expLMinusG;
        alphaEstArr(p,r) = tempAlphaSample;% / abs(temp);    

        tempAlphaSample = exp(pi*1j * betaSampled(p,r)) * expLMinusG;
        betaEstArr(p,r)  = tempAlphaSample;% / abs(temp);       
    end
end

alphaEst = (1/(2*pi)) * angle(mean(alphaEstArr, 2));
betaEst  = (2/(2*pi)) * angle(mean(betaEstArr,  2));

phiHat = [alphaEst * Fs, betaEst * 4*Fs^2];

%% Plots
close all

% Reconstruct signal
yHat = genExpPolyChirp3(Fs, Td, amp, phiHat);

% Time
tx = 0:1/Fs:Td-1/Fs;
figure(1)
plot(tx, imag(y)); hold on;
plot(tx, imag(yHat));
grid on; grid minor;
xlabel('Time (s)');
ylabel('Amp.');
legend('Original', 'Recon.');
title('Original vs Reconstructed superimposed chirps');

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
