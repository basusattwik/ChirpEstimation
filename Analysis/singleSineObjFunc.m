clc
% close all
clearvars

% This is not correct. Noise should most likely not be added to the
% measured signal. 

% It should be added during the LMC updates to the parameters
% This has the effect of Gaussian smoothing the gradients
%% Create a chirp

fs = 2000;
Td = 0.05;
Nc = 1; 
snr = 6;
N = Td * fs;

% Parameters
f0 = 60;
tind  = (0:1/fs:Td-1/fs).';

x  = exp(2*pi*1j * f0*tind);
y = addGaussianNoise(x, snr);

% 1. Grid of parameters
df = 0.1;
f  = 0:df:100;

J  = zeros(numel(f),1);
dJ = zeros(numel(f),1);
ddJ_tr = zeros(numel(f),1);
for find = 1:numel(f)

    H  = exp(2*pi*1j * f(find) * tind);
    P  = H * (H' * H)^(-1) * H';
    Po = eye(N,N) - P;
    dH = 2*pi*1j .* tind .* H;
    J(find)  = real(y' * (Po * y));
    dJ(find) = -2*real(y' * (Po * (dH * ((H' * H)^(-1) * (H' * y)))));
end

h = 1;
ddJ = zeros(numel(f),1);
for find = 2:numel(f)-1
    ddJ(find) = (J(find+1)-2*J(find)+J(find-1))/df^2;
end

%% Plot

sigma_crude = 100;
L = length(f);
alpha = (L - 1)/(2*sigma_crude);

% Filter gradient
w = gausswin(L, alpha);
w = w / sum(w);

dJ_filt = conv(w, dJ, "same");
J_filt  = conv(w, J,  "same");

%% Actual Smoothing algorithm

sigma_gauss  = 1;

numSmoothing = 50;

J_gs  = zeros(numel(f),1);
dJ_gs = zeros(numel(f),1);
ddJ_gs = zeros(numel(f),1);

J_temp_noisy   = 0;
dJ_temp_noisy  = 0;
ddJ_temp_noisy = 0;
I = eye(N);
% 
for find = 1:numel(f)

    for sind = 1:numSmoothing

        u = randn(1,1);
        f_noisy = f(find) + sigma_gauss * u;

        H_noisy  = exp(2*pi*1j * f_noisy * tind);
        P_noisy  = H_noisy * (H_noisy' * H_noisy)^(-1) * H_noisy';
        Po_noisy = I - P_noisy;
        dH_noisy = 2*pi*1j .* tind .* H_noisy;

        % Function at the perturbed point
        J_noisy = real(y' * (Po_noisy * y));
        J_temp_noisy = J_temp_noisy + J_noisy;

        % Gradient at perturbed point
        dJ_temp_noisy = dJ_temp_noisy + -2*real(y' * (Po_noisy * (dH_noisy * ((H_noisy' * H_noisy)^(-1) * (H_noisy' * y)))));

        % Hessian at perturbed point
        ddJ_temp_noisy = ddJ_temp_noisy + (1/sigma_gauss^2) * (u^2 - 1) * J_noisy;
    end

    J_gs(find)   = J_temp_noisy   / numSmoothing;
    dJ_gs(find)  = dJ_temp_noisy  / numSmoothing;
    ddJ_gs(find) = ddJ_temp_noisy / numSmoothing;

    J_temp_noisy   = 0;
    dJ_temp_noisy  = 0;
    ddJ_temp_noisy = 0;
end

% for find = 1:numel(f)
% 
%     f_curr  = f(find);
%     H  = exp(2*pi*1j * f_curr * tind);
%     P  = H * (H' * H)^(-1) * H';
%     Po = I - P;
%     dH = 2*pi*1j .* tind .* H;
% 
%     J_curr = real(y' * (Po * y));
% 
%     for sind = 1:numSmoothing
% 
%         u = randn(1,1);
%         f_noisy = f(find) + sigma_gauss * u;
% 
%         H_noisy  = exp(2*pi*1j * f_noisy * tind);
%         P_noisy  = H_noisy * (H_noisy' * H_noisy)^(-1) * H_noisy';
%         Po_noisy = I - P_noisy;
%         dH_noisy = 2*pi*1j .* tind .* H_noisy;
% 
%         % Function at the perturbed point
%         J_noisy = real(y' * (Po_noisy * y));
%         J_temp_noisy = J_temp_noisy + J_noisy;
% 
%         % Gradient at perturbed point
%         dJ_temp_noisy = dJ_temp_noisy + -2*real(y' * (Po_noisy * (dH_noisy * ((H_noisy' * H_noisy)^(-1) * (H_noisy' * y)))));
% 
%         % Hessian at perturbed point
%         ddJ_temp_noisy = ddJ_temp_noisy + (1/sigma_gauss^2) * (u^2 - 1) * (J_noisy - J_curr);
%     end
% 
%     J_gs(find)   = J_temp_noisy   / numSmoothing;
%     dJ_gs(find)  = dJ_temp_noisy  / numSmoothing;
%     ddJ_gs(find) = ddJ_temp_noisy / numSmoothing;
% 
%     J_temp_noisy   = 0;
%     dJ_temp_noisy  = 0;
%     ddJ_temp_noisy = 0;
% end

% Plots

% close all
figure('windowstyle','docked');
    plot(f,w);

figure('windowstyle','docked');
tiledlayout flow

ax1 = nexttile;
    plot(f, J,   'LineWidth', 1.5); hold on
    plot(f, J_gs, 'LineWidth',1.5)
    grid on; grid minor;
    xlabel('Frequency (Hz)');
    ylabel('Objective Function');
    legend('Actual', 'Gaussian Smoothed');
    title('Objective Function for a sine tone at 60 Hz');

ax2 = nexttile;
    plot(f, dJ,  'LineWidth', 1.5); hold on;
    plot(f, dJ_gs, 'LineWidth',1.5)
    xlabel('Frequency (Hz)');
    ylabel('Gradient');
    title('Gradient for a sine tone at 60 Hz');
    legend('Actual', 'Gaussian Smoothed');
    grid on; grid minor;

ax3 = nexttile;
    plot(f, ddJ,  'LineWidth', 1.5); hold on;
    plot(f, (ddJ_gs), 'LineWidth',1.5)
    xlabel('Frequency (Hz)');
    ylabel('Hessian');
    title('Hessian for a sine tone at 60 Hz');
    legend('Actual', 'Gaussian Smoothed');
    grid on; grid minor;
linkaxes([ax1, ax2, ax3]);

%% Helpers

function ym = addGaussianNoise(xm, snr)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    
    if ~isreal(xm)
        wm = randn(size(xm)) + 1j*randn(size(xm)); %  random noise signal
    else
        wm = randn(size(xm)); %  random noise signal
    end
    
    xpow = sum(abs(xm.^2)) ;
    wpow = sum(abs(wm.^2));   
    pn   = xpow / (10^(snr / 10)); % required SNR
    wsc  = (wm ./ sqrt(wpow)) .* sqrt(pn); % normalizing the noisy signal to get unity power then multiplying it by the new power to achieve the required SNR    
    ym   = xm + wsc; % add noise to signal
    
    % wnewpow = sum(abs(wscaled.^2));
    % snrtrue = 10*log10(ypow / wnewpow);
end