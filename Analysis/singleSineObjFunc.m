clc
% close all
clearvars

% This is not correct. Noise should most likely not be added to the
% measured signal. 

% It should be added during the LMC updates to the parameters
% This has the effect of Gaussian smoothing the gradients
%% Create a chirp

fs = 300;
Td = 1;
Nc = 1; 
snr = 60;
N = Td * fs;

% Parameters
f0 = 60;
tind  = (0:1/fs:Td-1/fs).';

x  = exp(2*pi*1j * f0*tind);
y = addGaussianNoise(x, snr);

% 1. Grid of parameters
df = 0.01;
f  = 50:df:70;

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

I = eye(N);
numSmoothing = 50;

J_gs  = zeros(numel(f),1);
dJ_gs = zeros(numel(f),1);

J_temp  = 0;
dJ_temp = 0;
sigma_gauss = 2;
for find = 1:numel(f)

    for sind = 1:numSmoothing

        f_noisy = f(find) + sigma_gauss * randn(1,1);

        H  = exp(2*pi*1j * f_noisy * tind);
        P  = H * (H' * H)^(-1) * H';
        Po = I - P;
        dH = 2*pi*1j .* tind .* H;

        J_temp = J_temp + real(y' * (Po * y));
        dJ_temp = dJ_temp + -2*real(y' * (Po * (dH * ((H' * H)^(-1) * (H' * y)))));

    end
    J_gs(find)  = J_temp / numSmoothing;
    dJ_gs(find) = dJ_temp / numSmoothing;
    dJ_temp = 0;
    J_temp = 0;
end

% Plots

close all
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
linkaxes([ax1, ax2]);

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