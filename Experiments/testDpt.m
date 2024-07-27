close all
clearvars
clc

% Example usage
N = 600;               % Number of samples
A = 1;                 % Amplitude of the signal
fs = 133000;           % Sampling Frequency
delta = 1/fs;          % Sampling interval
coeff = [1.69, 3.45e5, -6e7, 8e9];  % Coefficients for the polynomial phase

M = numel(coeff)-1;                 % Order of the polynomial phase
T = round(N/M);                 % Delay parameter

% Generate the signal
s = generatePolyPhaseSignal(N, M, A, coeff, delta);
scopy = s;

%% Sequential DPT

nFft = 2^nextpow2(N);
dw   = 2*pi/nFft;
wx   = 0:dw:2*pi-dw;
n    = 0:N-1;

a   = zeros(M,1);
dpt = zeros(M, nFft);
dp  = zeros(M, N);
freqEstMethod = 'dft';
for m = M:-1:1

    % Compute DPT of order M
    [dpt(m,:), dp(m,:)] = discretePolynomialTransform(s, m, T);
    
    % Compute the coefficient
    if strcmpi(freqEstMethod, 'dft')
        [~, maxInd] = max(abs(dpt(m,:)));  
        wArgMax = wx(maxInd);
    else
    [~, freqEst] = runFwdBwdMethod(dp(m,:).', N, 1, 'prony', 0.05);
    wArgMax = 2*pi*freqEst{1};
    end

    if wArgMax < pi
        a(m) = wArgMax / (factorial(m) * (T * delta)^(m-1)) * fs;
    else
        a(m) = (wArgMax - 2*pi) / (factorial(m) * (T * delta)^(m-1)) * fs;
    end

    % Peel off
    s = s .* exp(-1j * a(m) * (n * delta).^m);
end

%% Plots

figure(1)
plot(real(s)); %hold on;
% plot(imag(s));
grid on; grid minor;
xlabel('time(s)'); ylabel('Amplitude');
title('A Polynomial Phase Signal');

% Plot the magnitude of the DPTM
fx = 0:fs/N:fs-fs/N;
figure(2);
tiledlayout flow
for m = M:-1:1
    nexttile
    plot(wx, abs(dpt(m,:)));
    title('Magnitude of Discrete Polynomial-Phase Transform');
    xlabel('Frequency Index'); ylabel('Magnitude');
    grid on; grid minor;
end

%% Helpers

function [DPTM, DP] = discretePolynomialTransform(s, M, T)
    % s is the input signal
    % M is the order of the polynomial phase
    % T is the delay parameter
    
    N  = length(s);
    DP = s;  % Start with the first-order operator

    % Compute higher order DP operators recursively
    for m = 2:M
        tempDP = zeros(size(DP));
        for n = T+1:N
            tempDP(n) = DP(n) * conj(DP(n-T));
        end
        DP = tempDP;
    end

    % Compute the DPTM which is the DTFT of the last DP operator
    DPTM = fft(DP, 2^nextpow2(N));
end

function s = generatePolyPhaseSignal(N, M, A, coeff, delta)
    % N is the number of samples
    % M is the polynomial order
    % A is the amplitude of the signal
    % coeff is the vector of coefficients for the polynomial phase (a_m)
    % delta is the sampling interval
    
    n = 0:N-1;  % Time index
    phasePoly = zeros(1,N);
    
    for m = 0:M
        phasePoly = phasePoly + coeff(m+1) * (n * delta).^m;
    end
    
    s = A * exp(1j * phasePoly);  % Generate the complex signal
end
