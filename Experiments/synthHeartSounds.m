dclose all
clearvars
clc

%% Setup

fs         = 48000;  % Hz
durA2      = 60; % 60 msec
durP2      = 60; % 60 msec
timeAxisA2 = 0:1/fs:(durA2/1000)-(1/fs);
timeAxisP2 = 0:1/fs:(durP2/1000)-(1/fs);

% Chirp signal setup
freqRangeA2 = [20, 100];
freqRangeP2 = [30, 120];
levelA2     = 0.60;
levelP2     = 0.40;

% Bandpass noise setup
levelBpNoise = 0.04;
bpFiltOrd    = 8;
freqRangeBpNoise = [40, 200];
varNoiseA2   = 0.02;
varNoiseP2   = 0.01;

% ADSR curve setup
targetA2    = [0.99999; 0.01; 0];
smoothingA2 = [0.08; 0.06; 0.05];
durationA2  = [18; 24; 18]*1; % msec

targetP2    = [0.99999; 0.01; 0];
smoothingP2 = [0.08; 0.06; 0.05];
durationP2  = [18; 24; 18]*1; % msec

% Additional setup
foreHold  = 5;    % msec
postHold  = 1200; % msec
splitDur  = 20;   % msec

% Time Axes
totalTime     = durA2/1000 + foreHold/1000 + postHold/1000;
totalTimeAxis = 0:1/fs:totalTime-1/fs;

%% Generate Component Signals

% Generate Chirps
chirpA2 = levelA2 .* chirp(timeAxisA2, freqRangeA2(2), (durA2/1000), freqRangeA2(1), 'logarithmic');
chirpP2 = levelP2 .* chirp(timeAxisA2, freqRangeP2(2), (durP2/1000), freqRangeP2(1), 'logarithmic');

% Generate ADSR curves
adsrForA2 = genAdsr(fs, targetA2, smoothingP2, durationA2);
adsrForP2 = genAdsr(fs, targetP2, smoothingP2, durationP2);

% Synthetic Heart Beat
modChirpA2 = adsrForA2.' .* chirpA2;
modChirpP2 = adsrForP2.' .* chirpP2;

modChirpS2 = modChirpA2 + modChirpP2;

%% Get the polynomial fit for the instantaneous frequency

polyOrd = 5;
 
fInst_A2 = freqRangeA2(2) * ((freqRangeA2(1) / freqRangeA2(2))^(1/(durA2/1000))).^timeAxisA2; 
fInst_P2 = freqRangeP2(2) * ((freqRangeP2(1) / freqRangeP2(2))^(1/(durP2/1000))).^timeAxisP2; 

fInst_A2_polyFit = polyfit(timeAxisA2, fInst_A2, polyOrd);
fInst_P2_polyFit = polyfit(timeAxisP2, fInst_P2, polyOrd);

fInst_A2_poly = polyval(fInst_A2_polyFit, timeAxisA2);
fInst_P2_poly = polyval(fInst_P2_polyFit, timeAxisP2);

phaseInst_A2 = cumtrapz(timeAxisA2, 2*pi*fInst_A2);
phaseInst_P2 = cumtrapz(timeAxisP2, 2*pi*fInst_P2);

phaseInst_A2_polyFit = polyfit(timeAxisA2, phaseInst_A2, polyOrd);
phaseInst_P2_polyFit = polyfit(timeAxisP2, phaseInst_P2, polyOrd);

phaseInst_A2_poly = polyval(phaseInst_A2_polyFit, timeAxisA2);
phaseInst_P2_poly = polyval(phaseInst_P2_polyFit, timeAxisP2);

%% Generate Polynoial phase cosine to compare with the original chirp

nA2 = 0:durA2 * fs / 1000 - 1;
chirpA2_pps    = levelA2 * genPolynomialPhaseSignal(nA2, fs, polyOrd, fliplr(phaseInst_A2_polyFit));
modChirpA2_pps = adsrForA2.' .* chirpA2_pps;

nP2 = 0:durP2 * fs / 1000 - 1;
chirpP2_pps    = levelP2 * genPolynomialPhaseSignal(nP2, fs, polyOrd, fliplr(phaseInst_P2_polyFit));
modChirpP2_pps = adsrForP2.' .* chirpP2_pps;

% %% Final Measured Chirp signal
% 
% snr = 12;
% modChirpS2Meas = addWhiteGaussianNoise(modChirpS2, snr);

%% Discrete Polynomial Transform

M = numel(phaseInst_P2_polyFit)-1;
N = length(timeAxisP2);
T = round(N/M); 

% Get analytic signal
analytic_chirpP2 = hilbert(chirpP2_pps);

% Compute DPT
nFft  = 2^nextpow2(N);
dw    = 2*pi/nFft;
wx    = 0:dw:2*pi-dw;
delta = 1/fs;
n   = 0:N-1;
a   = zeros(1,M);
dpt = zeros(M, nFft);
dp  = zeros(M, N);

s = analytic_chirpP2;

freqEstMethod = 'dft';

for m = M:-1:1

    % Compute DPT of order M
    [dpt(m,:), dp(m,:)] = discretePolynomialTransform(s, m, T, nFft);
    
    % Compute the coefficient
    if strcmpi(freqEstMethod, 'dft')
        [~, maxInd] = max(abs(dpt(m,:)));  
        wArgMax = wx(maxInd);
    elseif strcmpi(freqEstMethod, 'harmonicRet')
        [~, freqEst] = runFwdBwdMethod(dp(m,:).', N, 2, 'prony', 0.05);
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

a = flip(a);

%% Plots
close all 

figure(1)
subplot(2,2,1)
    plot(timeAxisA2 * 1000, fInst_A2, 'LineWidth', 1.1); hold on;
    plot(timeAxisA2 * 1000, fInst_A2_poly, '-.',  'LineWidth', 1.1);
    xlabel('time (ms)'); ylabel('Frequency');
    grid on; grid minor;
    legend('Log', 'PolyFit');
    title('Instantaneous Frequency for Synthetic S2');
subplot(2,2,3)
    plot(timeAxisP2 * 1000, fInst_P2,  'LineWidth', 1.1); hold on;
    plot(timeAxisP2 * 1000, fInst_P2_poly, '-.', 'LineWidth', 1.1);
    xlabel('time (ms)'); ylabel('Frequency');
    grid on; grid minor;
    legend('Log', 'PolyFit');
    title('Instantaneous Frequency for Synthetic P2');
subplot(2,2,2)
    plot(timeAxisA2 * 1000, phaseInst_A2, 'LineWidth', 1.1); hold on;
    plot(timeAxisA2 * 1000, phaseInst_A2_poly, '-.',  'LineWidth', 1.1);
    xlabel('time (ms)'); ylabel('Inst. Phase');
    grid on; grid minor;
    legend('Log', 'PolyFit');
    title('Instantaneous Phase for Synthetic S2');
subplot(2,2,4)
    plot(timeAxisP2 * 1000, phaseInst_P2,  'LineWidth', 1.1); hold on;
    plot(timeAxisP2 * 1000, phaseInst_P2_poly, '-.', 'LineWidth', 1.1);
    xlabel('time (ms)'); ylabel('Inst. Phase');
    grid on; grid minor;
    legend('Log', 'PolyFit');
    title('Instantaneous Phase for Synthetic P2');

figure(3)
subplot(2,1,1)
    plot(timeAxisA2 * 1000, chirpA2, 'LineWidth', 1.1); hold on;
    plot(timeAxisA2 * 1000, chirpA2_pps, 'LineWidth', 1.1);
    xlabel('time (ms)'); ylabel('Amplitude');
    grid on; grid minor;
    legend('A2', 'A2 PPS');
    title('Synthetic S2');
subplot(2,1,2)
    plot(timeAxisP2 * 1000, chirpP2, 'LineWidth', 1.1); hold on;
    plot(timeAxisP2 * 1000, chirpP2_pps, 'LineWidth', 1.1);
    xlabel('time (ms)'); ylabel('Amplitude');
    grid on; grid minor;
    legend('P2', 'P2 PPS');
    title('Synthetic P2');

figure(4)
subplot(3,1,1)
    plot(timeAxisA2 * 1000, adsrForA2,  'LineWidth', 1.1); hold on;
    plot(timeAxisA2 * 1000, modChirpA2, 'LineWidth', 1.1);
    plot(timeAxisA2 * 1000, modChirpA2_pps, 'LineWidth', 1.1);
    xlabel('time (ms)'); ylabel('Amplitude');
    grid on; grid minor;
    legend('ADSR', 'A2', 'A2 PPS');
    title('ADSR curve and Synthetic S2');
subplot(3,1,2)
    plot(timeAxisP2 * 1000, adsrForP2,  'LineWidth', 1.1); hold on;
    plot(timeAxisP2 * 1000, modChirpP2, 'LineWidth', 1.1);
    plot(timeAxisP2 * 1000, modChirpP2_pps, 'LineWidth', 1.1);
    xlabel('time (ms)'); ylabel('Amplitude');
    grid on; grid minor;
    legend('ADSR', 'P2', 'P2 PPS');
    title('ADSR curve and Synthetic P2');
subplot(3,1,3)
    plot(timeAxisP2 * 1000, modChirpS2,  'LineWidth', 1.1);
    xlabel('time (ms)'); ylabel('Amplitude');
    grid on; grid minor;
    title('Synthetic S2');

N = nFft;
fx = 0:fs/N:fs-fs/N;
figure(5)
tiledlayout flow
for m = M:-1:1
    nexttile
    plot(wx, abs(dpt(m,:)));
    title('Magnitude of Discrete Polynomial-Phase Transform');
    xlabel('Frequency Index'); ylabel('Magnitude');
    grid on; grid minor;
end

%% Helper

function [x, phs] = genPolynomialPhaseSignal(n, fs, K, phi)
    
    N = length(n);
    phs = zeros(1, N);
    for k = 0:K
        phs = phs + phi(k+1) * (n / fs).^k;
    end

    x = cos(phs);
end

function a = genAdsr(fs,target,gain,duration)
% Input
% target - vector of attack, sustain, release target values
% gain - vector of attack, sustain, release gain values
% duration - vector of attack, sustain, release durations in ms
% Output
% a - vector of adsr envelope values
    totalDur = sum(duration./1000);
    a = zeros(fs * totalDur,1); % assume 1 second duration ADSR envelope
    duration = round(duration./1000.*fs); % envelope duration in samp
    
    % Attack phase
    start = 2;
    stop  = duration(1);
    for n = start:stop
        a(n) = target(1)*gain(1) + (1.0 - gain(1))*a(n-1);
    end
    
    % Sustain phase
    start = stop + 1;
    stop  = start + duration(2);
    for n = start:stop
        a(n) = target(2)*gain(2) + (1.0 - gain(2))*a(n-1);
    end
    
    % Release phase
    start = stop + 1;
    stop  = totalDur * fs;
    for n = start:stop
        a(n) = target(3)*gain(3) + (1.0 - gain(3))*a(n-1);
    end
end

function [DPTM, DP] = discretePolynomialTransform(s, M, T, nFft)
    % s is the input signal
    % M is the order of the polynomial phase
    % T is the delay parameter
    
    N = length(s);

    if nargin < 4
        nFft = N;
    end

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
    DPTM = fft(DP, nFft) / nFft;
end
