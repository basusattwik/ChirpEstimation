close all
clearvars
clc

%IDEA: Create chirp pulses modulated by an ADSR curve

%% Setup

fs = 3000;  % Hz
durA2 = 60; % msec
durP2 = 60; % msec
timeAxisA2 = 0:1/fs:(durA2/1000)-1/fs;
timeAxisP2 = 0:1/fs:(durP2/1000)-1/fs;

% Chirp signal setup
freqRangeA2 = [20, 100];
freqRangeP2 = [40, 140];
levelA2 = 0.60;
levelP2 = 0.40;

% Bandpass noise setup
levelBpNoise = 0.04;
bpFiltOrd    = 8;
freqRangeBpNoise = [40, 200];
varNoiseA2 = 0.02;
varNoiseP2 = 0.01;

% ADSR curve setup
targetA2    = [0.99999; 0.01; 0];
smoothingA2 = [0.08; 0.06; 0.05];
durationA2  = [18; 24;18]; %msec

targetP2    = [0.99999; 0.01; 0];
smoothingP2 = [0.08; 0.06; 0.05];
durationP2  = [18; 24;18]; %msec

% Additional setup
foreHold  = 5;    % msec
postHold  = 1200; % msec
splitDur  = 20;   % msec

% Time Axes
totalTime = durA2/1000 + foreHold/1000 + postHold/1000;
totalTimeAxis = 0:1/fs:totalTime-1/fs;

%% Generate Component Signals

% Generate Chirps
chirpA2 = levelA2 .* chirp(timeAxisA2, freqRangeA2(2), (durA2/1000), freqRangeA2(1), 'logarithmic');
chirpP2 = levelP2 .* chirp(timeAxisA2, freqRangeP2(2), (durP2/1000), freqRangeP2(1), 'logarithmic');

% Generate ADSR curves
adsrForA2 = genAdsr(fs, targetA2, smoothingP2, durationA2);
adsrForP2 = genAdsr(fs, targetP2, smoothingP2, durationP2);

% Add some bandpass filtered noise
[b,a]     = butter(bpFiltOrd, [freqRangeBpNoise(1)/(fs * 0.5), freqRangeBpNoise(2)/(fs * 0.5)], 'bandpass');
bpNoiseA2 = levelBpNoise * filter(b, a, varNoiseA2 * randn(1, durA2/1000*fs));
bpNoiseP2 = levelBpNoise * filter(b, a, varNoiseP2 * randn(1, durP2/1000*fs));

% Synthetic Heart Beat
modChirpA2 = adsrForA2.' .* (chirpA2 + bpNoiseA2);
modChirpP2 = adsrForA2.' .* (chirpP2 + bpNoiseP2);

fullEventA2 = [zeros(1,round(foreHold*fs/1000)), modChirpA2, zeros(1,round(postHold*fs/1000))];
fullEventP2 = [zeros(1,round(foreHold*fs/1000)), modChirpP2, zeros(1,round(postHold*fs/1000))];

repeatedEventA2 = repmat(fullEventA2, [1, 20]);
repeatedEventP2 = repmat(fullEventP2, [1, 20]);

additiveNoise = 0.005*randn(1, length(repeatedEventP2));
syntheticS2   = levelA2 .* repeatedEventA2 + delayseq(levelP2 .* repeatedEventP2.', splitDur/1000, fs).' + additiveNoise;

% Repeat the same thing many times
totalDuration = length(syntheticS2) / fs;
totalRepeatedTimeAxis = 0:1/fs:totalDuration-1/fs;

% Compute CWT
[cwtr, f] = cwt(syntheticS2, fs, 'VoicePerOctave', 32, 'amor');
scalogrm  = abs(cwtr);

%% Plots 

figure(1)
subplot(2,1,1)
    plot(timeAxisA2 * 1000, adsrForA2,  'LineWidth', 1.1); hold on;
    plot(timeAxisA2 * 1000, modChirpA2, 'LineWidth', 1.1);
    xlabel('time (ms)'); ylabel('Amplitude');
    grid on; grid minor;
    legend('ADSR', 'A2');
    title('ADSR curve and Synthetic S2')
subplot(2,1,2)
    plot(timeAxisP2 * 1000, adsrForP2,  'LineWidth', 1.1); hold on;
    plot(timeAxisP2 * 1000, modChirpP2, 'LineWidth', 1.1);
    xlabel('time (ms)'); ylabel('Amplitude');
    grid on; grid minor;
    legend('ADSR', 'P2');
    title('ADSR curve and Synthetic P2')

figure(2)
ax21 = subplot(3,1,1);
    plot(totalRepeatedTimeAxis, repeatedEventA2, 'LineWidth', 1.1);
    xlabel('time (s)'); ylabel('Amplitude');
    xlim([0, totalRepeatedTimeAxis(end)]);
    grid on; grid minor;
    title('Synthetic A2')
ax22 = subplot(3,1,2);
    plot(totalRepeatedTimeAxis, repeatedEventP2, 'LineWidth', 1.1);
    xlabel('time (s)'); ylabel('Amplitude');
    xlim([0, totalRepeatedTimeAxis(end)]);
    grid on; grid minor;
    title('Synthetic P2')
ax23 = subplot(3,1,3);
    plot(totalRepeatedTimeAxis, syntheticS2, 'LineWidth', 1.1);
    xlabel('time (s)'); ylabel('Amplitude');
    xlim([0, totalRepeatedTimeAxis(end)]);
    grid on; grid minor;
    title('Synthetic S2 i.e. A2 + P2')
linkaxes([ax21, ax22, ax23], 'x');

figure(3);
ax31 = subplot(2,1,1);
    plot(totalRepeatedTimeAxis, syntheticS2,  'LineWidth', 1.1);
    xlabel('time (s)'); ylabel('Amplitude');
    xlim([0, totalRepeatedTimeAxis(end)]);
    grid on; grid minor;
    title('One complete S2 event')
ax32 = subplot(2,1,2);
    imagesc(totalRepeatedTimeAxis, f, scalogrm);
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    xlim([0, totalRepeatedTimeAxis(end)]);ylim([20, f(1)]);
    title('Microphone 1');
    set(gca,'YDir','normal');
    set(gca,'YScale','log');
linkaxes([ax31, ax32], 'x');

%% Helpers

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