clc
close all
clearvars

%% Create a chirp

fs = 100;
Td = 1;
Nc = 2; 

N = Td * fs;
t = (0:1/fs:Td-1/fs).';
snr = 6;

% Parameters
p0_1 = 0;
f0_1 = 4.5;
p0_2 = 0;
f0_2 = 6.5;

fstart_1 = f0_1; 
fend_1   = 35.5;

fstart_2 = f0_2; 
fend_2   = 41.5;

m0_1   = (fend_1 - fstart_1) / Td;
m0_2   = (fend_2 - fstart_2) / Td;

phi_1  = [p0_1; f0_1; m0_1];
phi_2  = [p0_2; f0_2; m0_2];

xm = getChirp(fs, N, p0_1, f0_1, m0_1) + getChirp(fs, N, p0_2, f0_2, m0_2); 
ym = addGaussianNoise(xm, snr);

%% Cost function 

f = 0:0.8:12;
m = 20:0.8:45;

J = zeros(numel(f), numel(m), numel(f), numel(m));
H = zeros(N, Nc);

for find_1 = 1:numel(f)
    for mind_1 = 1:numel(m) 

        % Form basis matrix row 1
        H(:,1) = getChirp(fs, N, 0, f(find_1), m(mind_1));
        
        for find_2 = 1:numel(f)
            for mind_2 = 1:numel(m)

                % Form basis matrix row 2
                H(:,2) = getChirp(fs, N, 0, f(find_2), m(mind_2));
        
                % Compute cost function
                J(find_1, mind_1, find_2, mind_2) = - ym' * (H * pinv(H' * H) * H') * ym;
            end
        end
    end
end

J = real(J);

%% Plots

[F, M] = meshgrid(f, m);

f0_1_set = 6;
m0_1_set = 15;
f0_2_set = 9;
m0_2_set = 20;

figure('windowstyle','docked')
tiledlayout flow
nexttile
    plot(f, squeeze(J(:, m0_1_set, f0_2_set, m0_2_set)));
    xlabel('Frequency (Hz)');
    ylabel('Objective Function');
    title('Fixed chirp rate 1, frequency 2 and chirp rate 2');
    grid on; grid minor;
nexttile
    plot(m, squeeze(J(f0_1_set, :, f0_2_set, m0_2_set)));
    xlabel('Chirp Rate (Hz/z)');
    ylabel('Objective Function');
    title('Fixed frequency 1, frequency 2 and chirp rate 2');
    grid on; grid minor;
nexttile
    plot(f, squeeze(J(f0_1_set, m0_1_set, :, m0_2_set)));
    xlabel('Frequency (Hz)');
    ylabel('Objective Function');
    title('Fixed frequency 1, frequency 2 and chirp rate 2');
    grid on; grid minor;
nexttile
    plot(m, squeeze(J(f0_1_set, m0_1_set, f0_2_set, :)));
    xlabel('Chirp Rate (Hz/z)');
    ylabel('Objective Function');
    title('Fixed frequency 1, chirp rate 1, frequency 2');
    grid on; grid minor;

figure('windowstyle','docked')
tiledlayout flow
nexttile
    surf(F, M, squeeze(J(:, :, f0_2_set, m0_2_set)).', 'EdgeColor', 'interp', 'FaceColor', 'interp');
    xlabel('Frequency 1 (Hz)');
    ylabel('Chirp Rate 1 (Hz/s)');
    zlabel('Objective Function');
    title('Fixed frequency 2 and chirp rate 2');
    grid on; grid minor;
nexttile
    surf(F, M, squeeze(J(f0_1_set, m0_1_set, :, :)).', 'EdgeColor', 'interp', 'FaceColor', 'interp');
    xlabel('Frequency 2 (Hz)');
    ylabel('Chirp Rate 2 (Hz/s)');
    zlabel('Objective Function');
    title('Fixed frequency 1 and chirp rate 1');
    grid on; grid minor;
nexttile
    surf(m, m, squeeze(J(f0_1_set, :, f0_2_set, :)).', 'EdgeColor', 'interp', 'FaceColor', 'interp');
    xlabel('Chirp Rate 1 (Hz/s)');
    ylabel('Chirp Rate 2 (Hz/s)');
    zlabel('Objective Function');
    title('Fixed frequency 1 and frequency 2');
    grid on; grid minor;
nexttile
    surf(f, f, squeeze(J(:, m0_1_set, :, m0_2_set)).', 'EdgeColor', 'interp', 'FaceColor', 'interp');
    xlabel('Frequency 1 (Hz/s)');
    ylabel('Frequency 2 (Hz/s)');
    zlabel('Objective Function');
    title('Fixed chirp rate 1 and chirp rate 2');
    grid on; grid minor;
nexttile
    surf(m, f, squeeze(J(f0_1_set, :, :, m0_2_set)).', 'EdgeColor', 'interp', 'FaceColor', 'interp');
    xlabel('Frequency 1 (Hz/s)');
    ylabel('Frequency 2 (Hz/s)');
    zlabel('Objective Function');
    title('Fixed frequency 1 and chirp rate 2');
    grid on; grid minor;
nexttile
    surf(f, m, squeeze(J(:, m0_1_set, f0_2_set, :)).', 'EdgeColor', 'interp', 'FaceColor', 'interp');
    xlabel('Frequency 1 (Hz/s)');
    ylabel('Frequency 2 (Hz/s)');
    zlabel('Objective Function');
    title('Fixed chirp rate 1 and frequency 2');
    grid on; grid minor;

%% Helpers


function y = getChirp(fs, N, p0, f0, m0)

phi = [p0; f0; m0];
P = size(phi,1);
y = zeros(N,1);
    for nind = 1:N % -- loop over number of samples   
       
        % Get the exponential polynomial phase sinusoid
        npvec = ((nind-1) / fs).^(0:P-1).'; % vectors of powers of n/fs
        y(nind,1) = exp(2*pi*1j .* (phi.' * npvec));
    
    end % -- end loop over number of chirps
end

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
