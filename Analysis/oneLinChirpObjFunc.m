clc
close all
clearvars

% Create a chirp

fs = 200;
Td = 1;
Nc = 2; 

N = Td * fs;
snr = 3;

% Parameters
p0 = 0;
f0 = 35;

fstart = f0; 
fend   = 90;
m0  = (fend - fstart) / Td;
t   = (0:1/fs:Td-1/fs).';

phi = [p0; f0; m0];

xm = getChirp(fs, N, p0, f0, m0);
ym = addGaussianNoise(xm, snr);

plot(real(ym));

%% 1. Grid of parameters
f = 10:0.5:40;
m = 40:0.5:100;

J = zeros(numel(f),numel(m));
H = zeros(N, Nc);

for find = 1:numel(f)
    for mind = 1:numel(m)

        H  = getChirp(fs, N, 0, f(find), m(mind));
        P  = H * (H' * H)^(-1) * H';
    
        J(find, mind) = norm(ym)^2 - real(ym' * P * ym);
    end
end

%% Blur Objective function

J_blur = imgaussfilt(J, 20);

%% Plots

T = 1000;
% close all

[X, Y] = meshgrid(m, f);

figure('windowstyle','docked')
tiledlayout flow
ax1 = nexttile;
    surf(Y, X, (1/T .* J),'EdgeColor','none');
    xlabel('Frequency (Hz)');
    ylabel('Chirp Rate (Hz/s)');
    zlabel('Objective Function');

ax2 = nexttile;
    surf(Y, X, exp(-1/T .* J),'EdgeColor','none');
    xlabel('Frequency (Hz)');
    ylabel('Chirp Rate (Hz/s)');
    zlabel('Objective Function');

ax3 = nexttile;
    surf(Y, X, (1/T .* J_blur),'EdgeColor','none');
    xlabel('Frequency (Hz)');
    ylabel('Chirp Rate (Hz/s)');
    zlabel('Objective Function');

ax4 = nexttile;
    surf(Y, X, exp(-1/T .* J_blur),'EdgeColor','none');
    xlabel('Frequency (Hz)');
    ylabel('Chirp Rate (Hz/s)');
    zlabel('Objective Function');

Link = linkprop([ax1, ax2, ax3, ax4],{'CameraUpVector', 'CameraPosition'});
setappdata(gcf, 'StoreTheLink', Link);
%% Helpers

function y = getChirp(fs, N, p0, f0, m0)

    phi = [p0; f0; m0];
    P   = size(phi,1);
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
        wm = randn(size(xm)) + 1j * randn(size(xm)); %  random noise signal
    else
        wm = randn(size(xm)); %  random noise signal
    end
    
    xpow = sum(abs(xm.^2)) ;
    wpow = sum(abs(wm.^2));
    
    pn  = xpow / (10^(snr/10)); % required SNR
    wsc = (wm ./ sqrt(wpow)) .* sqrt(pn); % normalizing the noisy signal to get unity power then multiplying it by the new power to achieve the required SNR
    
    ym = xm + wsc; % add noise to signal
    
    % wnewpow = sum(abs(wscaled.^2));
    % snrtrue = 10*log10(ypow / wnewpow);

end

