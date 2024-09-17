clc
close all
clearvars

% Objective function plots but without the 1 - exp(-gamma n/fs) term

%% Create a chirp

fs = 300;
Td = 1;
Nc = 2; 

N = Td * fs;
t = (0:1/fs:Td-1/fs).';
snr = 6;

% Parameters
p0 = 0;
f0 = 12;

fstart = f0; 
fend   = 40;

beta = 5.4;

m0 = 28.2; (fend-fstart)/Td;

phi = [p0; f0; m0];

xm = getChirp(fs, N, p0, f0, m0, beta); 
ym = addGaussianNoise(xm, snr);

%% Cost function 

f = 0:0.6:30;
m = 0:0.6:48;
b = 0:0.6:10;

J = zeros(numel(f), numel(m), numel(b));
H = zeros(N, Nc);

for find = 1:numel(f)
    for mind = 1:numel(m) 
        for bind = 1:numel(b)
            % Form basis matrix 
            H = getChirp(fs, N, 0, f(find), m(mind), b(bind));
    
            % Compute cost function
            J(find, mind, bind) = -ym' * (H * pinv(H' * H) * H') * ym;
        end
    end
end

J = real(J);

%% Plots

[F,M] = meshgrid(f,m);
% [B,G] = meshgrid(b,g);

f0_set = 21;
m0_set = 48;
b_set  = 10;

figure('windowstyle','docked')
tiledlayout flow

nexttile
    surf(f, m, squeeze(J(:, :, b_set)).', 'EdgeColor', 'interp', 'FaceColor', 'interp');
    xlabel('Frequency (Hz)');
    ylabel('Chirp Rate (Hz/s)');
    zlabel('Objective Function');
    title('Fixed envelope parameter');
    grid on; grid minor;

nexttile
    surf(m, b, squeeze(J(f0_set, :, :)).', 'EdgeColor', 'interp', 'FaceColor', 'interp');
    xlabel('Chirp Rate(Hz/s)');
    ylabel('Beta');
    zlabel('Objective Function');
    title('Fixed frequency');
    grid on; grid minor;

nexttile
    surf(f, b, squeeze(J(:, m0_set, :)).', 'EdgeColor', 'interp', 'FaceColor', 'interp');
    xlabel('Frequency (Hz)');
    ylabel('Beta');
    zlabel('Objective Function');
    title('Fixed chirp rate parameter');
    grid on; grid minor;

%% Helpers

function yc = getChirp(fs, N, p0, f0, m0, beta)

phi = [p0; f0; m0];
P = size(phi,1);
A = zeros(N,1);
y = zeros(N,1);

    for nind = 1:N % -- loop over number of samples       
        % Envelope
        A(nind,1) = exp(-beta * nind / fs); 

        % Get the exponential polynomial phase sinusoid
        npvec     = ((nind-1) / fs).^(0:P-1).'; 
        y(nind,1) = exp(2*pi*1j .* (phi.' * npvec));
    
    end % -- end loop over number of chirps
    yc = A .* y;
end

function ym = addGaussianNoise(xm, snr)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    
    if ~isreal(xm)
        wm = randn(size(xm)) + 1j * randn(size(xm)); % random noise signal
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
