clc
close all
clearvars

%% Create a chirp

fs = 1000;
Td = 1;
Nc = 2; 

N = Td * fs;
t = (0:1/fs:Td-1/fs).';
snr = 6;

% Parameters
p0 = 0;
f0 = 4.8;

fstart = f0; 
fend   = 33;

beta  = 4;
gamma = 4.8;

m0 = (fend-fstart)/Td;

phi = [p0; f0; m0];

xm = getChirp(fs, N, p0, f0, m0, beta, gamma); 
ym = addGaussianNoise(xm, snr);

%% Cost function 

f = 0:0.6:12;
m = 20:0.6:45;
b = 0:0.8:10;
g = 0:0.8:10;

J = zeros(numel(f), numel(m), numel(b), numel(g));
H = zeros(N, Nc);

for find = 1:numel(f)
    for mind = 1:numel(m) 
        for bind = 1:numel(b)
            for gind = 1:numel(g)
                
                % Form basis matrix 
                H = getChirp(fs, N, 0, f(find), m(mind), b(bind), g(gind));
        
                % Compute cost function
                J(find, mind, bind, gind) = -ym' * (H * pinv(H' * H) * H') * ym;
            end
        end
    end
end

J = real(J);

% for bind = 1:numel(b)
%     for gind = 1:numel(g)

%% Plots

[F,M] = meshgrid(f,m);
[B,G] = meshgrid(b,g);

% [F, M, B, G] = ndgrid(f, m, b, g);

f0_set = 9;
m0_set = 15;
b_set  = 6;
g_set  = 7;

figure('windowstyle','docked')
tiledlayout flow
nexttile
    plot(t, real(ym)); hold on;
    plot(t, imag(ym));
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Chirp Signal vs Time');
    grid on; grid minor;

nexttile
    surf(F, M, squeeze(J(:, :, b_set, g_set).'), 'EdgeColor', 'interp', 'FaceColor', 'interp');
    xlabel('Frequency (Hz)');
    ylabel('Chirp Rate (Hz/s)');
    zlabel('Objective Function');
    title('Fixed envelope parameters');
    grid on; grid minor;

nexttile
    surf(g, b, squeeze(J(f0_set, m0_set, :, :)), 'EdgeColor', 'interp', 'FaceColor', 'interp');
    xlabel('Gamma');
    ylabel('Beta');
    zlabel('Objective Function');
    title('Fixed frequency and chirp rate');
    grid on; grid minor;

nexttile
    surf(m, g, squeeze(J(f0_set, :, b_set, :)).', 'EdgeColor', 'interp', 'FaceColor', 'interp');
    xlabel('Frequency (Hz)');
    ylabel('Gamma');
    zlabel('Objective Function');
    title('Fixed frequency and beta');
    grid on; grid minor;

nexttile
    surf(f, b, squeeze(J(:, m0_set, :, g_set)).', 'EdgeColor', 'interp', 'FaceColor', 'interp');
    xlabel('Frequency (Hz)');
    ylabel('Beta');
    zlabel('Objective Function');
    title('Fixed chirp rate and gamma');
    grid on; grid minor;

nexttile
    contour(b, g, squeeze(J(f0_set, m0_set, :, :)))
    xlabel('Beta');
    ylabel('Gamma');
    title('Fixed frequency and chirp rate contour');
    grid on; grid minor;

nexttile
    contour(m, f, squeeze(J(:, :, b_set, g_set)))
    xlabel('Chirp Rate (Hz/s)');
    ylabel('Frequency (Hz)');
    title('Fixed beta and gamma contour');
    grid on; grid minor;

%% Helpers

function yc = getChirp(fs, N, p0, f0, m0, beta, gamma)

phi = [p0; f0; m0];
P = size(phi,1);
A = zeros(N,1);
y = zeros(N,1);

    for nind = 1:N % -- loop over number of samples       
        % Envelope
        A(nind,1) = exp(-beta * nind / fs) * (1 - exp(-gamma * nind / fs)); 

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
