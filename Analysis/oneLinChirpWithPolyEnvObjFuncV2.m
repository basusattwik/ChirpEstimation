clc
close all
clearvars

%% Create a chirp

fs = 500;
Td = 1;
Nc = 1; 

% c = [0.374455411205074	-0.920152301725562	0.536018714179969].';
c = [1 0 0].';

N = Td * fs;
t = (0:1/fs:Td-1/fs).';
snr = 60;

% Parameters
p0 = 0;
f0 = 60;

% fstart = f0; 
% fend = 33;
m0 = 50; %(fend-fstart)/Td;

phi = [p0; f0; m0];

xm = getEnvFunc(fs, N, c) .* getChirp(fs, N, p0, f0, m0);
ym = addGaussianNoise(xm, snr);

nvec = (0:N-1).'/fs;

%% Cost function 

f = 20:0.05:100;
m = 20:0.05:100;
J = zeros(numel(f), numel(m));
H = zeros(N, 3);

for find = 1:numel(f)
    for mind = 1:numel(m) 

        % Get chirp
        chr = getChirp(fs, N, 0, f(find), m(mind));

        % Form basis matrix col 2
        H(:,1) = chr;
        H(:,2) = chr .* nvec;
        H(:,3) = chr .* nvec.^2;

        P  = H * ((H' * H) \ H');
        Po = eye(N,N) - P;

        % Compute cost function
        J(find, mind) = real(ym' * (Po * ym));
    end
end

save("JTensor_oneLinChirpNoEvn.mat", "J");

J = real(J);

%% Plots
load("JTensor_oneLinChirpNoEvn.mat", "J");

[F,M] = meshgrid(f,m);

% figure('windowstyle','docked')
%     plot(t, real(ym)); hold on;
%     plot(t, imag(ym));
%     xlabel('Time (s)');
%     ylabel('Amplitude');
%     title('Chirp Signal vs Time');
%     grid on; grid minor;
figure('windowstyle','docked')
    surf(F, M, J.', 'EdgeColor', 'interp', 'FaceColor', 'interp');
    xlabel('Frequency (Hz)');
    ylabel('Chirp Rate (Hz/s)');
    zlabel('Objective Function');
    title('Fixed envelope parameters');
    grid on; grid minor;
figure('windowstyle','docked')
    contour(F, M, J.', 'LineWidth', 1.2);
    xlabel('Frequency (Hz)');
    ylabel('Chirp Rate (Hz/s)');
    zlabel('Objective Function');
    title('Fixed envelope parameters');
    grid on; grid minor;

%% Helpers

function y = getChirp(fs, N, p0, f0, m0)

    phi = [p0; f0; m0];
    P  = size(phi,1);
    y  = zeros(N,1);
    
    for nind = 1:N       
        % Get the exponential polynomial phase sinusoid
        npvec     = ((nind-1) / fs).^(0:P-1).'; 
        y(nind,1) = exp(2*pi*1j .* (phi.' * npvec));
    
    end 
end

function A = getEnvFunc(fs, N, c)

    K  = size(c,1);
    A  = zeros(N,1);
    
    for nind = 1:N       
        % Get the exponential polynomial phase sinusoid
        nkvec     = ((nind-1) / fs).^(0:K-1).'; 
        A(nind,1) = c.' * nkvec;
    
    end 
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
