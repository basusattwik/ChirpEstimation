clc
close all
clearvars

%% Create a chirp

fs = 1000;
Td = 1;
Nc = 2; 

c = [0.374455411205074	-0.920152301725562	0.536018714179969].';

N = Td * fs;
t = (0:1/fs:Td-1/fs).';
snr = 6;

% Parameters
p0 = 0;
f0 = 5;

fstart = f0; 
fend   = 33;

beta  = 3;
gamma = 5;
delta = beta + gamma;

m0 = 28;%(fend-fstart)/Td;

phi = [p0; f0; m0];

xm = getEnvFunc(fs, N, c) .* getChirp(fs, N, p0, f0, m0);
ym = addGaussianNoise(xm, snr);

%% Cost function 

f   = 0:1:12;
m   = 20:1:45;
ac1 = 0:0.05:1;
ac2 = -1:0.05:0;
ac3 = 0:0.05:1;

J = zeros(numel(f), numel(m), numel(ac1), numel(ac2), numel(ac3));
H = zeros(N, 2*Nc);

for find = 1:numel(f)
    for mind = 1:numel(m) 

        % Get chirp
        chr = getChirp(fs, N, 0, f(find), m(mind));

        for c1ind = 1:numel(ac1)
            for c2ind = 1:numel(ac2) 
                for c3ind = 1:numel(ac3)   

                    c = [ac1(c1ind), ac2(c2ind), ac3(c3ind)].';
    
                    % Form basis matrix col 2
                    H(:,1) = chr .* getEnvFunc(fs, N, c);
                    H(:,2) = chr .* getEnvFunc(fs, N, c);
                    H(:,3) = chr .* getEnvFunc(fs, N, c);
            
                    % Compute cost function
                    J(find, mind, c1ind, c2ind, c3ind) = -ym' * (H * pinv(H' * H) * H') * ym;
                end
            end
        end
    end
end

save("JTensor.mat", "J");

J = real(J);
%% Plots

[F,M] = meshgrid(f,m);

f0_set  = 6;
m0_set  = 9;
c1_set  = 8;
c2_set  = 3;
c3_set  = 12;

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
    surf(F, M, squeeze(J(:, :, c1_set, c2_set, c3_set).'), 'EdgeColor', 'interp', 'FaceColor', 'interp');
    xlabel('Frequency (Hz)');
    ylabel('Chirp Rate (Hz/s)');
    zlabel('Objective Function');
    title('Fixed envelope parameters');
    grid on; grid minor;

nexttile
    surf(squeeze(J(f0_set, m0_set, c1_set, :, :)), 'EdgeColor', 'interp', 'FaceColor', 'interp');
    xlabel('Delta');
    ylabel('Beta');
    zlabel('Objective Function');
    title('Fixed frequency and chirp rate');
    grid on; grid minor;

% nexttile
%     surf(m, d, squeeze(J(f0_set, :, c1_set, :)).', 'EdgeColor', 'interp', 'FaceColor', 'interp');
%     xlabel('Frequency (Hz)');
%     ylabel('Gamma');
%     zlabel('Objective Function');
%     title('Fixed frequency and beta');
%     grid on; grid minor;
% 
% nexttile
%     surf(f, b, squeeze(J(:, m0_set, :, c2_set)).', 'EdgeColor', 'interp', 'FaceColor', 'interp');
%     xlabel('Frequency (Hz)');
%     ylabel('Beta');
%     zlabel('Objective Function');
%     title('Fixed chirp rate and gamma');
%     grid on; grid minor;
% 
% nexttile
%     contour(b, d, squeeze(J(f0_set, m0_set, :, :)).')
%     xlabel('Beta');
%     ylabel('Gamma');
%     title('Fixed frequency and chirp rate contour');
%     grid on; grid minor;
% 
% nexttile
%     contour(f, m, squeeze(J(:, :, c1_set, c2_set)).')
%     xlabel('Chirp Rate (Hz/s)');
%     ylabel('Frequency (Hz)');
%     title('Fixed beta and gamma contour');
%     grid on; grid minor;

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

% function pc = doPolyFit(x, Fx, K)
% 
% pc = polyfit(x, Fx, K);
% 
% end

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
