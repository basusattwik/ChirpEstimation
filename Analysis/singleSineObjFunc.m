clc
% close all
clearvars

% This is not correct. Noise should most likely not be added to the
% measured signal. 

% It should be added during the LMC updates to the parameters
% This has the effect of Gaussian smoothing the gradients
%% Create a chirp

fs = 200;
Td = 1;
Nc = 1; 

N = Td * fs;


% Parameters
f0 = 64.5;
tind  = (0:1/fs:Td-1/fs).';

y1 = exp(2*pi*1j * f0*tind) + 3.99 * complex(randn(N,1), randn(N,1));
y2 = exp(2*pi*1j * f0*tind) + 3.99 * complex(randn(N,1), randn(N,1)) + 0*complex(1.5 * randn(N,1), 1.5 * randn(N, 1));
% m = -1.1;
% M = 0;
% y2 = exp(2*pi*1j * f0*tind) + 0.01 * complex(randn(N,1), randn(N,1)) + complex(unifrnd(m, M, N, 1), unifrnd(m, M, N, 1));

% %% Add Gaussian Noise to signal
% 
% if ~isreal(obj.ym)
%     obj.wm = randn(size(obj.ym)) + 1j * randn(size(obj.ym)); %  random noise signal
% else
%     obj.wm = randn(size(obj.ym)); %  random noise signal
% end
% 
% ypow = sum(abs(obj.ym.^2)) ;
% wpow = sum(abs(obj.wm.^2));
% 
% pn  = ypow / (10^(obj.snr/10)); % required SNR
% wsc = (obj.wm ./ sqrt(wpow)) .* sqrt(pn); % normalizing the noisy signal to get unity power then multiplying it by the new power to achieve the required SNR
% 
% obj.ym = obj.ym + wsc; % add noise to signal


% 1. Grid of parameters

f = 0:0.005:fs/2;

J1  = zeros(numel(f),1);
dJ1 = zeros(numel(f),1);
J2  = zeros(numel(f),1);
dJ2 = zeros(numel(f),1);
for find = 1:numel(f)

    H1  = exp(2*pi*1j * f(find) * tind);
    P1  = H1 * (H1' * H1)^(-1) * H1';
    Po1 = eye(N,N) - P1;
    dH1 = 2*pi*1j .* tind .* H1;
    J1(find)  = real(y1' * Po1 * y1);
    dJ1(find) = -2*real(y1' * Po1 * dH1 * (H1' * H1)^(-1) * H1' * y1);

    H2  = exp(2*pi*1j * f(find) * tind);
    P2  = H2 * (H2' * H2)^(-1) * H2';
    Po2 = eye(N,N) - P2;
    dH2 = 2*pi*1j .* tind .* H2;
    J2(find)  = real(y2' * Po2 * y2);
    dJ2(find) = -2*real(y2' * Po2 * dH2 * (H2' * H2)^(-1) * H2' * y2);
end

%% Gaussian smoothing J

 w = gausswin(512, 3);

Jcw = filter(w, 1, J1);
dJcw = filter(w, 1, dJ1);

figure(1)
plot(f, Jcw); 
hold on;
plot(f, J1)

figure(2)
plot(f, dJcw); 
hold on;
plot(f, dJcw)

%% Get distributions

T = [1];%, 1, 0.1];
eP1    = zeros(numel(f), numel(T));
eP2     = zeros(numel(f), numel(T));

for tind = 1:numel(T)
    
    eP1(:, tind) = exp(-1/T(tind) .* J1);
    eP2(:, tind) = exp(-1/T(tind) .* J2);
end

% close all
figure('windowstyle','docked')
tiledlayout flow
nexttile
    plot(f, J1); hold on;
    plot(f, J2);
    grid on; 
nexttile
    plot(f, dJ1); hold on;
    plot(f, dJ2);
    grid on;
nexttile
    for tind = 1:numel(T)
        hold on;
        plot(f, eP1(:,tind), 'DisplayName', ['T = ', num2str(T(tind)), ' , \lambda = ', num2str(1/T(tind))]);
        plot(f, eP2(:,tind), 'DisplayName', ['T = ', num2str(T(tind)), ' , \lambda = ', num2str(1/T(tind))]);
    end
    legend('show');
    grid on; 