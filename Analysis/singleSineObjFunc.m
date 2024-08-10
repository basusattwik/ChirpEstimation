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

y = exp(2*pi*1j * f0*tind) + 0.01 * complex(randn(N,1), randn(N,1));

% 1. Grid of parameters

f = 0:0.005:fs/2;

J  = zeros(numel(f),1);
dJ = zeros(numel(f),1);
for find = 1:numel(f)

    H  = exp(2*pi*1j * f(find) * tind);
    P  = H * (H' * H)^(-1) * H';
    Po = eye(N,N) - P;
    dH = 2*pi*1j .* tind .* H;
    J(find)  = real(y' * Po * y);
    dJ(find) = -2*real(y' * Po * dH * (H' * H)^(-1) * H' * y);
end

%% Gaussian smoothing J

 w = gausswin(512, 3);

Jcw = filter(w, 1, J);
dJcw = filter(w, 1, dJ);

figure(1)
plot(f, Jcw); 
hold on;
plot(f, J)

figure(2)
plot(f, dJcw); 
hold on;
plot(f, dJcw)

%% Get distributions

T = [1];%, 1, 0.1];
eP1    = zeros(numel(f), numel(T));

for tind = 1:numel(T) 
    eP1(:, tind) = exp(-1/T(tind) .* J);
end

% close all
figure('windowstyle','docked')
tiledlayout flow
nexttile
    plot(f, J);

    grid on; 
nexttile
    plot(f, dJ); hold on;
    grid on;
nexttile
    for tind = 1:numel(T)
        hold on;
        plot(f, eP1(:,tind), 'DisplayName', ['T = ', num2str(T(tind)), ' , \lambda = ', num2str(1/T(tind))]);
    end
    legend('show');
    grid on; 