clc
% close all
clearvars

% Trying to see what happens if we vary chirp rate and frequency for any
% arbitrary phase offset. seems like there is a global minimum. 

%% Create a chirp

fs = 500;
Td = 0.8;
Nc = 1; 

N = Td * fs;


% Parameters
f0 = 128;
m0 = 40;
p0 = deg2rad(60);
tind = (0:1/fs:Td-1/fs).';

y = exp(1j * (2*pi*f0*tind + 2*pi * m0 * tind.^2 + p0)) + 0.1 * complex(randn(N,1), randn(N,1));

% 1. Grid of parameters

f = 0:0.1:fs/2;
p = 0:pi/30:pi;
m = 20:0.1:50;

J = zeros(numel(f),numel(p));
I = eye(N,N);
for find = 1:numel(f)

    for mind = 1:numel(m)
         
        H  = exp(1j * (2*pi*f(find)*tind + 2*pi*m(mind)*(tind).^2));
        P  = H * pinv(H' * H) * H';
        Po = I - P;

        J(find, mind) = real(y' * Po* y);
    end
end

%% Gaussian smoothing J

%% Get distributions

[F,M] = meshgrid(f,m);

figure('windowstyle','docked')
surf(F, M, J.','EdgeColor','none');
xlabel('Frequency (Hz)');
ylabel('Chirp Rate (Hz/z)');
grid on; grid minor;