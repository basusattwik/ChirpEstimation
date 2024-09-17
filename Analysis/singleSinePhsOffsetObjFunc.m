clc
% close all
clearvars

%% Create a chirp

fs = 20000;
Td = 1;
Nc = 1; 

N = Td * fs;


% Parameters
f0 = 7500;
p0 = deg2rad(60);

df = fs/N;
tind = (0:1/fs:Td-1/fs).';

y = exp(1j * (2*pi*f0*tind + p0)) + 0.0001 * complex(randn(N,1), randn(N,1));
% y = sin(2*pi*f0*tind + p0) + 0.1 * randn(N,1);

% 1. Grid of parameters

f = 7000:df:8000;
p = 0:pi/3:pi;

J = zeros(numel(f),numel(p));
I = eye(N,N);
for find = 1:numel(f)

    for pind = 1:numel(p)
         
        H  = exp(1j * (2*pi*f(find)*tind + p(pind)));
        % H  = sin(2*pi*f(find)*tind + p(pind));
        P  = H * (H' * H)^(-1) * H';
        Po = I - P;

        J(find, pind) = exp(-real(y' * Po * y));
    end
end

%% Gaussian smoothing J

%% Get distributions

[F,P] = meshgrid(f,p);

figure('windowstyle','docked')
surf(F, P, J.','EdgeColor','none');
xlabel('Frequency (Hz)');
ylabel('Phase offset (rad)');
grid on; grid minor;

