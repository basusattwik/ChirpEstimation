clc
close all
clearvars

%% Create a chirp

fs = 100;
Td = 1;
Nc = 2; 

N = Td * fs;


% Parameters
f01 = 31.5;
f02 = 22.5;
t   =   (0:1/fs:Td-1/fs).';
ym  = exp(2*pi*1j * f01 * t) + exp(2*pi*1j * f02 * t) + 1 * complex(randn(N,1), rand(N,1));

% 1. Grid of parameters
f = 0:0.1:fs/2;

J = zeros(numel(f),1);
H = zeros(N, Nc);

for i = 1:numel(f)
    H(:,1) = exp(2*pi*1j * f(i) * t);
    for j = 1:numel(f)   
        H(:,2) = exp(2*pi*1j * f(j) * t);       
        
        P  = H * inv(H' * H) * H';
        Po = eye(N,N) - P;
        J(i, j) = real(ym' * Po * ym);
    end
end


%% Plots

close all

[X, Y] = meshgrid(f, f);
surf(X, Y, J,'EdgeColor','none');
xlabel('Frequency');
ylabel('Frequency');
zlabel('Objective Function');
