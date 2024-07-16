close all
clearvars
clc

%% Signal model and Cost

% Actual signal
fs = 200; % Hz
Td = 1; % 1 sec
t  = 0:1/fs:Td-1/fs;
f0 = 30;                                                                                                                                                                                                                                                                                                                                                                                                                                            
m0 = 10;
x  = sin(2*pi * (f0*t + 0.5*m0*t.^2)).';

% Add noise
lx  = length(x);
var = 0.001;
y   = x + var*randn(lx,1);

% Cost function
xhat = @(f,m,t) sin(2*pi * (f*t + 0.5*m*t.^2)).'; % Model

fx  = 0:0.05:fs/2;
mx  = 0:0.05:20;
Jfm   = zeros(numel(fx), numel(mx));
feJfm = zeros(numel(fx), numel(mx));
meJfm = zeros(numel(fx), numel(mx));
lambda = 3;
for find = 1:numel(fx)
    for mind = 1:numel(mx)
        error = y - xhat(fx(find),mx(mind),t);
        Jfm(find, mind)  = -(1/lx) * (error.' * error); % looking at negative for maximization

        feJfm(find, mind) = fx(find) * exp(lambda * Jfm(find, mind));
        meJfm(find, mind) = mx(mind) * exp(lambda * Jfm(find, mind));
        
    end
end

eJfm = exp(lambda * Jfm);

%% Pincus Maximization

lambda = 5;

den  = trapz(fx, trapz(mx, eJfm, 2));
fopt = trapz(fx, trapz(mx, feJfm, 2)) / den;
mopt = trapz(fx, trapz(mx, meJfm, 2)) / den;

%% Plots
close all

figure(1)
plot(t, x); hold on;
plot(t, y);
plot(t, xhat(fopt, mopt, t));
xlabel('Time (s)');
ylabel('Amplitude');
grid on; grid minor;
legend('Actual', 'Measured', 'Recons.');
title('Comparison of original vs reconstructed signals');

figure(2)
surf(mx, fx, Jfm, 'EdgeColor','none'); hold on;
xlabel('Rate (Hz/s)');
ylabel('Frequency (Hz)');
zlabel('Cost');
title('Negative cost function at different frequencies');
grid on; grid minor;

figure(3)
surf(mx, fx, eJfm, 'EdgeColor','none'); hold on;
xlabel('Rate (Hz/s)');
ylabel('Frequency (Hz)');
zlabel('Cost');
title(['$e^{\lambda J(f,m)}$ for $\lambda$ = ', num2str(lambda)], 'Interpreter','latex');
grid on; grid minor;