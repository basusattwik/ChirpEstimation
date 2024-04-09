close all
clearvars
clc

%% Signal model and Cost

% Actual signal
fs = 200; % Hz
Td = 1; % 1 sec
t  = 0:1/fs:Td-1/fs;
f0 = 40;
m0 = 10;
x  = sin(2*pi * (f0*t + 0.5*m0*t.^2)).';

% Add noise
lx  = length(x);
var = 0.01;
y   = x + var*randn(lx,1);

% Cost function
xhat = @(f,m,t) sin(2*pi * (f*t + 0.5*m*t.^2)).'; % Model

fx  = 0:0.01:fs/2;
mx  = 0:0.01:20;
Jfm  = zeros(numel(fx), numel(mx));
feJfm = zeros(numel(fx), numel(mx));
meJfm = zeros(numel(fx), numel(mx));
lambda = 10;
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

lambda = 100;

% Various relevant functions for Pincus
% ej  = @(f, m, lambda) exp(lambda * Jfm);     % Denominator
% fej = @(f, m, lambda) f .* ej(f, m, lambda); % Numerator
% mej = @(f, m, lambda) m .* ej(f, m, lambda); % Numerator
den  = trapz(fx, trapz(mx, eJfm, 2));

fopt = trapz(fx, trapz(mx, feJfm, 2)) / den;
mopt = trapz(fx, trapz(mx, meJfm, 2)) / den;

%% Plots
close all

figure(1)
plot(t, x); hold on;
plot(t, y);
% plot(t, xhat(fopt,t));
xlabel('Time (s)');
ylabel('Amplitude');
grid on; grid minor;
% legend('Actual', 'Measured', 'Recons.');
title('Comparison of original vs reconstructed signals');

% lambdaArr = [0.5, 1, 2, 10, 20];
% 
% figure(2)
% t1 = tiledlayout('flow');
% nexttile
% plot(fx, Jfm, 'LineWidth', 1.1);
% xlabel('$f$', 'Interpreter', 'latex');
% ylabel('$J(f)$', 'Interpreter','latex');
% title('$J(f) = \| y - sin(2\pi ft) \|^2$', 'Interpreter', 'latex', 'FontSize', 12);
% grid on; grid minor;
% for i = 1:numel(lambdaArr)
%     nexttile
%     plot(fx, ej(fx, lambdaArr(i)), 'LineWidth', 1.1);
%     xlabel('$x$', 'Interpreter', 'latex');
%     ylabel('$e^{\lambda J(f)}$', 'Interpreter','latex');
%     title(['$\lambda =\ $', num2str(lambdaArr(i))], 'Interpreter', 'latex', 'FontSize', 12);
%     grid on; grid minor;
% end
% title(t1, 'Effect of changing $\lambda$ on $e^{\lambda J(f)}$', 'Interpreter', 'latex', 'FontSize', 14);
% 
figure(3)
surf(mx, fx, Jfm, 'EdgeColor','none'); hold on;
% stem(fopt, 0, 'filled', 'v');
xlabel('Rate (Hz/s)');
ylabel('Frequency (Hz)');
zlabel('Cost');
title('Negative cost function at different frequencies');
grid on; grid minor;