close all
clearvars
clc

%% Signal model and Cost

% Actual signal
fs = 200; % Hz
Td = 1; % 1 sec
t  = 0:1/fs:Td-1/fs;
f0 = 40;
x  = sin(2*pi*f0*t).';

% Add noise
lx  = length(x);
var = 0.01;
y   = x + var*randn(lx,1);

% Cost function
xhat = @(f,t) sin(2*pi*f*t).'; % Model

fx = 0:0.001:fs/2;
Jf = zeros(size(fx));
for find = 1:numel(fx)
    error    = y - xhat(fx(find),t);
    Jf(find) = -(1/lx) * (error.' * error); % looking at negative for maximization
end

%% Pincus Maximization

lambda = 100;

% Various relevant functions for Pincus
ej  = @(f, lambda) exp(lambda * Jf);   % Denominator
fej = @(f, lambda) f .* ej(f, lambda); % Numerator

% Closed form Pincus optimizer
fopt = trapz(fx, fej(fx, lambda)) / trapz(fx, ej(fx, lambda));

%% Plots

figure(1)
plot(t, x); hold on;
plot(t, y);
plot(t, xhat(fopt,t));
xlabel('Time (s)');
ylabel('Amplitude');
grid on; grid minor;
legend('Actual', 'Measured', 'Recons.');
title('Comparison of original vs reconstructed signals');

lambdaArr = [0.5, 1, 2, 10, 20];

figure(2)
t1 = tiledlayout('flow');
nexttile
plot(fx, Jf, 'LineWidth', 1.1);
xlabel('$f$', 'Interpreter', 'latex');
ylabel('$J(f)$', 'Interpreter','latex');
title('$J(f) = \| y - sin(2\pi ft) \|^2$', 'Interpreter', 'latex', 'FontSize', 12);
grid on; grid minor;
for i = 1:numel(lambdaArr)
    nexttile
    plot(fx, ej(fx, lambdaArr(i)), 'LineWidth', 1.1);
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$e^{\lambda J(f)}$', 'Interpreter','latex');
    title(['$\lambda =\ $', num2str(lambdaArr(i))], 'Interpreter', 'latex', 'FontSize', 12);
    grid on; grid minor;
end
title(t1, 'Effect of changing $\lambda$ on $e^{\lambda J(f)}$', 'Interpreter', 'latex', 'FontSize', 14);

figure(3)
plot(fx, Jf); hold on;
stem(fopt, 0, 'filled', 'v');
xlabel('Frequency (Hz)');
ylabel('Cost');
title('Negative cost function at different frequencies');
grid on; grid minor;