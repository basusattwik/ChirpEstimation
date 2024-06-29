close all
clearvars
clc

% This script performs Pincus Maximization on a 1D function with two local
% maxima and one unique global maxima. A direct closed form integral is
% evaluated in the one method and is followed by an importance sampling
% approach for optimization. 

%% General parameters

% Random seed
seed = 12345;
rng(seed);

% Domain of integration
x = -10:1e-4:10;

% Lambdas for Pincus
lambda = [0.05, 0.5, 1, 5, 10];

% Number of draws for Smirnov and Importance Sampling
ndraws = 100000;

% Mean and variance of Gaussian proposal distribution
mu    = 0.5;
sigma = 1;
c = 1;
d = 0;

%% Setup the function handles

% Function to be maximized
f = @(x) exp(-0.3*x - 0.5*x.^2) .* (1 + (sin(3*x).^2));

% Various relevant functions for Pincus
ef  = @(x, lambda) exp(lambda * f(x)); % Denominator
xef = @(x, lambda) x .* ef(x, lambda); % Numerator
g   = @(x, lambda) ef(x, lambda) / trapz(x, ef(x, lambda)); % pseudo density function for importance sampling

% Closed form Pincus optimizer
xopt = @(x, lambda) trapz(x, xef(x, lambda)) / trapz(x, ef(x, lambda)); 

% CDF for Smirnov Transform
G = @(x, lambda) cumtrapz(x, g(x, lambda)); % pseudo cdf of g(x)

% Choose proposal distribution q(x) as a Gaussian for Importance Sampling
q = @(x, mu, sigma, c, d) c * 1/(sqrt(2*pi)*sigma) * exp(-(x-mu).^2/(2*sigma^2)) + d; % Proposal distribution

ncases = numel(lambda);

%% Direct Pincus Maximization 

% Run loop over all lambdas and see which ones give the right answer
xo_direct = zeros(ncases,1);
for i = 1:ncases
    xo_direct(i,1) = xopt(x, lambda(i)); % Pincus closed form for optimizer
end

%% Optimizing via inverse CDF sampling (Smirnov Transform)

u = rand(ndraws,1);

lenx = length(x);
xo_smirnov = zeros(ncases, 1);
for i = 1:ncases
    xi = interp1(G(x, lambda(i)) + 1e-5*rand(1,lenx), x, u, 'spline');
    xo_smirnov(i,1) = mean(xi);
end

%% Maximization of the closed form Pincus Integral through Importance Sampling

xo_impSamp = zeros(ncases, 1);
for i = 1:ncases
    xn = c * normrnd(mu, sigma, ndraws, 1) + d;
    wt = g(xn, lambda(i)) ./ q(xn, mu, sigma, c, d);
    xo_impSamp(i,1) = mean((wt ./ sum(wt)) .* xn);
end

% impWeight = 

%% Plots

xplot = -5:1e-4:5;

figure(1)
t1 = tiledlayout('flow');
nexttile
plot(xplot, f(xplot), 'LineWidth', 1.1);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$f(x)$', 'Interpreter','latex');
title('$e^{-0.3x - 0.5x^2}(1 + sin^2(3x))$', 'Interpreter', 'latex', 'FontSize', 12);
grid on; grid minor;
for i = 1:ncases
    nexttile
    plot(xplot, ef(xplot,lambda(i)), 'LineWidth', 1.1);
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$f(x)$', 'Interpreter','latex');
    title(['$e^{\lambda f(x)}, \lambda =\ $', num2str(lambda(i))], 'Interpreter', 'latex', 'FontSize', 12);
    grid on; grid minor;
end
title(t1, 'Effect of changing $\lambda$ on $e^{\lambda f(x)}$', 'Interpreter', 'latex', 'FontSize', 14);

figure(2)
t1 = tiledlayout('flow');
nexttile
plot(xplot, f(xplot), 'LineWidth', 1.1);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$f(x)$', 'Interpreter','latex');
title('$e^{-0.3x - 0.5x^2}(1 + sin^2(3x))$', 'Interpreter', 'latex', 'FontSize', 12);
grid on; grid minor;
for i = 1:ncases
    nexttile
    plot(xplot, xef(xplot,lambda(i)), 'LineWidth', 1.1);
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$f(x)$', 'Interpreter','latex');
    title(['$xe^{\lambda f(x)}, \lambda =\ $', num2str(lambda(i))], 'Interpreter', 'latex', 'FontSize', 12);
    grid on; grid minor;
end
title(t1, 'Effect of changing $\lambda$ on $xe^{\lambda f(x)}$', 'Interpreter', 'latex', 'FontSize', 14);

figure(3)
t1 = tiledlayout('flow');
for i = 1:ncases
    nexttile
    plot(G(xplot, lambda(i)), xplot, 'LineWidth', 1.1);
    xlabel('$u \sim U[0,1]$', 'Interpreter', 'latex');
    ylabel('$G^{-1}(u)$', 'Interpreter','latex');
    title(['$\lambda =\ $', num2str(lambda(i))], 'Interpreter', 'latex', 'FontSize', 12);
    grid on; grid minor;
end
title(t1, 'Effect of changing $\lambda$ on CDF for Smirnov', 'Interpreter', 'latex', 'FontSize', 14);

figure(4)
t1 = tiledlayout('flow');
for i = 1:ncases
    nexttile
    plot(xplot, g(xplot,lambda(i)), 'LineWidth', 1.1, 'DisplayName', '$g(x)$'); hold on;
    plot(xplot, q(xplot, mu, sigma, c, d), 'LineWidth', 1.1, 'DisplayName', '$q(x)$');
    plot(x, g(x,lambda(i)) ./  q(x, mu, sigma, c, d), 'LineWidth', 1.1, 'DisplayName', '$g(x) / q(x)$');
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$g(x)$', 'Interpreter','latex');
    title(['$g(x), \lambda =\ $', num2str(lambda(i))], 'Interpreter', 'latex', 'FontSize', 12);
    grid on; grid minor;
    hl = legend('show');
    set(hl, 'Interpreter','latex');
end
title(t1, '$g(x)$ vs $q(x)$', 'Interpreter', 'latex', 'FontSize', 14);

figure(5)
subplot(3,1,1)
    plot(xplot, f(xplot), 'LineWidth', 1.1, 'DisplayName', '$f(x)$'); hold on;
    for i = 1:ncases
        stem(xo_direct(i), f(xo_direct(i)), 'filled', 'DisplayName', ['$\lambda =\ $', num2str(lambda(i))]);
    end
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$f(x)$', 'Interpreter','latex');
    title('Direct Integration of Closed Form');
    grid on; grid minor;
    hl = legend('show');
    set(hl, 'Interpreter','latex')
subplot(3,1,2)
    plot(xplot, f(xplot), 'LineWidth', 1.1, 'DisplayName', '$f(x)$'); hold on;
    for i = 1:ncases
        stem(xo_smirnov(i), f(xo_smirnov(i)), 'filled', 'DisplayName', ['$\lambda =\ $', num2str(lambda(i))]);
    end
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$f(x)$', 'Interpreter','latex');
    title('Smirnov Transform');
    grid on; grid minor;
    hl = legend('show');
    set(hl, 'Interpreter','latex')
subplot(3,1,3)
    plot(xplot, f(xplot), 'LineWidth', 1.1, 'DisplayName', '$f(x)$'); hold on;
    for i = 1:ncases
        stem(xo_impSamp(i), f(xo_impSamp(i)), 'filled', 'DisplayName', ['$\lambda =\ $', num2str(lambda(i))]);
    end
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$f(x)$', 'Interpreter','latex');
    title('Importance Sampling');
    grid on; grid minor;
    hl = legend('show');
    set(hl, 'Interpreter','latex')
sgtitle('Effect of changing \lambda on the optimizer');
