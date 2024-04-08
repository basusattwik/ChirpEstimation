close all
clearvars
clc

% This script performs Pincus Maximization on a 1D function with two local
% maxima and one unique global maxima. A direct closed form integral is
% evaluated in the one method and is followed by an importance sampling
% approach for optimization. 

%% Setup the function handles

% Domain of integration
x    = -10:1e-4:10;
lenx = length(x);

% Various relevant functions
f    = @(x) exp(-0.3*x - 0.5*x.^2) .* (1 + (sin(3*x).^2)); % Function to be maximized
ef   = @(x, lambda) exp(lambda * f(x)); % Denominator
xef  = @(x, lambda) x .* ef(x, lambda); % Numerator
g    = @(x, lambda) ef(x, lambda) / trapz(x, ef(x, lambda)); % pseudo density function for importance sampling
G    = @(x, lambda) cumtrapz(x, g(x, lambda)); % pseudo cdf of g(x)

% Closed form optimizer
xopt = @(x, lambda) trapz(x, xef(x, lambda)) / trapz(x, ef(x, lambda)); 

%% Direct Pincus Maximization 

lambda = [0.05, 0.5, 1, 10, 50];

% Run loop over all lambdas and see which ones give the right answer
ntries = numel(lambda);
xo_direct = zeros(ntries,1);
for i = 1:ntries
    xo_direct(i,1) = xopt(x, lambda(i)); % Pincus closed form for optimizer
end

%% Optimizing via inverse CDF sampling (Smirnov Transform)

ndraws = 100000;
u  = rand(ndraws,1);

xo_smirnov = zeros(ntries, 1);
for i = 1:ntries
    xi = interp1(G(x, lambda(i)) + 1e-5*rand(1,lenx), x, u, 'spline');
    xo_smirnov(i,1) = mean(xi);
end

%% Maximization of the closed form Pincus Integral through Importance Sampling

% Choose proposal distribution as a Gaussian 
q  = @(x, mu, sigma) c * 1/(sqrt(2*pi)*sigma) * exp(-(x-mu).^2/(2*sigma^2));

%% Plots

figure(1)
t1 = tiledlayout('flow');
nexttile
plot(x, f(x), 'LineWidth', 1.1);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$f(x)$', 'Interpreter','latex');
title('$e^{-0.3x - 0.5x^2}(1 + sin^2(3x))$', 'Interpreter', 'latex', 'FontSize', 12);
grid on; grid minor;
for i = 1:ntries
    nexttile
    plot(x, ef(x,lambda(i)), 'LineWidth', 1.1);
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$f(x)$', 'Interpreter','latex');
    title(['$e^{\lambda f(x)}, \lambda =\ $', num2str(lambda(i))], 'Interpreter', 'latex', 'FontSize', 12);
    grid on; grid minor;
end
title(t1, 'Effect of changing $\lambda$ on $e^{\lambda f(x)}$', 'Interpreter', 'latex', 'FontSize', 14);

figure(2)
t1 = tiledlayout('flow');
nexttile
plot(x, f(x), 'LineWidth', 1.1);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$f(x)$', 'Interpreter','latex');
title('$e^{-0.3x - 0.5x^2}(1 + sin^2(3x))$', 'Interpreter', 'latex', 'FontSize', 12);
grid on; grid minor;
for i = 1:ntries
    nexttile
    plot(x, xef(x,lambda(i)), 'LineWidth', 1.1);
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$f(x)$', 'Interpreter','latex');
    title(['$xe^{\lambda f(x)}, \lambda =\ $', num2str(lambda(i))], 'Interpreter', 'latex', 'FontSize', 12);
    grid on; grid minor;
end
title(t1, 'Effect of changing $\lambda$ on $xe^{\lambda f(x)}$', 'Interpreter', 'latex', 'FontSize', 14);

figure(3)
subplot(3,1,1)
    plot(x, f(x), 'LineWidth', 1.1, 'DisplayName', '$f(x)$'); hold on;
    for i = 1:ntries
        stem(xo_direct(i), f(xo_direct(i)), 'filled', 'DisplayName', ['$\lambda =\ $', num2str(lambda(i))]);
    end
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$f(x)$', 'Interpreter','latex');
    title('Direct Integration of Closed Form');
    grid on; grid minor;
    hl = legend('show');
    set(hl, 'Interpreter','latex')
subplot(3,1,2)
    plot(x, f(x), 'LineWidth', 1.1, 'DisplayName', '$f(x)$'); hold on;
    for i = 1:ntries
        stem(xo_smirnov(i), f(xo_smirnov(i)), 'filled', 'DisplayName', ['$\lambda =\ $', num2str(lambda(i))]);
    end
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$f(x)$', 'Interpreter','latex');
    title('Smirnov Transform');
    grid on; grid minor;
    hl = legend('show');
    set(hl, 'Interpreter','latex')
sgtitle('Effect of changing \lambda on the optimizer');
