clearvars
close all
clc

% rng(1,"twister");

%% Setup a multimodal function

N    = 1000;
beta = 1;

syms x
% F  = -(0.2 * exp(-(x-6).^2) + 0.6 * exp(-x.^2) + 0.2 * exp(-(x+6).^2));
% F  = -abs(sinc(x));
F   = (x - 4).^2;
dF  = diff(F);
P   = exp(-beta * F);
Fx  = matlabFunction(F);
dFx = matlabFunction(dF);
Px  = matlabFunction(P);

%% Langevin Monte Carlo Updates

L = 1000;

xt  = unifrnd(-10, 10, [1, 1]);
eta = 0.1;

x_stores = zeros(L,1);

for l = 1:L
    xt = xt - eta * dFx(xt) + sqrt(2 * eta / beta) * randn(1,1);
    x_stores(l,1) = xt;
end

%% Simulated Tempering

beta = 1:-0.05:0.01;
nSwap   = 5;
nChains = numel(beta);

%% Find optimum

x = linspace(-10, 10, N);

[~, minInd] = min(Fx(x_stores));
x_opt = x_stores(minInd);

gradF = dFx(x_stores);

figure(1)
subplot(5,1,1)
    plot(x, Fx(x));
    grid on; grid minor;
    xlabel('x');
    ylabel('$f(x)$', 'Interpreter','latex');
    title('Function');
subplot(5,1,2)
    plot(x, Px(x));
    grid on; grid minor;
    xlabel('x');
    ylabel('$f(x)$', 'Interpreter','latex');
    title('Gibbs Density');
subplot(5,1,3)
    plot(x_stores, 1:L);
    grid on; grid minor;
    xlabel('x');
    ylabel('Iterations');
    xlim([-10 10]);
    set(gca, 'ydir','reverse')
    title('Sample Path');
subplot(5,1,4)
    plot(1:L, gradF);
    grid on; grid minor;
    xlabel('x');
    ylabel('$\nabla f(x)$', 'Interpreter','latex');
    title('Gradient');
subplot(5,1,5)
    histogram(x_stores);%, numel(U));
    xlim([-10,10]);
    title('Histogram');
    grid on; 

%% Add heating to distribution

function pTx = getTemperedDistribution(f, beta)
    pTx = exp(-1/beta * f(x));
end