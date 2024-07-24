clearvars
close all
clc

% rng(1,"twister");

%% Setup a multimodal function

c = 0.6;

syms x
F   = -((1-c)/8 * exp(-(x-8).^2 / 0.4)+ (1-c)/4 * exp(-(x-4).^2 / 0.4) + ...
        c * exp(-x.^2 / 0.4) + ...
       (1-c)/4 * exp(-(x+4).^2 / 0.4) + (1-c)/8 * exp(-(x+8).^2 / 0.4));
dF  = diff(F);
Fx  = matlabFunction(F);
dFx = matlabFunction(dF);
Px   = @(T, x) exp(-1/T * Fx(x));

%% Langevin Monte Carlo Updates

L = 10000;
d = 0.45;
xt  = unifrnd(-10, 10, [1, 1]);

xinit = xt;
eta = 0.01;
% b = 0.1;
% gamma = 0.95;

x_stores = zeros(L,1);
invTempSave = zeros(L,1);
etaSave = zeros(L,1);
T = 1;
for l = 1:L
    % Simulated Tempering
    if mod(l,50) == 0
        T = d / log(l); %beta + exp(l * 0.000001);
    end

    % Stepsize Annealing
    % if mod(l,500) == 0
    %     eta = eta * (b + l)^(-gamma);
    % end

    % Langevin Update
    xp = xt - eta  * (1/T) * dFx(xt) + sqrt(2 * eta) * randn(1,1);

    % Metropolis Adjustment
    alpha = min(1, Px(T, xp) / Px(T, xt));
    if rand(1,1) <= alpha
        xt = xp;
    end

    x_stores(l,1)    = xt;    
    invTempSave(l,1) = 1/T;
    etaSave(l,1) = eta;

end

%% Find optimum

N = 10000;
x = linspace(-10, 10, N);

[~, minInd] = min(Fx(x_stores));
x_opt = x_stores(minInd);

gradF = dFx(x_stores);

Pend = exp(-(invTempSave(end)) * F);
Pini = exp(-(invTempSave(1)) * F);
Pxend  = matlabFunction(Pend);
Pxini  = matlabFunction(Pini);

figure('windowstyle','docked')
tiledlayout flow
nexttile
    plot(x, Fx(x));
    grid on; grid minor;
    xlabel('x');
    ylabel('$f(x)$', 'Interpreter','latex');
    title('Function');
nexttile
    yyaxis right
    plot(x, Pxend(x)); hold on;
    yyaxis left
    plot(x, Pxini(x));
    grid on; grid minor;
    xlabel('x');
    ylabel('$f(x)$', 'Interpreter','latex');
    title('Gibbs Density');
    legend('Initial', 'Final');
nexttile
    plot(x_stores, 1:L); hold on;
    stem(xinit, 0);
    grid on; grid minor;
    xlabel('x');
    ylabel('Iterations');
    xlim([-10 10]);
    set(gca, 'ydir','reverse')
    title('Sample Path');
nexttile
    histogram(x_stores);%, numel(U));
    xlim([-10,10]);
    title('Histogram');
    grid on; 
nexttile
    plot(1:L, invTempSave);
    grid on; grid minor;
    xlabel('Iterations');
    ylabel('beta');
    title('Inverse Temperature');
nexttile
    plot(1:L, gradF);
    grid on; grid minor;
    xlabel('x');
    ylabel('$\nabla f(x)$', 'Interpreter','latex');
    title('Gradient');
% nexttile
%     plot(1:L, etaSave);
%     grid on; grid minor;
%     xlabel('x');
%     ylabel('Stepsize');
%     title('Stepsize Annealing');


disp(['Initial inv. temp = ', num2str(invTempSave(1))]);
disp(['Final inv. temp   = ', num2str(invTempSave(end))]);