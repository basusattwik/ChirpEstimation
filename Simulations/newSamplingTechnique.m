close all
clearvars
clc

%% Setup a multi-modal function

N = 200;

v = linspace(-10, 10, N);
f = @(x) exp(-(x-4).^2/2) + exp(-(x+4).^2/6); % handle to function

A = [];
U = [];

% 1. Sample uniformly from a bounded region (x)
ui = unifrnd(v(1), v(end), 2*N, 1);
U  = [U; ui];
for i = 1:15
    
    % 2. Empirical pmf & cdf
    pu = f(U) / sum(f(U));
    cu = cumsum(pu);
    
    %3. inverse CDF sampling from pu ... replace with SIR
    x = interp1(cu + 0.00001*rand(size(cu)), U, rand(N,1), "nearest");

    A = [A; x];
    
    setA = unique(A); % eliminate duplicates
    for j = 1:numel(setA)    
        Atemp = setA(1:end ~= j); % remove setA(j) from the array for comparison
        var   = min(abs(setA(j) - Atemp));

        % DD says I need to use duplicates as well...
    
        % Get two samples from a normal centered at A(n) with variance
        % equal to the distance to the closest point to A(n)
        y = normrnd(setA(j), var, 2, 1);
        A = [A; y];       
    end
    U = sort([U; A]);
    A = [];
end

%% Plot

figure(1)
subplot(2,1,1)
plot(v, f(v));
grid on; grid minor;
xlabel('x');
ylabel('$f(x)$', 'Interpreter','latex');
title('Gaussian Mixture');

subplot(2,1,2)
histogram(U);%, numel(U));
xlim([-10,10]);
title('Histogram');
grid on; 
