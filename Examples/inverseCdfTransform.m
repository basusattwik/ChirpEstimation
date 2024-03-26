close all
clc
clearvars

mu = 0;
v  = 1;
dx = 0.001;
x  = -5:dx:5;

% Get CDF
fx = (1 / (sqrt(2*pi) * v)) * exp(-0.5 * ((x - mu) / v).^2);
Fx = cumtrapz(fx) * dx;

% Get inverse CDF
N   = 100000;
y   = unifrnd(0,1,[N,1]);
Fxi = interp1(Fx, x, y, 'spline','extrap');

figure
subplot(3,1,1)
    plot(x, fx);
    grid on;
    title('Gaussian PDF');
subplot(3,1,2)
    plot(x, Fx)
    grid on
    title('Gaussian CDF');
subplot(3,1,3)
    histogram(Fxi);
    grid on
    title('Histogram after Inv. CDF Sampling');

