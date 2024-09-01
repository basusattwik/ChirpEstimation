
clearvars
% close all

fx = @ (x, b, g) exp(-b * x) .* (1 - exp(-g * x));

fs = 500;
T  = 1;
N  = fs * T;
t  = 0:1/fs:T-1/fs;
b  = 8;
g  = 5;

f = fx(t,b,g);

c = polyfit(t, f, 10);
p = polyval(c, t);

figure('WindowStyle','docked')

plot(t, f); hold on;
plot(t, p);

