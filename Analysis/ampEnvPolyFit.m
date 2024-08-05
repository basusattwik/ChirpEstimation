
clearvars
close all

fx = @ (x, b, g) exp(-b * x) .* (1 - exp(-g * x));

fs = 100;
T  = 1;
N  = fs * T;
t  = 0:1/fs:T-1/fs;
b  = 4;
g  = 20;

f = fx(t,b,g);

c = polyfit(t, f, 2);
p = polyval(c, t);

figure('WindowStyle','docked')

plot(t, f); hold on;
plot(t, p);