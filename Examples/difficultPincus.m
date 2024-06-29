close all
clearvars
clc

%% Function 

rho = [1, 2, 3, 5, 10, 20];

f  = @(x,y) 2*exp(-((x-2).^2/30 + (y-2).^2/10)) + 1.8*exp(-((x-8).^2/0.08 + (y-2).^2/0.08));
ef = @(x,y,rho) exp(rho * f(x,y));

x = -20:0.05:20;
y = -20:0.05:20;

[X,Y] = meshgrid(x,y);
F = f(X,Y);

EF  = zeros(numel(x), numel(x), numel(rho));
Fx = zeros(numel(x), numel(rho)); 
for i = 1:numel(rho)
    EF(:,:,i) = ef(X,Y,rho(i));
    Fx(:,i) = trapz(y, EF(:,:,i), 1); % Marginalize
end

% Pincus maximization
% xopt = trapz(x, Fx(x, lambda)) / trapz(x, ef(x, lambda)); 

%% Plots

figure(1)
surf(X,Y,F,'EdgeColor','none');
xlabel('x');
ylabel('y');
zlabel('f(x,y)');
title('Function');

figure(2)
tiledlayout('flow', 'TileSpacing', 'tight');
for i = 1:numel(rho)
    nexttile
    surf(x, y, EF(:,:,i), 'EdgeColor','none');
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16);
    zlabel('$e^{\rho f(x,y)}$', 'Interpreter', 'latex', 'FontSize', 16);
    title(['$\rho = \ $', num2str(rho(i))], 'Interpreter', 'latex', 'FontSize', 14);
end

figure(3)
tiledlayout('flow', 'TileSpacing', 'tight');
for i = 1:numel(rho)
    nexttile
    plot(x,Fx(:,i));
    grid on; 
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('$\int e^{\rho f(x,y)} dy$', 'Interpreter', 'latex', 'FontSize', 16);
    title(['$\rho = \ $', num2str(rho(i))], 'Interpreter', 'latex', 'FontSize', 14);
end