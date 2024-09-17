close all

% Parameters
fs = 1000;      % Sampling frequency (Hz)
T = 1.2;         % Duration of signal (seconds)
t = 0:1/fs:T-1/fs;    % Time vector

% Chirp parameters
f0 = 10;        % Initial frequency (Hz) 50, 60, -90, 100 | 10, 40, -70, 120
a1 = 40;        % Quadratic coefficient
a2 = -70;         % Cubic coefficient
a3 = 105;         % 4th-order coefficient

% 4th-order phase
phase = 2 * pi * (f0 * t + a1 * t.^2 + a2 * t.^3 + a3 * t.^4);
finst = (f0 + 2 * a1 * t  + 3 * a2 * t.^2 + 4 * a3 * t.^3);

amp = exp(-4 *t);
% Generate the chirp signal
x = amp .* cos(phase);

% Plot the signal
hfig = figure;
plot(t, x, 'LineWidth', 2.0); hold on;
% plot(t, fs/2, '-.');
% title('Chirp');
xlabel('Time (s)');
ylabel('Amplitude');
grid on; 
fname = 'myfigure';
makeNiceFigure(hfig, fname)

hfig = figure;
plot(t, finst, 'LineWidth', 2.0); hold on;
% plot(t, fs/2, '-.');
title('Instantaneous Frequency');
xlabel('Time (s)');
ylabel('IF (Hz)');
grid on; 
makeNiceFigure(hfig, fname);
% 
% hfig = figure;  % save the figure handle in a variable
% t = 0:0.02:10; x = t.*sin(2*pi*t)+ 2*rand(1,length(t)); % data
% plot(t,x,'k-','LineWidth',1.5,'DisplayName','$\Omega(t)$');
% xlabel('time $t$ (s)')
% ylabel('$\Omega$ (V)')


makeNiceFigure(hfig, fname)

%% Make nice figure

function [] = makeNiceFigure(hfig, fname)

picturewidth = 22; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',26) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
print(hfig, fname,'-dpng','-painters')

end
