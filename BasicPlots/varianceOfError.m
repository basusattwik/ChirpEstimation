close all
clearvars
clc

%% Create dummy arrays


snr_axis = [0, 3, 6, 9, 12, 15, 18];

numSnr = numel(snr_axis);

baseline_1 = linspace(0, -10, numSnr);
baseline_2 = linspace(-2, -12, numSnr);
algorithm  = linspace(-3, -12, numSnr);

hfig = figure;
plot(snr_axis, baseline_1, '-*', 'MarkerSize', 11, 'LineWidth', 1.5, 'DisplayName', 'LMC'); hold on;
plot(snr_axis, baseline_2, '-s', 'MarkerSize', 11,  'LineWidth', 1.5, 'DisplayName', 'NA-LMC');
plot(snr_axis, algorithm,  '-+','MarkerSize', 11, 'LineWidth', 1.5, 'DisplayName', 'GS-LMC');
grid on; grid minor;
xlabel('SNR (dB)');
ylabel('Variance of Error ($\log_{10}$), $\varphi_{1,1}$');
xlim([0, snr_axis(end)]);
xticks(snr_axis);
legend('show');

fname = "voe";
makeNiceFigure(hfig, fname);
%% Helper

function [] = makeNiceFigure(hfig, fname)

picturewidth = 22; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
print(hfig, fname,'-dpng','-painters')

end