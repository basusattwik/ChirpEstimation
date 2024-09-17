% Data from the image, organized into a matrix
data = [
    -7.65054, -12.3999, -15.9525, -21.0674, -6.61883, -11.7038, -17.1383, -23.2619;
    -6.11949, -10.2244, -15.2463, -21.0987, -6.03528, -10.699, -15.5844, -20.989;
    -5.69664, -10.0131, -15.6111, -21.1292, -5.91137, -10.1419, -15.5433, -22.2436
];

% Transpose the data matrix so that columns form clusters of 3 bars
data = -10*data'; % Now each column will form a cluster

% Create a figure with the appropriate size for ICASSP format
hfig = figure('Units', 'centimeters', 'Position', [5, 5, 9, 6]); % 9 cm width, 6 cm height

% Bar chart with grouped bars
b = bar(data, 'grouped');

% Set custom colors for the 3 bars in each cluster
b(1).FaceColor = [0 0.4470 0.7410]; % Default blue
b(2).FaceColor = [0.8500 0.3250 0.0980]; % Default red-orange
b(3).FaceColor = [0.4940 0.1840 0.5560]; % Purple

% Customizing the axes labels and title with Times font
set(gca, 'XTickLabel', {'$\varphi_{11}$', '$\varphi_{12}$', '$\varphi_{13}$', '$\varphi_{14}$', ...
    '$\varphi_{21}$', '$\varphi_{22}$', '$\varphi_{23}$', '$\varphi_{24}$'}, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
xlabel('Phase parameters', 'FontName', 'Times');
ylabel('Variance of Errror $(-10\log_{10})$', 'FontName', 'Times', 'Interpreter','latex');
% title('Bar Chart with Clusters of 3 Bars', 'FontName', 'Times');

% Add legend for the bar groups with Times font
legend('CG-LMC', 'LMC', 'NA-LMC', 'Location', 'Best', 'FontName', 'Times');

% Display grid
grid on; grid minor;

% Adjust the figure to remove unnecessary white space and make it publication-ready
set(gca, 'LooseInset', max(get(gca, 'TightInset'), 0.02));

% Optional: save the figure as a high-resolution PNG or EPS file for publication
print('ICASSP_bar_chart', '-dpng', '-r300');  % Save as PNG with 300 DPI
% print('ICASSP_bar_chart', '-depsc', '-r300'); % Save as EPS if required

fname = 'myFig';
makeNiceFigure(hfig, fname);

%% Make nice figure


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