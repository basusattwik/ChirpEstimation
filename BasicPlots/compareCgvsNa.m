clc

% Load two data sets
numIter = 1000;
nx = 1:numIter;
numParticles = 9;
xAxMax = 500;

cglmcData = load('Data/Output/ForPaper/MATFiles/CG_LMC_Test_plt_100_dB_StRun_1.mat');
nalmcData = load('Data/Output/ForPaper/MATFiles/NA_LMC_plt_100_dB_StRun_1.mat');

cglmcSigma = cglmcData.saveSimData.classData.saveNoiseVar;
nalmcSigma = nalmcData.saveSimData.classData.noiseVar;

cglmcTrace = cglmcData.saveSimData.classData.saveHessTrc;

cglmcObjFunc = cglmcData.saveSimData.classData.saveObjFunc;
nalmcObjFunc = nalmcData.saveSimData.classData.saveObjFunc;

% ----------------------------------- %
% ----- For Paper Start ------------- %

black       = [0, 0, 0];       % Black
% darkGrey  = [0.4, 0.4, 0.4]; % Dark grey
% lightGrey = [0.7, 0.7, 0.7]; % Light grey
blue  = [0 0.4470 0.7410];
red   = [0.8500 0.3250 0.0980];
green = [0.4660 0.6740 0.1880];

plotColors = {blue, blue, blue, red, red, red, green, green, green};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     %
% --- Plotting Objective Function --- %
%                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure('windowstyle','docked');
% tl = tiledlayout('flow');
hfig = figure;
t = tiledlayout("flow");
naLmcLineStyle = ["--", "-", "-."];
nexttile;
    for pind = 1:numParticles
        objp = cglmcObjFunc(pind,:);
        plot(nx, objp, 'LineWidth', 1.3, 'DisplayName', ['particle ', num2str(pind)], 'Color', plotColors{1,pind}); hold on;
    end
    for pind = [1,3,2]
        objp = nalmcObjFunc(pind,:);
        % plot(nx, objp, 'LineWidth', 1.2, 'DisplayName', 'NA-LMC', 'Color', black, 'LineStyle', '--'); hold on;
        plot(nx, objp, 'LineWidth', 0.9, 'Color', black, 'LineStyle', naLmcLineStyle(pind)); hold on;
    end
    
    hold off; grid on; grid minor;
    xlabel('Iterations', 'FontSize',20,'FontName', 'Times'); 
    ylabel('Objective Func $J$', 'FontSize',20,'FontName', 'Times', 'Interpreter','latex');
    % title('Objective Function vs Iterations', 'FontSize', 14);
    xlim([0, xAxMax]);
    set(findall(hfig,'-property','FontSize'),'FontSize', 20, 'FontName', 'Times')
    % lgd = legend('show', 'Location', 'best');
    lgd1 = legend('CG-LMC, Near: Run 1', 'CG-LMC, Near: Run 2', 'CG-LMC, Near: Run 3', 'CG-LMC, Far: Run 1', ...
   'CG-LMC, Far: Run 2', 'CG-LMC, Far: Run 3', 'CG-LMC, Very Far: Run 1', 'CG-LMC, Very Far: Run 2', 'CG-LMC, Very Far: Run 3', 'NA-LMC, Near', 'NA-LMC, Far', 'NA-LMC Very Far', 'Location', 'best');
    fontsize(lgd1, 13, 'points');
    % fontsize(lgd1, 13, 'points');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 %
% --- Plotting Noise Variance --- %
%                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noiseVarInit   = 1;
noiseVarFinal  = 0.001;
numIterNoise   = 10;
samplesPerIter = 100;  % Number of samples per noise variance
noiseVarRatio  = (noiseVarInit / noiseVarFinal)^(1/(numIterNoise - 1)); % GP common ratio

noiseVarNaLmc = zeros(1, numIterNoise);
for ni = 1:numIterNoise
    noiseVarNaLmc(ni) = noiseVarInit / noiseVarRatio^(ni-1);
end

% Array to hold each value for 100 samples
noiseVarSamples = repmat(noiseVarNaLmc', 1, samplesPerIter); % Repeat each value 100 times
noiseVarSamples = noiseVarSamples';  % Transpose to match dimension
noiseVarSamples = noiseVarSamples(:)';  % Flatten the matrix into a single array


nexttile
for pind = 1:numParticles
    % plot(nx, cglmcSigma (pind,:), 'LineWidth', 1.1, 'DisplayName', ['particle ', num2str(pind)], 'Color', plotColors{1,pind}); hold on;
    plot(nx, cglmcSigma (pind,:), 'LineWidth', 1.3, 'Color', plotColors{1,pind}, 'DisplayName',''); hold on;
    grid on; grid minor;
    xlabel('Iterations', 'FontSize', 20, 'FontName', 'Times'); 
    ylabel('$$\sigma$$', 'FontSize',20,'FontName', 'Times', 'Interpreter','latex'); 
    % title('Noise Variance ', 'FontSize', 14);     
end
xlim([0, xAxMax]);
plot(nx, noiseVarSamples, 'LineWidth', 0.9, 'Color', black, 'LineStyle', '--');
legend('show', 'Location', 'best'); 
set(findall(hfig,'-property','FontSize'),'FontSize',20,'FontName', 'Times')
lgd2 = legend('', '', '', '', '', '', '', '', '', 'NA-LMC', 'Location', 'best');
fontsize(lgd2, 13, 'points');

% legend('show', 'Location', 'best');
% legend('NA-LMC', 'Location', 'best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                %
% --- Plotting Hessian Trace --- %
%                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nexttile
    for pind = 1:numParticles
        % plot(nx, log10(cglmcTrace(pind,:)), 'LineWidth', 1.1, 'DisplayName', ['particle ', num2str(pind)], 'Color', plotColors{1,pind}); hold on;
        plot(nx, log10(cglmcTrace(pind,:)), 'LineWidth', 1.3, 'Color', plotColors{1,pind}); hold on;
        grid on; grid minor;
        xlabel('Iterations', 'FontSize',20,'FontName', 'Times'); 
        ylabel('Abs. Trace\{Hessian\} ($\log_{10}$)', 'FontSize',20,'FontName', 'Times', 'Interpreter','latex'); 
        % title('Trace of Hessian', 'FontSize', 14);     
    end
    % legend('CG-LMC, Near: Run 1', 'CG-LMC, Near: Run 2', 'CG-LMC, Near: Run 3', 'CG-LMC, Far: Run 1', ...
    %    'CG-LMC, Far: Run 2', 'CG-LMC, Far: Run 3', 'CG-LMC, Very Far: Run 1', 'CG-LMC, Very Far: Run 2', 'CG-LMC, Very Far: Run 3', 'Location', 'best');
    xlim([0, xAxMax]);
    % legend('show', 'Location', 'best');
set(findall(hfig,'-property','FontSize'),'FontSize',20,'FontName', 'Times')
t.TileSpacing = 'tight';
t.Padding = 'compact';
% ----- For Paper End-- ------------- %
% ----------------------------------- %

fontsize(lgd1, 13, 'points');
fontsize(lgd2, 13, 'points');

%% Make nice figure

fname = 'myFig';
makeNiceFigure(hfig, fname);


function [] = makeNiceFigure(hfig, fname)

picturewidth = 22; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
% set(findall(hfig,'-property','FontSize'),'FontSize',20) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
print(hfig, fname,'-dpng','-painters')

end
