classdef classGenPlots_SL < handle
    %CLASSGENPLOTS Summary of this class goes here
    %   Detailed explanation goes here

    properties

        numIter;
        numParticles;
        numParams;
        chirps = struct('sampleRate', 1000, ...
                        'numSamples', 100, ...
                        'numChirps',  1,  ...
                        'envelope',   [], ...
                        'chirpComponents', [], ...
                        'chirpMixed', [], ...
                        'chirpRecon', [], ...
                        'instFreq',   []);
        param;
        stepSize;
        grad;
        gradNorm;
        hessTrace;
        avgGradNorm;
        objFunc;
        temp;
        bAccept;
        noiseVar;
        paramProp; 

        % Compute
        instFreq;

        % Axis
        yAxisMin;
        yAxisMax;
    end

    methods
        function obj = classGenPlots_SL(simClass)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here

            obj.chirps.sampleRate = simClass.cpe{1,1}.fs;
            obj.chirps.numSamples = simClass.cpe{1,1}.N;
            obj.chirps.numChirps  = simClass.cpe{1,1}.Nc;
            obj.chirps.chirpMixed = simClass.cpe{1,1}.ym;
            obj.chirps.chirpRecon = simClass.cpe{1,simClass.bestParticleInd}.ymRecon;
            obj.chirps.chirpComponents = simClass.cpe{1,1}.xm;
            obj.chirps.instFreq   = simClass.cpe{1,1}.fim;

            obj.param       = simClass.saveParams;
            obj.grad        = simClass.saveGrads;
            obj.stepSize    = simClass.saveStepSize;
            obj.gradNorm    = simClass.saveGradNorm;
            obj.avgGradNorm = simClass.saveAvgGradNorm;
            obj.objFunc     = simClass.saveObjFunc;
            obj.temp        = simClass.saveTemper;
            obj.bAccept     = simClass.savebAccept;
            obj.noiseVar    = simClass.saveNoiseVar;
            obj.paramProp   = simClass.saveProposedParams;
            obj.hessTrace   = simClass.saveHessTrc;

            obj.numParticles = simClass.numParticles;
            obj.numIter = simClass.numIterLmc;
            obj.numParams = simClass.numParams;
            % obj.yAxisMin  = simClass.initValMinMax(1);
            % obj.yAxisMax  = simClass.initValMinMax(2);
        end

        function genAllPlots(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            fs = obj.chirps.sampleRate;
            N  = obj.chirps.numSamples;
            Nc = obj.chirps.numChirps;
            xm = obj.chirps.chirpComponents;
            ym = obj.chirps.chirpMixed;
            ymRecon = obj.chirps.chirpRecon;
            fim = obj.chirps.instFreq;
            bAc = obj.bAccept;

            Td = N / fs;
            tx = 0:1/fs:Td-1/fs;
            nx = 1:obj.numIter;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                         %
            % --- Plotting Chirps --- %
            %                         %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%

            figure('windowstyle','docked');
            tiledlayout flow;

            % Instantaneous Frequencies
            nexttile
                for cind = 1:Nc
                    plot(tx, fim(:,cind), 'LineWidth', 2, 'DisplayName', ['chirp ', num2str(cind)]); hold on;
                end
                hold off; grid on; grid minor;
                xlabel('Time (s)', 'FontSize', 12); 
                ylabel('Frequency (Hz)', 'FontSize', 12);
                title('Instantaneous Frequency vs Time', 'FontSize', 14);
                lgd = legend('show', 'Location', 'best');
                fontsize(lgd, 12, 'points');
            
             % Individual chirp components (no envelope)
            for cind = 1:Nc
                nexttile
                plot(tx, real(xm(:,cind)), 'LineWidth', 1.2);
                grid on; grid minor;
                xlabel('Time (s)', 'FontSize', 12);
                ylabel('Frequency (Hz)', 'FontSize', 12);
                title(['Chirp ', num2str(cind), ' vs Time'], 'FontSize', 14);
            end

            % Mixed chirp with envelope and added noise
            nexttile
                plot(tx, real(ym), 'LineWidth', 1.3);
                grid on; grid minor;
                xlabel('Time (s)', 'FontSize', 12);
                ylabel('Frequency (Hz)', 'FontSize', 12);
                title('Multicomponent Chirp vs Time', 'FontSize', 14);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                            %
            % --- Plotting Grad Norm --- %
            %                            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            figure('windowstyle','docked');
            tiledlayout flow;

            nexttile
                for pind = 1:obj.numParticles
                    gradp = obj.avgGradNorm(pind,:);
                    plot(nx, gradp, 'LineWidth', 1.1, 'DisplayName', ['particle ', num2str(pind)]); hold on;
                end
                hold off; grid on; grid minor;
                xlabel('Iterations', 'FontSize', 12); 
                ylabel('$$\nabla_{\varphi} J$$', 'FontSize', 20, 'Interpreter','latex');
                title('Average Gradient Norm vs Iterations', 'FontSize', 14);
                lgd = legend('show', 'Location', 'best');
                fontsize(lgd, 12, 'points');

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                                %
            % --- Plotting Accept/Reject --- %
            %                                %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            nexttile
                for pind = 1:obj.numParticles
                    bAVec = bAc(pind,:);
                    plot(nx, bAVec, 'LineWidth', 1.1, 'DisplayName', ['particle ', num2str(pind)]); hold on;
                end
                hold off; grid on; grid minor;
                xlabel('Iterations', 'FontSize', 12); 
                ylabel('Accept/Reject', 'FontSize', 20);
                title('Metropolis Accept/Reject vs Iterations', 'FontSize', 14);
                lgd = legend('show', 'Location', 'best');
                fontsize(lgd, 12, 'points');

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                                        %
            % --- Plotting Stepsize Trajectories --- %
            %                                        %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            figure('windowstyle','docked');
            tl = tiledlayout('flow');
            for prind = 1:obj.numParams
                nexttile
                    for pind = 1:obj.numParticles                   
                        mu = reshape(squeeze(obj.stepSize(prind, pind, :)), [], 1); hold on;
                        plot(nx, mu, 'LineWidth', 1.1, 'DisplayName', ['particle ', num2str(pind)]);
                        grid on; grid minor;
                        xlabel('Iterations', 'FontSize', 12); 
                        ylabel(['$$\eta_',num2str(prind),'$$'], 'FontSize', 20, 'Interpreter','latex');              
                    end
                lgd = legend('show', 'Location', 'best');
                fontsize(lgd, 12, 'points');
                title(['Parameter ', num2str(prind)], 'FontSize', 14);                 
            end   
            title(tl, 'Stepsize vs Iterations', 'FontSize', 16); 


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                              %
            % --- Plotting Temperature --- %
            %                              %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            figure('windowstyle','docked');
                plot(obj.temp(:,1), 'LineWidth', 1.2, 'DisplayName', ['particle ', num2str(pind)]);
                grid on; grid minor;
                xlabel('LMC Iterations', 'FontSize', 12); 
                ylabel('Temperature', 'FontSize', 12);
                title('Temperature vs Iterations', 'FontSize', 14);
                lgd = legend('show', 'Location', 'best');
                fontsize(lgd, 12, 'points');


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                                         %
            % --- Plotting Parameter Trajectories --- %
            %                                         %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            for prind = 1:obj.numParams
                ax(prind+1) = nexttile;
                for pind = 1:obj.numParticles                   
                    theta = reshape(squeeze(obj.param(prind, pind, :)), [], 1); 
                    thetaProp = reshape(squeeze(obj.paramProp(prind, pind, :, :)), [], 1);
                    hold on;
                    plot(nx, theta, 'LineWidth', 1.1, 'DisplayName', ['particle ', num2str(pind)]); 
                    % plot(nx, thetaProp,'LineWidth', 0.1);
                    grid on; grid minor;
                    xlabel('Iterations', 'FontSize', 12); 
                    ylabel(['$$\varphi_',num2str(prind),'$$'], 'FontSize', 20, 'Interpreter','latex');
                    % ylim([obj.yAxisMin, obj.yAxisMax]);
                end
                lgd = legend('show', 'Location', 'best');
                fontsize(lgd, 12, 'points');
                title(['Parameter ', num2str(prind)], 'FontSize', 14);                 
            end   
            title(tl, 'Parameter Trajectories vs Iterations', 'FontSize', 16);
            linkaxes(ax, 'x');

            % ----------------------------------- %
            % ----- For Paper Start ------------- %

            black  = [0, 0, 0];       % Black
            blue  = [0 0.4470 0.7410];
            red   = [0.8500 0.3250 0.0980];
            green = [0.4660 0.6740 0.1880];

            plotColors = {blue, blue, blue, red, red, red, green, green, green};

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                                 %
            % --- Plotting Noise Variance --- %
            %                                 %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            noiseVarInit  = 1;
            noiseVarFinal = 0.001;
            numIterNoise  = 10;
            samplesPerIter = 100;  % Number of samples per noise variance
            noiseVarRatio = (noiseVarInit / noiseVarFinal)^(1/(numIterNoise - 1)); % GP common ratio
            
            noiseVarNaLmc = zeros(1, numIterNoise);
            for ni = 1:numIterNoise
                noiseVarNaLmc(ni) = noiseVarInit / noiseVarRatio^(ni-1);
            end
            
            % Array to hold each value for 100 samples
            noiseVarSamples = repmat(noiseVarNaLmc', 1, samplesPerIter); % Repeat each value 100 times
            noiseVarSamples = noiseVarSamples';  % Transpose to match dimension
            noiseVarSamples = noiseVarSamples(:)';  % Flatten the matrix into a single array

            figure('windowstyle','docked');
            tiledlayout flow 
            nexttile
            for pind = 1:obj.numParticles
                plot(nx, obj.noiseVar(pind,:), 'LineWidth', 1.1, 'DisplayName', ['particle ', num2str(pind)], 'Color', plotColors{1,pind}); hold on;
                grid on; grid minor;
                xlabel('Iterations', 'FontSize', 12); 
                ylabel('$$\sigma$$', 'FontSize', 20, 'Interpreter','latex'); 
                title('Noise Variance ', 'FontSize', 14);     
            end
            plot(nx, noiseVarSamples, 'LineWidth', 1.1, 'DisplayName', 'NA-LMC', 'Color', black, 'LineStyle', '--');
            legend('show', 'Location', 'best'); 

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                                %
            % --- Plotting Hessian Trace --- %
            %                                %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            nexttile
                for pind = 1:obj.numParticles
                    plot(nx, log10(obj.hessTrace(pind,:)), 'LineWidth', 1.1, 'DisplayName', ['particle ', num2str(pind)], 'Color', plotColors{1,pind}); hold on;
                    grid on; grid minor;
                    xlabel('Iterations', 'FontSize', 12); 
                    ylabel('Trace of Hessian $\log_{10}$', 'FontSize', 20, 'Interpreter','latex'); 
                    % title('Trace of Hessian', 'FontSize', 14);     
                end
                legend('show', 'Location', 'best');

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                                     %
            % --- Plotting Objective Function --- %
            %                                     %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            nexttile;
                for pind = 1:obj.numParticles
                    objp = obj.objFunc(pind,:);
                    plot(nx, objp, 'LineWidth', 1.2, 'DisplayName', ['particle ', num2str(pind)], 'Color', plotColors{1,pind}); hold on;
                end
                
                hold off; grid on; grid minor;
                xlabel('Iterations', 'FontSize', 12); 
                ylabel('Objective Func', 'FontSize', 12);
                % title('Objective Function vs Iterations', 'FontSize', 14);
                lgd = legend('show', 'Location', 'best');
                fontsize(lgd, 12, 'points');

            % ----- For Paper End-- ------------- %
            % ----------------------------------- %
      
        end
    end
end