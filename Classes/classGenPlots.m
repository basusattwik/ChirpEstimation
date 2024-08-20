classdef classGenPlots < handle
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
                        'instFreq',   []);
        param;
        stepSize;
        grad;
        gradNorm;
        avgGradNorm;
        objFunc;
        temp;
        bAccept;

        % Compute
        instFreq;
    end

    methods
        function obj = classGenPlots(simClass)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here

            obj.chirps.sampleRate = simClass.cpe{1,1}.fs;
            obj.chirps.numSamples = simClass.cpe{1,1}.N;
            obj.chirps.numChirps  = simClass.cpe{1,1}.Nc;
            obj.chirps.chirpMixed = simClass.cpe{1,1}.ym;
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

            obj.numParticles = simClass.numParticles;
            obj.numIter = simClass.numIterLmc * simClass.numIterNoise;
            obj.numParams = simClass.numParams;
        end

        function genAllPlots(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            fs = obj.chirps.sampleRate;
            N  = obj.chirps.numSamples;
            Nc = obj.chirps.numChirps;
            xm = obj.chirps.chirpComponents;
            ym = obj.chirps.chirpMixed;
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
                    gradp = reshape(squeeze(obj.avgGradNorm(pind,:,:)), [], 1);
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
                    bAVec = reshape(squeeze(bAc(pind,:,:)), [], 1);
                    plot(nx, bAVec, 'LineWidth', 1.1, 'DisplayName', ['particle ', num2str(pind)]); hold on;
                end
                hold off; grid on; grid minor;
                xlabel('Iterations', 'FontSize', 12); 
                ylabel('Accept/Reject', 'FontSize', 20);
                title('Metropolis Accept/Reject vs Iterations', 'FontSize', 14);
                lgd = legend('show', 'Location', 'best');
                fontsize(lgd, 12, 'points');

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                                     %
            % --- Plotting Objective Function --- %
            %                                     %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            figure('windowstyle','docked');
                for pind = 1:obj.numParticles
                    objp    = reshape(squeeze(obj.objFunc(pind,:,:)), [], 1);
                    avgObjp = smoothdata(objp, 'sgolay', 100);
                    plot(nx, avgObjp, 'LineWidth', 1.2, 'DisplayName', ['particle ', num2str(pind)]); hold on;
                end
                hold off; grid on; grid minor;
                xlabel('Iterations', 'FontSize', 12); 
                ylabel('Objective Func', 'FontSize', 12);
                title('Objective Function vs Iterations', 'FontSize', 14);
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
                    mu = reshape(squeeze(obj.stepSize(prind, pind, :, :)), [], 1); hold on;
                    avgmu = smoothdata(mu, 'sgolay', 100);
                    plot(nx, avgmu, 'LineWidth', 1.1, 'DisplayName', ['particle ', num2str(pind)]);
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
                plot(obj.temp, 'LineWidth', 1.2, 'DisplayName', ['particle ', num2str(pind)]);
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

            figure('windowstyle','docked');
            tl = tiledlayout('flow');
            for prind = 1:obj.numParams
                nexttile
                for pind = 1:obj.numParticles                   
                    theta = reshape(squeeze(obj.param(prind, pind, :, :)), [], 1); hold on;
                    plot(nx, theta, 'LineWidth', 1.1, 'DisplayName', ['particle ', num2str(pind)]);
                    grid on; grid minor;
                    xlabel('Iterations', 'FontSize', 12); 
                    ylabel(['$$\varphi_',num2str(prind),'$$'], 'FontSize', 20, 'Interpreter','latex');              
                end
                lgd = legend('show', 'Location', 'best');
                fontsize(lgd, 12, 'points');
                title(['Parameter ', num2str(prind)], 'FontSize', 14);                 
            end   
            title(tl, 'Parameter Trajectories vs Iterations', 'FontSize', 16); 

        end
    end
end