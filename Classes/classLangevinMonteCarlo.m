classdef classLangevinMonteCarlo < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties

        % Hyperparameters
        stepSize;
        stepSizeInit;
        stepSizeConst;
        stepSizeMax;
        initTemp;
        tempConst;
        avgConst;
        numParticles;
        numIterNoise;
        numIterLmc;
        numIterLmcAndNoise;
        numIterGradSmooth;
        noiseVarInit;
        noiseVarFinal;

        % Parameters 
        optParam;
        numParams;
        noiseVar; 
        noiseVarRatio;

        % States
        temper;
        param;
        initParam;
        grads;
        objectiveJ; 
        gradNorm;
        avgGradNorm;
        avgGrads; 

        % Bools
        bDisplayPlots;
        bStopLmc;

        % External
        cpe; % classChirpParamEst

        % For analysis
        saveObjFunc;
        saveParams;
        saveGrads;
        saveGradNorm;
        saveAvgGradNorm;
        saveTemper;
        savebAccept;
        saveStepSize;
    end

    methods
        function obj = classLangevinMonteCarlo(tuning, setup)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.numParams = setup.numParams;
            obj.bStopLmc  = false;

            % Set tuning for hyperparams
            obj.stepSizeConst = tuning.stepSizeConst;
            obj.tempConst     = tuning.tempConst;
            obj.avgConst      = tuning.avgConst;
            obj.numParticles  = tuning.numParticles;
            obj.numIterNoise  = tuning.numIterNoise;
            obj.numIterLmc    = tuning.numIterLmc;
            obj.noiseVarInit  = tuning.noiseVarInit;
            obj.noiseVarFinal = tuning.noiseVarFinal; 
            obj.bDisplayPlots = tuning.bDisplayPlots;
            obj.numIterGradSmooth = tuning.numIterGradSmooth;

            obj.stepSizeInit = zeros(obj.numParams, obj.numParticles);
            obj.stepSize     = zeros(obj.numParams, obj.numParticles);
            obj.stepSizeMax  = zeros(obj.numParams, obj.numParticles);
            obj.temper       = zeros(1, 1);

            for pind = 1:obj.numParticles
                obj.stepSizeInit(:,pind) = tuning.stepSize(1:obj.numParams);
                for npr = 1:obj.numParams
                    obj.stepSizeMax(npr,pind) = tuning.stepSizeMax;
                end
            end

            % Setup sequence of noise variances
            obj.noiseVarRatio = (obj.noiseVarInit / obj.noiseVarFinal)^(1/(obj.numIterNoise - 1)); % GP common ratio
            
            obj.noiseVar = zeros(1, obj.numIterNoise);
            for ni = 1:obj.numIterNoise
                obj.noiseVar(1,ni) = obj.noiseVarInit / (obj.noiseVarRatio)^(ni-1);
            end

            obj.numIterLmcAndNoise = obj.numIterLmc * obj.numIterNoise;

            % Init state params and gradients
            obj.param       = unifrnd(10, 60, obj.numParams, obj.numParticles);
            obj.grads       = zeros(obj.numParams, obj.numParticles); 
            obj.avgGrads    = zeros(obj.numParams, obj.numParticles);
            obj.avgGradNorm = zeros(1, obj.numParticles);
            obj.gradNorm    = zeros(1, obj.numParticles);

            obj.initParam = obj.param; % just store this for use during analysis

            % Preallocate arrays to store state data for later analysis
            obj.saveParams      = zeros(obj.numParams,    obj.numParticles, obj.numIterLmc,   obj.numIterNoise);
            obj.saveGrads       = zeros(obj.numParams,    obj.numParticles, obj.numIterLmc,   obj.numIterNoise);
            obj.saveStepSize    = zeros(obj.numParams,    obj.numIterLmc,   obj.numParticles, obj.numIterNoise);
            obj.saveGradNorm    = zeros(obj.numParticles, obj.numIterLmc,   obj.numIterNoise);
            obj.saveAvgGradNorm = zeros(obj.numParticles, obj.numIterLmc,   obj.numIterNoise);
            obj.savebAccept     = zeros(obj.numParticles, obj.numIterLmc,   obj.numIterNoise);
            obj.saveObjFunc     = zeros(obj.numParticles, obj.numIterLmc,   obj.numIterNoise);
            obj.saveTemper      = zeros(obj.numIterLmc, obj.numIterNoise);

            %  Call constructor of cpe class. Initialize several instances
            %  of the class... one for each particle
            obj.cpe = cell(1,obj.numParticles);
            for pind = 1:obj.numParticles
                obj.cpe{1,pind} = classChirpParamEst(setup);
            end
        end

        % Core Simulation Function
        function obj = runLmcCore(obj)
            %RUNLMCCORE Summary of this method goes here
            %   Detailed explanation goes here

            try 

                % Create a wait bar display
                waitbarIterCount = 1;
                wbar = waitbar(0, 'Sit tight!', 'Name','Running Langevin Monte Carlo', ...
                                  'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');
                totalNumIter = obj.numIterNoise * obj.numIterLmc * obj.numParticles;
                
                % Initialize a live plotter
                if obj.bDisplayPlots
                    ax = initLivePlots(obj);
                end

                totalGrads = zeros(obj.numParams, 1);
                for nind = 1:obj.numIterNoise

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %                             %
                    % --- Simulated Tempering --- %
                    %                             %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    % Save stepsize numerator
                    stepSizeNum = obj.stepSizeInit .* (obj.noiseVar(1, nind) / obj.noiseVarFinal)^2;

                    for tind = 1:obj.numIterLmc 
                        for pind = 1:obj.numParticles 

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %                                                      %
                            % --- Run MLE Model for Chirp Parameter Estimation --- %
                            %                                                      %
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                            % Gradient smoothing
                            for sind = 1:obj.numIterGradSmooth
                                obj.cpe{1,pind} = obj.cpe{1,pind}.runCpeCore(obj.param(:,pind) + obj.noiseVar(1, nind) .* randn(obj.numParams,1)); % This gives the gradients wrt all params for all particles
                                totalGrads = totalGrads + obj.cpe{1,pind}.dJ_phi.'; 
                            end
                            obj.grads(:,pind) = (totalGrads / obj.numIterGradSmooth);
                            totalGrads(:) = 0;

                            % Objective function for current param for current particle
                            obj.cpe{1,pind}.J = obj.cpe{1,pind}.evalObjectiveFunc(obj.param(:,pind));  
                            if obj.cpe{1,pind}.J <= obj.cpe{1,pind}.minObjTol
                                obj.cpe{1,pind}.bMinFound = true; % Check if J is close to 0... that means we have converged
                            end
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %                            %
                            % --- Stepsize Annealing --- %
                            %                            %
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
                            % Compute gradient smoothing         
                            obj.avgGrads(:,pind) = obj.avgConst .* obj.grads(:,pind) + (1 - obj.avgConst) .* obj.avgGrads(:,pind);
                            
                            % Keep track of average gradient norm
                            obj.gradNorm(1,pind)    = norm(obj.grads(:,pind))^2;
                            obj.avgGradNorm(1,pind) = obj.avgConst .* abs(obj.gradNorm(1,pind)) + (1 - obj.avgConst) .* obj.avgGradNorm(1,pind);
                            
                            % Update stepsize
                            obj.stepSize(:,pind) = stepSizeNum(:,pind) ./ (1e-5 + obj.stepSizeConst .* obj.avgGrads(:,pind).^2);
                            
                            % Check if stepsize exceeds limits
                            for npr = 1:obj.numParams
                                if obj.stepSize(npr,pind) > obj.stepSizeMax(npr,pind)
                                    obj.stepSize(npr,pind) = obj.stepSizeMax(npr,pind);
                                end
                            end
    
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %                              %
                            % --- Langevin Monte Carlo --- %
                            %                              %
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                         
    
                            % Adjust the temperature
                            obj.temper = 1;%obj.tempConst * log10(1 + (tind));

                            % Do Langevin updates on all params to get a new proposed point
                            paramProp = obj.param(:,pind) - obj.stepSize(:,pind) .* obj.grads(:,pind) + ...
                                                            sqrt(2 * obj.stepSize(:,pind) ./ obj.temper) .* randn(obj.numParams,1) + ...
                                                            obj.noiseVar(1, nind) .* randn(obj.numParams,1);   

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %                                   %
                            % --- Run Metropolis Correction --- %
                            %                                   %
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                            alpha = min(1, exp(-obj.temper * (obj.cpe{1,pind}.evalObjectiveFunc(paramProp) - obj.cpe{1,pind}.J)));
                            if rand(1,1) <= alpha
                                obj.param(:,pind) = paramProp;
                                bAccept = true;
                            else
                                bAccept = false;
                            end

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %                                %
                            % --- Save Data for Analysis --- %
                            %                                %
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                            
                            % Fill analysis arrays
                            obj.saveParams(:,pind, tind, nind)    = obj.param(:,pind);
                            obj.saveGrads(:, pind, tind, nind)    = obj.grads(:,pind);
                            obj.saveStepSize(:,pind,tind,nind)    = obj.stepSize(:,pind); 
                            obj.saveGradNorm(pind, tind, nind)    = obj.gradNorm(1,pind);
                            obj.saveAvgGradNorm(pind, tind, nind) = obj.avgGradNorm(1,pind);
                            obj.saveObjFunc(pind, tind, nind)     = obj.cpe{1,pind}.J;
                            obj.savebAccept(pind,tind,nind)       = bAccept;

                            %%%%%%%%%%%%%%%%%%%
                            %                 %
                            % --- Waitbar --- %
                            %                 %
                            %%%%%%%%%%%%%%%%%%%

                            % Increment iteration count for waitbar
                            waitbarIterCount = waitbarIterCount + 1;

                            % Update waitbar and message
                            if getappdata(wbar, 'canceling')
                                disp('Simulation cancelled!')
                                delete(wbar);
                                return
                            end    
                            if mod(waitbarIterCount, 10) == 0
                                waitbar(waitbarIterCount / totalNumIter, wbar, ['Noise Iteration: ', num2str(nind), '/', num2str(obj.numIterNoise)]);
                            end

                        end % end numParticles
                      
                        % Save temperature
                        obj.saveTemper(tind,nind) = obj.temper;

                        % Stopping criterion
                        for pind = 1:obj.numParticles
                            if obj.cpe{1,pind}.bMinFound
                                obj.bStopLmc = true;
                                break; 
                            end
                        end
                    end % end numIterLmc

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %                                                     %
                    % --- Starting point for the next noise iteration --- %
                    %                                                     %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

                    % Get the time index for which the param led to the lowest objective func value
                    for pind = 1:obj.numParticles
                        [~, minObjLmcInd] = min(obj.saveObjFunc(pind, 1:tind, nind));
                        obj.param(:,pind) = obj.saveParams(:,pind, minObjLmcInd, nind); % This is the starting point for the next round of LMC updates
                    end

                    % Update Live Plot
                    if obj.bDisplayPlots
                        updateLivePlots(obj, tind, nind, ax);
                    end

                    % Stop iterations since min has been found
                    if obj.bStopLmc
                        disp(['Minimum has been found already. LMC stopped after ', num2str(nind), ' iterations.']);
                        break;
                    end
                    
                end % end numIterNoise

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                                           %
                % --- Find the optimum parameter vector --- %
                %                                           %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                objFuncVals = zeros(obj.numParticles, 1);
                for pind = 1:obj.numParticles
                    objFuncVals(pind, 1) = obj.cpe{1,pind}.evalObjectiveFunc(obj.param(:,pind));
                end

                [~, minObjFuncInd] = min(objFuncVals);
                obj.optParam = obj.param(:, minObjFuncInd);
                disp(['The optimum chirp parameters are =  ', num2str(obj.optParam.')]);

                % Finally, compute the scalar gains for CPE
                obj.cpe{1,minObjFuncInd} = obj.cpe{1,minObjFuncInd}.compScalarGains(obj.optParam);   

                % Store best particle index
                % Compute b vector (amp and phase too)

                % Housekeeping
                disp('Simulation complete!')
                delete(wbar);

            catch me
                delete(wbar); % Close wait bar when simulation errors out
                throw(me);
            end
        end

        % Initialize live plots of the objective function
        function ax = initLivePlots(obj)
            figure('windowstyle','docked');
                ax = plot(1:obj.numIterLmcAndNoise, -ones(obj.numIterLmcAndNoise,obj.numParticles), 'LineWidth', 1.5); 
                grid on; grid minor;
                xlabel('Iterations', 'FontSize', 12); 
                ylabel('Objective Func', 'FontSize', 12); 
                xlim([0, obj.numIterLmcAndNoise]);
                ylim([0 inf]);
                title('Objective Function vs Iterations', 'FontSize', 14);
        end

        % Generate a live plot of the objective function as the sim progresses
        function updateLivePlots(obj, tind, nind, ax)
            for pind = 1:obj.numParticles
                objp    = reshape(squeeze(obj.saveObjFunc(pind,1:tind, 1:nind)), tind * nind, 1);
                avgObjp = smoothdata(objp, 'sgolay', 20);
                set(ax(pind), 'XData', 1:tind*nind, 'YData', avgObjp);
                set(ax(pind), 'DisplayName', ['particle ', num2str(pind)]);
                lgd = legend('show', 'Location', 'best');
                fontsize(lgd, 12, 'points');
            end
        end
    end
end
