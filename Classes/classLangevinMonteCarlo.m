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
        grads;
        objectiveJ; 
        gradNorm;
        avgGradNorm;
        avgGrads; % This is a vector NumParam x 1

        % External
        cpe; % classChirpParamEst

        % For analysis
        saveObjFunc;
        saveParams;
        saveGrads;
        saveGradNorm;
        saveTemp;
        savebAccept;
        saveStepSize;
    end

    methods
        function obj = classLangevinMonteCarlo(tuning, setup)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.numParams = setup.numParams;

            % Set tuning for hyperparams
            obj.stepSizeConst = tuning.stepSizeConst;
            obj.tempConst     = tuning.tempConst;
            obj.avgConst      = tuning.avgConst;
            obj.numParticles  = tuning.numParticles;
            obj.numIterNoise  = tuning.numIterNoise;
            obj.numIterLmc    = tuning.numIterLmc;
            obj.noiseVarInit  = tuning.noiseVarInit;
            obj.noiseVarFinal = tuning.noiseVarFinal; 

            obj.stepSizeInit = zeros(obj.numParams, obj.numParticles);
            obj.stepSize     = zeros(obj.numParams, obj.numParticles);
            obj.stepSizeMax  = zeros(obj.numParams, obj.numParticles);
            obj.temper       = zeros(1, obj.numIterNoise);

            for np = 1:obj.numParticles
                obj.stepSizeInit(:,np) = tuning.stepSize;
                for npr = 1:obj.numParams
                    obj.stepSizeMax(npr,np) = tuning.stepSizeMax;
                end
            end

            % Setup sequence of noise variances
            obj.noiseVarRatio = (obj.noiseVarInit / obj.noiseVarFinal)^(1/(obj.numIterNoise - 1)); % GP common ratio
            
            obj.noiseVar = zeros(1, obj.numIterNoise);
            for ni = 1:obj.numIterNoise
                obj.noiseVar(1,ni) = obj.noiseVarInit / (obj.noiseVarRatio)^(ni-1);
            end

            % Init state params and gradients
            obj.param                  = unifrnd(0, 50, obj.numParams, obj.numParticles);
            obj.grads                  = zeros(obj.numParams, obj.numParticles); 
            obj.avgGrads               = zeros(obj.numParams, obj.numParticles);
            obj.avgGradNorm            = zeros(1, obj.numParticles);
            obj.gradNorm               = zeros(1, obj.numParticles);

            % Preallocate arrays to store state data for later analysis
            obj.saveParams   = zeros(obj.numParams,    obj.numParticles, obj.numIterLmc,   obj.numIterNoise);
            obj.saveGrads    = zeros(obj.numParams,    obj.numParticles, obj.numIterLmc,   obj.numIterNoise);
            obj.saveStepSize = zeros(obj.numParams,    obj.numIterLmc,   obj.numParticles, obj.numIterNoise);
            obj.saveGradNorm = zeros(obj.numParticles, obj.numIterLmc,   obj.numIterNoise);
            obj.savebAccept  = zeros(obj.numParticles, obj.numIterLmc,   obj.numIterNoise);
            obj.saveObjFunc  = zeros(obj.numParticles, obj.numIterLmc,   obj.numIterNoise);
            obj.saveTemp     = zeros(1, obj.numIterNoise);

            %  Call constructor of cpe class. Initialize several instances
            %  of the class... one for each particle
            obj.cpe = cell(1,obj.numParticles);
            for np = 1:obj.numParticles
                obj.cpe{1,np} = classChirpParamEst(setup);
            end
        end

        function obj = runLmcCore(obj)
            %RUNLMCCORE Summary of this method goes here
            %   Detailed explanation goes here

            % Create a wait bar display
            wbar = waitbar(0, 'Sit tight!', 'Name','Running Langevin Monte Carlo', ...
                              'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');
            totalNumIter = obj.numIterNoise * obj.numIterLmc * obj.numParticles;
            iterCount    = 1;

            try 
                for nind = 1:obj.numIterNoise

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %                           %
                    % --- Simulated Tempering ---
                    %                           %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    % Adjust the temperature
                    obj.temper(1,nind) = obj.tempConst * log10(1 + (nind));

                    % Save stepsize numerator
                    stepSizeNum = obj.stepSizeInit .* (obj.noiseVar(1, nind) / obj.noiseVarFinal)^2;

                    for tind = 1:obj.numIterLmc 

                        for pind = 1:obj.numParticles 

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %                                                    %
                            % --- Run MLE Model for Chirp Parameter Estimation ---
                            %                                                    %
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                            % Get all the gradients
                            obj.cpe{1,pind}   = obj.cpe{1,pind}.runCpeCore(obj.param(:,pind)); % This gives the gradients wrt all params for all particles
                            obj.grads(:,pind) = obj.cpe{1,pind}.dJ_phi.'; 

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %                          %
                            % --- Stepsize Annealing ---
                            %                          %
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
                            % Compute gradient smoothing         
                            obj.avgGrads(:,pind) = obj.avgConst .* obj.grads(:,pind) + (1 - obj.avgConst) .* obj.avgGrads(:,pind);
    
                            % Update stepsize
                            obj.stepSize(:,pind) = stepSizeNum(:,pind) ./ (1 + obj.stepSizeConst .* obj.avgGrads(:,pind).^2);
                            
                            % Check if stepsize exceeds limits
                            for npr = 1:obj.numParams
                                if obj.stepSize(npr,pind) > obj.stepSizeMax(npr,pind)
                                    obj.stepSize(npr,pind) = obj.stepSizeMax(npr,pind);
                                end
                            end
            
                            % Keep track of average gradient norm
                            obj.gradNorm(1,pind)    = norm(obj.grads(:,pind))^2;
                            obj.avgGradNorm(1,pind) = obj.avgConst .* abs(obj.gradNorm(1,pind)) + (1-obj.avgConst) .* obj.avgGradNorm(1,pind);
    
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %
                            % --- Metropolis Adjusted Langevin Monte Carlo ---
                            %
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           
    
                            % Do Langevin updates on all params to get a new proposed point
                            paramProp = obj.param(:,pind) - obj.stepSize(:,pind) .* obj.grads(:,pind) + ...
                                                            sqrt(2 * obj.stepSize(:,pind) ./ obj.temper(1,pind)) .* randn(obj.numParams,1) + ...
                                                            obj.noiseVar(1, nind) .* randn(obj.numParams,1);                
            
                            % Run Metropolis Adjustment
                            alpha = min(1, exp(-(obj.temper(1,pind)) * (obj.cpe{1,pind}.evalObjectiveFunc(paramProp) - obj.cpe{1,pind}.J)));
    
                            if rand(1,1) <= alpha
                                obj.param(:,pind) = paramProp;
                                obj.savebAccept(pind,tind,nind) = 1;
                            end

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %
                            % --- Save Data for Analysis ---
                            %
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                            
                            % Fill analysis arrays
                            obj.saveParams(:,pind, tind, nind) = obj.param(:,pind);
                            obj.saveGrads(:, pind, tind, nind) = obj.grads(:,pind);
                            obj.saveStepSize(:,pind,tind,nind) = obj.stepSize(:,pind); 
                            obj.saveGradNorm(pind, tind, nind) = obj.gradNorm(1,pind);
                            obj.saveObjFunc(pind,  tind, nind) = obj.cpe{1,pind}.J;


                            %%%%%%%%%%%%%%%%%
                            %
                            % --- Waitbar ---
                            %
                            %%%%%%%%%%%%%%%%%

                            % Increment iteration count for waitbar
                            iterCount = iterCount + 1;

                            % Update waitbar and message
                            if getappdata(wbar, 'canceling')
                                disp('Simulation cancelled!')
                                break
                            end    
                            if mod(iterCount, 10) == 0
                                waitbar(iterCount / totalNumIter);
                            end

                        end % end numParticles

                    end % end numIterLmc

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %
                    % --- Starting point for the next noise iteration ---
                    %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

                    % Get the time index for which the param led to the lowest objective func value
                    for pind = 1:obj.numParticles
                        [~, minObjLmcInd] = min(obj.saveObjFunc(pind, :, nind));
                        obj.param(:,pind) = obj.saveParams(:,pind, minObjLmcInd, nind); % This is the starting point for the next round of LMC updates
                    end

                end % end numIterNoise

                % Save temperature
                obj.saveTemp(1,nind) = obj.temper(1,nind);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 
                % --- Find the optimum parameter vector ---
                % 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                objFuncVals = zeros(obj.numParticles, 1);
                for pind = 1:obj.numParticles
                    objFuncVals(pind, 1) = obj.cpe{1,pind}.evalObjectiveFunc(obj.param(:,pind));
                end

                [~, minObjFuncInd] = min(objFuncVals);
                obj.optParam = obj.param(:, minObjFuncInd);
                disp(obj.optParam)

                % Finally, compute the scalar gains for CPE
                obj.cpe{1,minObjFuncInd} = obj.cpe{1,minObjFuncInd}.compScalarGains(obj.optParam); 
                
                % Housekeeping!
                delete(wbar);
                disp('Simulation complete!');

            catch me
                delete(wbar); % Close wait bar when simulation errors out
                throw(me);
            end
        end
    end
end
