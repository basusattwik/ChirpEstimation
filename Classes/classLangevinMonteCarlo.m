classdef classLangevinMonteCarlo < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        
        % Hyperparameters
        stepSize;
        numIter; 
        temper;
        tempSwap;     
        avgConst;
        stepSizeConst;
        stepSizeMax;

        % Parameters 
        numParticles;
        numParams;
        param;
        grads;
        objectiveJ; 
        paramOptPersonalBest; % Optimum/Almost-optimum parameters
        paramOptCollectiveBest;
        numSamplesToUse;
        gradNorm;
        avgGradNorm;
        avgGrads; % This is a vector NumParam x 1

        % External
        cpe; % classChirpParamEst

        % For analysis
        saveObjectiveFunc;
        saveParams;
        saveGrads;
        saveGradNorm;
        saveInvTemp;
        savebAccept;
        saveStepSize;
    end

    methods
        function obj = classLangevinMonteCarlo(tuning, setup)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.numParams = setup.numParams;

            % Set tuning for hyperparams
            obj.numIter  = tuning.numIter;
            obj.tempSwap = tuning.tempSwap;
            obj.avgConst = tuning.avgConst;
            obj.numSamplesToUse = tuning.numSamplesToUse;
            obj.stepSizeConst   = tuning.stepSizeConst;
            obj.numParticles    = tuning.numParticles;

            obj.stepSize    = zeros(obj.numParams, obj.numParticles);
            obj.stepSizeMax = zeros(obj.numParams, obj.numParticles);
            obj.temper      = zeros(1, obj.numParticles);
            for np = 1:obj.numParticles
                obj.temper(1,np)   = tuning.initInvTemp;
                obj.stepSize(:,np) = tuning.stepSize;
                for npr = 1:obj.numParams
                    obj.stepSizeMax(npr,np) = tuning.stepSizeMax;
                end
            end

            % Analysis buffers
            obj.saveParams   = zeros(obj.numParams, obj.numParticles, obj.numIter);
            obj.saveGrads    = zeros(obj.numParams, obj.numParticles, obj.numIter);
            obj.saveGradNorm = zeros(obj.numParticles, obj.numIter);
            obj.saveInvTemp  = zeros(obj.numParticles, obj.numIter);
            obj.savebAccept  = zeros(obj.numParticles, obj.numIter);
            obj.saveStepSize = zeros(obj.numParams, obj.numParticles, obj.numIter);
            obj.saveObjectiveFunc = zeros(obj.numParticles, obj.numIter);

            % Init params and gradient
            paramInit = [5, 2, 1.5, 40].';
            obj.param = [repmat(paramInit, 1, obj.numParticles)];
            % obj.param = unifrnd(0, 50, obj.numParams, obj.numParticles);
            % 10*rand(obj.numParams, obj.numParticles); % what is a good init [5, 7, 5, 4, 0, 40, 0.5, 0, 65, 1].'; % 
            obj.grads = zeros(obj.numParams, obj.numParticles); 
            obj.avgGrads    = zeros(obj.numParams, obj.numParticles);
            obj.avgGradNorm = zeros(1, obj.numParticles);
            obj.gradNorm    = zeros(1, obj.numParticles);
            obj.paramOptPersonalBest   = zeros(obj.numParams, obj.numParticles);
            obj.paramOptCollectiveBest = zeros(obj.numParams, 1);

            %  Call constructor of cpe class. Initialize several instances
            %  of the class for each particle in the swarm
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

            try 
                for itr = 1:obj.numIter

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %
                    % --- Run MLE Model for Chirp Parameter Estimation ---
                    %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    % Call the process function for CPE. TODO: Parallelize?
                    for np = 1:obj.numParticles
                        obj.cpe{1,np} = obj.cpe{1,np}.runCpeCore(obj.param(:,np)); % This gives the gradients wrt all params for all particles
                    end

                    for np = 1:obj.numParticles

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %
                        % --- Simulated Tempering ---
                        %
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
                        % Run Simulated Tempering
                        if mod(itr, obj.tempSwap) == 0
                            obj.temper(1,np) = log(itr) / 0.7; % remove magic number
                        end
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %
                        % --- Stepsize Annealing ---
                        %
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
                        % Do Stepsize Annealing         
                        obj.avgGrads(:,np) = (1 - obj.avgConst) .* abs(obj.grads(:,np)) + obj.avgConst .* obj.avgGrads(:,np);
                        obj.stepSize(:,np) = obj.stepSize(:,np) ./ (1 + obj.stepSizeConst .* obj.avgGrads(:,np));
                        for npr = 1:obj.numParams
                            if obj.stepSize(npr,np) > obj.stepSizeMax(npr,np)
                                obj.stepSize(npr,np) = obj.stepSizeMax(npr,np);
                            end
                        end
        
                        % Keep track of average gradient norm
                        obj.gradNorm(1,np)    = norm(obj.grads(:,np))^2;
                        obj.avgGradNorm(1,np) = (1 - obj.avgConst) .* abs(obj.gradNorm(1,np)) + obj.avgConst .* obj.avgGradNorm(1,np);

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %
                        % --- Langevin Monte Carlo ---
                        %
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        % Do Langevin updates on all params
                        obj.grads(:,np) = [obj.cpe{1,np}.dJ_beta.'; obj.cpe{1,np}.dJ_gamma.'; obj.cpe{1,np}.dJ_phi.'];               
                        paramProp = obj.param(:,np) - obj.stepSize(:,np) .* obj.grads(:,np) + sqrt(2 * obj.stepSize(:,np) ./ obj.temper(1,np)) .* randn(obj.numParams, 1);     
                       
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %
                        % --- Apply Metropolis step ---
                        %
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
                        % Run Metropolis Adjustment
                        alpha = min(1, exp(-obj.temper * (obj.cpe{1,np}.evalObjectiveFunc(paramProp) - obj.cpe{1,np}.J)));
                        if rand(1,1) <= alpha
                            obj.param(:,np) = paramProp;
                            obj.savebAccept(np,itr) = 1;
                        end
        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %
                        % --- Save Data for Analysis ---
                        %
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                        
                        % Fill analysis arrays
                        obj.saveParams(:,np,itr) = obj.param(:,np);
                        obj.saveGrads(:,np,itr)  = obj.grads(:,np);
                        obj.saveGradNorm(np,itr) = obj.gradNorm(1,np);
                        obj.saveInvTemp(np,itr)  = obj.temper(1,np);
                        obj.saveStepSize(:,np,itr)    = obj.stepSize(:,np); 
                        obj.saveObjectiveFunc(np,itr) = obj.cpe{1,np}.J;

                    end % end numParticles

                    % Update waitbar and message
                    if getappdata(wbar, 'canceling')
                        disp('Simulation cancelled!')
                        break
                    end    
                    if mod(itr, 10) == 0
                        waitbar(itr / obj.numIter);
                    end
                end % end numIter

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                % --- Find the optimum parameter vector ---
                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Find "personal" best for each particle
                for np = 1:obj.numParticles
                    paramSet = obj.saveParams(:, np, obj.numIter-obj.numSamplesToUse+1:end); % choose parameters obtained in the last few iterations
                    
                    objectiveFuncValsForSamples = zeros(obj.numSamplesToUse, 1);
                    for nind = 1:obj.numSamplesToUse
                        objectiveFuncValsForSamples(nind, 1) = obj.cpe{1,np}.evalObjectiveFunc(paramSet(:, :, nind));
                    end 

                    [~, paramOptInd] = min(objectiveFuncValsForSamples);
                    obj.paramOptPersonalBest(:,np) = paramSet(:, :, paramOptInd);
                end

                % Find the "collective" best param from the "personal"
                % bests of all the particles
                objectiveFuncVals = zeros(obj.numParticles, 1);
                for np = 1:obj.numParticles
                    objectiveFuncVals(np, 1) = obj.cpe{1,np}.evalObjectiveFunc(obj.paramOptPersonalBest(:,np));
                end
                [~, paramOptParticleInd]   = min(objectiveFuncVals);
                obj.paramOptCollectiveBest = obj.paramOptPersonalBest(:, paramOptParticleInd);

    
                % Finally, compute the scalar gains
                % obj.cpe = obj.cpe.compScalarGains(); % This should now take the optimum param from above
                
                delete(wbar);
                disp('Simulation complete!');

            catch me
                delete(wbar); % Close wait bar when simulation errors out
                throw(me);
            end
        end
    end
end
