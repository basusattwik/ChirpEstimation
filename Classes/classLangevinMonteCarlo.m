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
        numParams;
        param;
        grads;
        objectiveJ; 
        paramOpt; % Optimum/Almost-optimum parameters
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
            obj.stepSize = tuning.stepSize;
            obj.stepSizeMax = tuning.stepSizeMax;
            obj.numIter  = tuning.numIter;
            obj.temper   = tuning.initInvTemp;
            obj.tempSwap = tuning.tempSwap;
            obj.avgConst = tuning.avgConst;
            obj.numSamplesToUse = tuning.numSamplesToUse;
            obj.stepSizeConst   = tuning.stepSizeConst;

            % Analysis buffers
            obj.saveParams   = zeros(obj.numParams, obj.numIter);
            obj.saveGrads    = zeros(obj.numParams, obj.numIter);
            obj.saveGradNorm = zeros(obj.numIter, 1);
            obj.saveInvTemp  = zeros(obj.numIter, 1);
            obj.savebAccept  = zeros(obj.numIter, 1);
            obj.saveStepSize = zeros(obj.numParams, obj.numIter);
            obj.saveObjectiveFunc = zeros(obj.numIter, 1);

            % Init params and gradient
            obj.param = [5, 7, 5, 4, 0, 40, 0.5, 0, 65, 1].'; % randn(obj.numParams, 1); % what is a good init [0, 0, 20000, 20000, 0, 40, 0, 65].'
            obj.grads = zeros(obj.numParams, 1); 
            obj.avgGrads    = zeros(obj.numParams, 1);
            obj.avgGradNorm = 0;
            obj.gradNorm    = 0;

            % Call constructor of cpe class
            obj.cpe = classChirpParamEst(setup);
        end

        function obj = runLmcCore(obj)
            %RUNLMCCORE Summary of this method goes here
            %   Detailed explanation goes here
           

            % Create a wait bar display
            wbar = waitbar(0, 'Please wait', 'Name','Running Langevin Monte Carlo', ...
                              'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');

            try 
                for itr = 1:obj.numIter
    
                    % parfor np = 1:numParticles
                    % end
                    % Call the process function for CPE
                    obj.cpe = obj.cpe.runCpeCore(obj.param); % This gives the gradients wrt all params
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % --- Langevin Monte Carlo ---
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                    % Do Langevin updates on all params
                    obj.grads = [obj.cpe.dJ_beta.'; obj.cpe.dJ_gamma.'; obj.cpe.dJ_phi.'];               
                    paramProp = obj.param - obj.stepSize .* obj.grads + sqrt(2 * obj.stepSize / obj.temper) .* randn(obj.numParams, 1);     
                   
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % --- Apply Metropolis step ---
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                    % Run Metropolis Adjustment
                    alpha = min(1, exp(-obj.temper * (obj.cpe.evalObjectiveFunc(paramProp) - obj.cpe.J)));
                    if rand(1,1) <= alpha
                        obj.param = paramProp;
                        obj.savebAccept(itr,1) = 1;
                    end
    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % --- Simulated Tempering ---
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                    % Run Simulated Tempering
                    if mod(itr, obj.tempSwap) == 0
                        obj.temper = log(itr) / 0.9; % remove magic number
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % --- Stepsize Annealing ---
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                    % Do Stepsize Annealing         
                    obj.avgGrads = (1 - obj.avgConst) .* abs(obj.grads) + obj.avgConst .* obj.avgGrads;
                    obj.stepSize = obj.stepSize ./ (1 + obj.stepSizeConst .* obj.avgGrads);
                    if obj.stepSize > obj.stepSizeMax
                        obj.stepSize = obj.stepSizeMax;
                    end
    
                    % Keep track of average gradient norm
                    obj.gradNorm    = norm(obj.grads);
                    obj.avgGradNorm = (1 - obj.avgConst) .* abs(obj.gradNorm) + obj.avgConst .* obj.avgGradNorm;
    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % --- Save Data for Analysis ---
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                    
                    % Fill analysis arrays
                    obj.saveParams(:,itr)   = obj.param;
                    obj.saveGrads(:,itr)    = obj.grads;
                    obj.saveGradNorm(itr,1) = obj.gradNorm;
                    obj.saveInvTemp(itr,1)  = obj.temper;
                    obj.saveStepSize(:,itr) = obj.stepSize; 
                    obj.saveObjectiveFunc(itr,1) = obj.cpe.J;
    
                    % Update waitbar and message
                    if getappdata(wbar, 'canceling')
                        disp('Simulation cancelled!')
                        break
                    end    
                    if mod(itr, 10) == 0
                        waitbar(itr / obj.numIter);
                    end
                end
    
                % Find the optimum parameter vector
                paramSet = obj.saveParams(:, obj.numIter - obj.numSamplesToUse:end);
                objectiveFuncVals = zeros(obj.numSamplesToUse, 1);
                for nind = 1:obj.numSamplesToUse
                    objectiveFuncVals(nind, 1) = obj.cpe.evalObjectiveFunc(real(paramSet(:, nind)));
                end
                [~, paramOptInd] = min(objectiveFuncVals);
                obj.paramOpt = paramSet(:, paramOptInd);
    
                % Finally, compute the scalar gains
                obj.cpe = obj.cpe.compScalarGains(); % This should now take the optimum param from above
                
                delete(wbar);
                disp('Simulation complete!');

            catch me
                delete(wbar); % Close wait bar when simulation errors out
                rethrow(me);
            end
        end
    end
end
