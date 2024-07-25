classdef classLangevinMonteCarlo < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        
        % Hyperparameters
        stepSize;
        numIter; 
        temper;
        tempSwap;

        % Parameters 
        numParams;
        theta;
        dtheta;
        objectiveJ; 
        noise; 
        paramOpt; % Optimum/Almost-optimum parameters
        numSamplesToUse;

        % External
        cpe; % classChirpParamEst

        % For analysis
        objectiveFunc;
        params;
        grads;
        gradNorm;
        invTemp;
        bAccept;
    end

    methods
        function obj = classLangevinMonteCarlo(tuning, setup)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.numParams = setup.numParams;

            % Set tuning for hyperparams
            obj.stepSize = tuning.stepSize;
            obj.numIter  = tuning.numIter;
            obj.temper   = tuning.initInvTemp;
            obj.tempSwap = tuning.tempSwap;
            obj.numSamplesToUse = tuning.numSamplesToUse;

            % Analysis buffers
            obj.objectiveFunc = zeros(obj.numIter, 1);
            obj.params   = zeros(obj.numParams, obj.numIter);
            obj.grads    = zeros(obj.numParams, obj.numIter);
            obj.gradNorm = zeros(obj.numIter, 1);
            obj.invTemp  = zeros(obj.numIter, 1);
            obj.bAccept  = zeros(obj.numIter, 1);

            % Init params and gradient
            obj.theta  = [5, 7, 5, 4, 0, 40, 0, 65].'; %randn(obj.numParams, 1); % what is a good init [0, 0, 20000, 20000, 0, 40, 0, 65].'
            obj.dtheta = zeros(obj.numParams, 1); 

            % Call constructor of cpe class
            obj.cpe = classChirpParamEst(setup);

        end

        function obj = runLmcCore(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
           

            % Create a wait bar display
            wbar = waitbar(0, 'Please wait', 'Name','Running Langevin Monte Carlo', ...
                              'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');
            for itr = 1:obj.numIter

                % parfor np = 1:numParticles
                % end
                % Call the process function for CPE
                obj.cpe = obj.cpe.runCpeCore(obj.theta); % This gives the gradients wrt all params

                % Do Langevin updates on all params
                obj.dtheta = [obj.cpe.dJ_beta.'; obj.cpe.dJ_gamma.'; obj.cpe.dJ_phi.'];               
                thetaProp  = obj.theta - obj.stepSize .* obj.dtheta + sqrt(2 * obj.stepSize / obj.temper) .* randn(obj.numParams, 1);     
               
                % Run Metropolis Adjustment
                alpha = min(1, exp(-obj.temper * (obj.cpe.evalObjectiveFunc(thetaProp) - obj.cpe.J)));
                if rand(1,1) <= alpha
                    obj.theta = thetaProp;
                    obj.bAccept(itr,1) = 1;
                end

                % Run Simulated Tempering
                if mod(itr, obj.tempSwap) == 0
                    obj.temper = log(itr) / 0.9; %beta + exp(l * 0.000001);
                end
                
                % Do Stepsize Annealing
                


                % Fill analysis arrays
                obj.objectiveFunc(itr,1) = obj.cpe.J;
                obj.params(:,itr)   = obj.theta;
                obj.grads(:,itr)    = obj.dtheta;
                obj.gradNorm(itr,1) = norm(obj.dtheta);
                obj.invTemp(itr,1)  = obj.temper;

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
            paramSet = obj.params(:, obj.numIter - obj.numSamplesToUse:end);
            objectiveFuncVals = zeros(obj.numSamplesToUse, 1);
            for nind = 1:obj.numSamplesToUse
                objectiveFuncVals(nind, 1) = obj.cpe.evalObjectiveFunc(real(paramSet(:, nind)));
            end
            [~, paramOptInd] = min(objectiveFuncVals);
            obj.paramOpt = paramSet(:, paramOptInd);

            % Finally, compute the scalar gains
            obj.cpe = obj.cpe.compScalarGains(); % This should now take the optimum param from above
            
            % Insert analysis buffers for storing data
            % ... ToDo ...

            delete(wbar);
            disp('Simulation complete!');
        end
    end
end
