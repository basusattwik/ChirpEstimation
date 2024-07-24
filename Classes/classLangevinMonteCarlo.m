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
        costJ; 
        noise; 

        % External
        cpe; % classChirpParamEst

        % For analysis
        costFunc;
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

            % Analysis buffers
            obj.costFunc = zeros(obj.numIter, 1);
            obj.params   = zeros(obj.numParams, obj.numIter);
            obj.grads    = zeros(obj.numParams, obj.numIter);
            obj.gradNorm = zeros(obj.numIter, 1);
            obj.invTemp  = zeros(obj.numIter, 1);
            obj.bAccept  = zeros(obj.numIter, 1);

            % Init params and gradient
            obj.theta  = randn(obj.numParams, 1); % what is a good init [0, 0, 20000, 20000, 0, 40, 0, 65].'
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

                % Call the process function for CPE
                obj.cpe = obj.cpe.runCpeCore(obj.theta); % This gives the gradients wrt all params

                % Do Langevin updates on all params
                obj.dtheta = [obj.cpe.dJ_beta.'; obj.cpe.dJ_gamma.'; obj.cpe.dJ_phi.'];               
                thetaProp  = obj.theta - obj.stepSize .* obj.dtheta + sqrt(2 * obj.stepSize / obj.temper) .* randn(obj.numParams, 1);     
               
                % Run Metropolis Adjustment
                alpha = min(1, exp(-obj.temper * (evalCostFunc(obj.cpe, thetaProp) - obj.cpe.J)));
                if rand(1,1) <= alpha
                    obj.theta = thetaProp;
                    obj.bAccept(itr,1) = 1;
                end

                % Run Simulated Tempering
                if mod(itr, obj.tempSwap) == 0
                    obj.temper = log(itr) / 0.2; %beta + exp(l * 0.000001);
                end
                
                % Do Stepsize Annealing
                


                % Fill analysis arrays
                obj.costFunc(itr,1) = obj.cpe.J;
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

            % Finally, compute the scalar gains
            obj.cpe = obj.cpe.compScalarGains();
            
            % Insert analysis buffers for storing data
            % ... ToDo ...

            close(wbar);
            delete(wbar);
            disp('Simulation complete!');
        end
    end
end
