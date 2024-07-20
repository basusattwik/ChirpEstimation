classdef classLangevinMonteCarlo < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        
        % Hyperparameters
        stepSize;
        numIter; 
        temper;

        % Parameters 
        numParams;
        theta;
        dtheta;
        costJ; 
        noise; 

        cpe;% = classChirpParamEst;
    end

    methods
        function obj = classLangevinMonteCarlo(tuning, setup)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.numParams = setup.numParams;

            % Init params and gradient
            obj.theta  = randn(obj.numParams, 1);
            obj.dtheta = zeros(obj.numParams, 1); 


            % Set tuning for hyperparams
            obj.stepSize = tuning.stepSize;
            obj.numIter  = tuning.numIter;
            obj.temper   = tuning.temper;

            obj.cpe = classChirpParamEst(setup);

        end

        function obj = runLmcCore(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
           

            for itr = 1:obj.numIter

                % Call the process function for CPE
                obj.cpe = obj.cpe.runCpeCore(obj.theta); % This gives the gradients wrt all params

                % Do Langevin updates on all params
                obj.dtheta = [obj.cpe.dJ_beta.'; obj.cpe.dJ_gamma.'; obj.cpe.dJ_phi.'];
                obj.theta  = obj.theta - obj.stepSize * obj.dtheta + ...
                               sqrt(2 * obj.stepSize / obj.temper) * randn(size(obj.dtheta));

            end

            obj.cpe = obj.cpe.compScalarGains();
        end
    end
end
