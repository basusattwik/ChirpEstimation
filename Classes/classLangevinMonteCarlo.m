classdef classLangevinMonteCarlo
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        
        % Hyperparameters
        stepSize;
        numIter; 
        lambda;

        % Parameters 
        theta;
        dtheta;
        costJ; 
        noise; 

    end

    methods
        function obj = classLangevinMonteCarlo(inputArg1,inputArg2)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end

        function outputArg = runLmcProcess(obj,cpe)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            obj.dtheta = [cpe.dJ_beta; cpe.dJ_gamma, cpe.dJ_phi];

            for itr = 1:numItr
                obj.theta  = obj.theta - obj.stepSize * obj.dtheta + sqrt(obj.stepSize / obj.lambda) * randn(size(obj.dtheta));
            end
        end
    end
end
