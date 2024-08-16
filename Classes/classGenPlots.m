classdef classGenPlots < handle
    %CLASSGENPLOTS Summary of this class goes here
    %   Detailed explanation goes here

    properties
        param;
        stepSize;
        grad;
        gradNorm;
        avgGradNorm;
        objFunc;
    end

    methods
        function obj = classGenPlots(simData, setup, tuning)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.param       = simData.saveParams;
            obj.grad        = simData.saveGrads;
            obj.stepSize    = simData.saveStepSize;
            obj.gradNorm    = simData.saveGradNorm;
            obj.avgGradNorm = simData.saveAvgGradNorm;
            obj.objFunc     = simData.saveObjFunc;
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end