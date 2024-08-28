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
        numIterSmooth;
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
        bEnableGaussSmooth;
        bEnableLangevin;

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
        saveProposedParams;
    end

    methods
        function obj = classLangevinMonteCarlo(tuning, setup)
            %CLASSLANGEVINMONTECARLO Construct an instance of this class
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
            obj.numIterSmooth = tuning.numIterSmooth;
            obj.bEnableGaussSmooth  = tuning.bGaussSmooth;
            obj.bEnableLangevin = tuning.bEnableLangevin;

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
            obj.param       = unifrnd(-100, 100, obj.numParams, obj.numParticles);%
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
            obj.saveProposedParams = zeros(obj.numParams, obj.numParticles, obj.numIterLmc, obj.numIterNoise);

            %  Call constructor of cpe class. Initialize several instances
            %  of the class... one for each particle
            obj.cpe = cell(1,obj.numParticles);
            for pind = 1:obj.numParticles
                obj.cpe{1,pind} = classChirpParamEst(setup);
            end
        end

        % Function declarations are provided below. Definitions are in
        % separate .m files in the class folder. 
        obj = runLmcCore(obj);
        obj = runLmcCore_singleLoop(obj);
        ax  = initLivePlots(obj);
        []  = updateLivePlots(obj, tind, nind, ax);
       
    end
end
