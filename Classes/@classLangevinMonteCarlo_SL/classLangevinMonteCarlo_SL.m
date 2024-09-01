classdef classLangevinMonteCarlo_SL < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties

        % Hyperparameters
        stepSize;
        stepSizeInit;
        stepSizeConst;
        stepSizeMax;
        stepSizeMin;
        initTemp;
        tempConst;
        avgConst;
        numParticles;
        numIterLmc;
        numIterSmooth;
        noiseVarInit; % should be high
        noiseVarMin;
        stepNoiseVar;
        initValMinMax;

        % Parameters 
        optParam;
        numParams;
        bestParticleInd;

        % States
        temper;
        param;
        initParam;
        grads;
        hess; 
        objectiveJ; 
        gradNorm;
        avgGradNorm;
        avgGrads; 
        noiseVar;
        absHessianTrace;

        % Constant Matrices
        Id; % Identity matrix

        % Bools
        bAccept;
        bDisplayPlots;
        bStopLmc;
        bEnableGaussSmooth;
        bEnableLangevin;
        bInitParamGiven;
        bStopSim;

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
        saveNoiseVar;
        saveProposedParams;
        saveHessTrc;
        
    end

    methods
        function obj = classLangevinMonteCarlo_SL(tuning, setup)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.numParams = setup.numParams;
            obj.bStopLmc  = false;
            
            % Set tuning for hyperparams
            obj.stepSizeConst = tuning.stepSizeConst;
            obj.tempConst     = tuning.tempConst;
            obj.avgConst      = tuning.avgConst;
            obj.numParticles  = tuning.numParticles;
            obj.numIterLmc    = tuning.numIterLmc;
            obj.noiseVarInit  = tuning.noiseVarInit;
            obj.noiseVarMin   = tuning.noiseVarMin;
            obj.bDisplayPlots = tuning.bDisplayPlots;
            obj.numIterSmooth = tuning.numIterSmooth;
            obj.bEnableGaussSmooth = tuning.bGaussSmooth;
            obj.bEnableLangevin = tuning.bEnableLangevin;
            obj.initValMinMax = tuning.initValMinMax;

            % Init tuning param arrays
            obj.stepSizeInit = zeros(obj.numParams, obj.numParticles);
            obj.stepSize     = zeros(obj.numParams, obj.numParticles);
            obj.stepSizeMax  = zeros(obj.numParams, obj.numParticles);
            obj.temper       = zeros(1, 1);
            obj.noiseVar     = zeros(1, obj.numParticles);
            obj.bestParticleInd  = [];
            
            for pind = 1:obj.numParticles
                obj.stepSizeInit(:,pind) = tuning.stepSize(1:obj.numParams);
                obj.stepNoiseVar(:,pind) = tuning.stepNoiseVar;
                obj.noiseVar(1,pind) = tuning.noiseVarInit;
                % for npr = 1:obj.numParams
                %     obj.stepSizeMax(npr,pind) = tuning.stepSizeMax;
                % end
            end

            obj.stepSizeMax = tuning.stepSizeMax;
            obj.stepSizeMin = tuning.stepSizeMin;

            % Init state params and gradients
            if ~isempty(tuning.initParams)
                obj.param = tuning.initParams;
            else
                obj.param = unifrnd(obj.initValMinMax(1), obj.initValMinMax(2), obj.numParams, obj.numParticles); %[30, 57, 59.3];
            end

            % obj.param       = normrnd(obj.initValMinMax(1), obj.initValMinMax(2), obj.numParams, obj.numParticles);
            obj.grads       = zeros(obj.numParams, obj.numParticles);
            obj.hess        = zeros(obj.numParams, obj.numParams, obj.numParticles);
            obj.avgGrads    = zeros(obj.numParams, obj.numParticles);
            obj.avgGradNorm = zeros(1, obj.numParticles);
            obj.gradNorm    = zeros(1, obj.numParticles);
            obj.Id          = eye(obj.numParams, obj.numParams);
            obj.absHessianTrace = zeros(1, obj.numParticles);

            obj.initParam = obj.param; % just store this for use during analysis
            obj.bAccept   = false(1,obj.numParticles);

            % Preallocate arrays to store state data for later analysis
            obj.saveParams      = zeros(obj.numParams,    obj.numParticles, obj.numIterLmc);
            obj.saveGrads       = zeros(obj.numParams,    obj.numParticles, obj.numIterLmc);
            obj.saveStepSize    = zeros(obj.numParams,    obj.numParticles, obj.numIterLmc);
            obj.saveGradNorm    = zeros(obj.numParticles, obj.numIterLmc);
            obj.saveAvgGradNorm = zeros(obj.numParticles, obj.numIterLmc);
            obj.savebAccept     = zeros(obj.numParticles, obj.numIterLmc);
            obj.saveObjFunc     = zeros(obj.numParticles, obj.numIterLmc);
            obj.saveTemper      = zeros(obj.numIterLmc, 1);
            obj.saveNoiseVar    = zeros(obj.numParticles, obj.numIterLmc);
            obj.saveProposedParams = zeros(obj.numParams,    obj.numParticles, obj.numIterLmc);
            obj.saveHessTrc     = zeros(obj.numParticles, obj.numIterLmc);

            obj.bStopSim = false;
            %  Call constructor of cpe class. Initialize several instances
            %  of the class... one for each particle
            obj.cpe = cell(1,obj.numParticles);
            for pind = 1:obj.numParticles
                obj.cpe{1,pind} = classChirpParamEst(setup);
            end
        end

        % Function declarations are provided below. Definitions are in
        % separate .m files in the class folder. 
        obj = runLmcCore_SL(obj);
        ax  = initLivePlots_SL(obj);
        []  = updateLivePlots_SL(obj, tind, ax);

    end
end
