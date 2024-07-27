classdef classChirpParamEst < handle
    %CLASSCHIRPPARAMEST Summary of this class goes here
    %   Detailed explanation goes here

    properties

        % System properties
        fs = 1000; % Sampling rate (Hz)
        Td = 1;    % Duration (sec)
        Nc = 2;    % Number of Chirps
        N  = 1000; % Total number of samples (fs * Td)
        bAmpGamma = false;
        bAmpEnv   = true;

        % Indexing
        n;        % Sample indexing vector from 0:N-1
        c;        % Chirp index
        p;        % Polynomial phase coeff index
        k;        % Index of the current polynomial phase coeff (in the entire set)
        K;        % Total number of polynomial phase coeff
        Pc;       % Array containing number of polynomial phase coeff in each chirp

        % Signals for iteration steps
        A;        % Chirp amplitude envelopes (N x Nc);
        e;        % Chirps without amplitude envelope (N x Nc)
        u;        % Chirps with amplitude envelepe (N x Nc)
        x;        % Clean chirp (N x 1)
        xg;       % Basis for gamma (Nx1)
        w;        % White Gaussian noise at specified snr (N x 1)

        % Generated signals
        Am;       % Chirp amplitude envelopes (N x Nc);
        em;       % Chirps without amplitude envelope (N x Nc)
        um;       % Chirps with amplitude envelepe (N x Nc)
        xm;       % Clean mixture of chirps (N x 1)
        ym;       % Noisy mixture of chirps (N x 1)
        wm;       % White Gaussian noise at specified snr (N x 1)

        % Chirp properties (settings)
        env;      % Amplitude env for each chirp (N x Nc)
        alpha;    % Scalar gains (1 x Nc)
        beta;     % Amplitude envelope parameter (1 x Nc)
        gamma;    % The other amplitude envelope parameter (1 x Nc)
        phi;      % Phase polynomial parameters (1 x Nc but cell array)
        snr;      % Signal-to-Noise ratio in dB

        % Terms used in iterative estimation
        J;        % Objective function value (scalar)
        alphaEst; % Scalar gains (1 x Nc)
        betaEst;  % Amplitude envelope parameter (1 x Nc)
        gammaEst; % The other amplitude envelope parameter (1 x Nc)
        phiEst;   % Phase polynomial parameters array (1 x K)
        phiEstCell; % Phase polynomial parameters in cell (1 x Nc)
        H;        % Basis matrix (N x Nc)
        Hhat;     % An important intermediate term (N x Nc)
        P;        % Signal projection matrix
        Po;       % Noise projection matrix

        % Gradients
        dH_phi;   % Gradient of H wrt phi   (N x Nc x K)
        dH_beta;  % Gradient of H wrt beta  (N x Nc)
        dH_gamma; % Gradient of H wrt gamma (N x Nc)
        dJ_phi;   % Gradient of J wrt phi   (1 x K)
        dJ_beta;  % Gradient of J wrt beta  (1 x Nc)
        dJ_gamma; % Gradient of J wrt gamma (1 x Nc)
        
    end

    methods
        function obj = classChirpParamEst(cpeSetting)
            %CLASSCHIRPPARAMEST Construct an instance of this class
            %   Detailed explanation goes here
            
            % Load chirp param file along with tuning
            obj.fs    = cpeSetting.fs;
            obj.Td    = cpeSetting.Td;
            obj.Nc    = cpeSetting.Nc;
            obj.alpha = cpeSetting.alpha;
            obj.beta  = cpeSetting.beta;
            obj.gamma = cpeSetting.gamma;
            obj.phi   = cpeSetting.phi;
            obj.snr   = cpeSetting.snr;

            % Error checking basics
            if size(obj.phi,2) ~= obj.Nc || size(obj.alpha,2) ~= obj.Nc || size(obj.beta,2) ~= obj.Nc || size(obj.gamma,2) ~= obj.Nc
                error('Number of chirps does not match number of parameters provided');
            end
            
            obj.N = obj.fs * obj.Td;
            obj.n = (0:obj.N-1).';
            obj.c = 1;
            obj.p = 1;
            obj.k = 1;

            % Fill out Pc array and compute total phase params K
            obj.Pc  = zeros(obj.Nc, 1);
            for c = 1:obj.Nc
                obj.Pc(c,1) = size(obj.phi{1,c}, 1);
            end
            obj.K = sum(obj.Pc);

            % Init chirp signals based on the settings file
            obj = obj.resetArrays();
            obj = obj.initChirpSignals();
        end

        % Function declarations are provided below. Definitions are in
        % separate .m files in the class folder. 

        % Initialization
        obj = genChirpSignal(obj);
        obj = addGaussianNoise(obj);
        obj = initChirpSignals(obj);
        obj = resetArrays(obj);

        % Computational steps
        obj = compObjectiveFunc(obj);
        obj = evalObjectiveFunc(obj, params);
        obj = compBasisSignals(obj);
        obj = compBasisMatrix(obj);
        obj = compProjMatrix(obj);
        obj = compAllGradients(obj);
        obj = compScalarGains(obj);

        % The process function
        obj = runCpeCore(obj, params);

        % Helper functions
        obj = convertGradJArray2Cell(obj);

    end
end