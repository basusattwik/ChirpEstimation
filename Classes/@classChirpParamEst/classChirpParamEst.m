classdef classChirpParamEst
    %CLASSCHIRPPARAMEST Summary of this class goes here
    %   Detailed explanation goes here

    properties

        % System properties
        fs = 1000; % Sampling rate (Hz)
        Td = 1;    % Duration (sec)
        Nc = 2;    % Number of Chirps
        N  = 1000; % Total number of samples (fs * Td)
        mode = 'gen'; % Generating chirps or computing

        % Signals for iteration steps
        A;  % Chirp amplitude envelopes (N x Nc);
        e;  % Chirps without amplitude envelope (N x Nc)
        u;  % Chirps with amplitude envelepe (N x Nc)
        x;  % Clean chirp (N x 1)
        w;  % White Gaussian noise at specified snr (N x 1)

        % Generated signals
        Am;  % Chirp amplitude envelopes (N x Nc);
        em;  % Chirps without amplitude envelope (N x Nc)
        um;  % Chirps with amplitude envelepe (N x Nc)
        xm;  % Clean mixture of chirps (N x 1)
        ym;  % Noisy mixture of chirps (N x 1)
        wm;  % White Gaussian noise at specified snr (N x 1)

        % Chirp properties
        env;   % Amplitude env for each chirp (N x Nc)
        alpha; % Scalar gains (1 x Nc)
        beta;  % Amplitude envelope parameter (1 x Nc)
        gamma; % The other amplitude envelope parameter (1 x Nc)
        phi;   % Phase polynomial parameters (1 x Nc but cell array)
        snr;   % Signal-to-Noise ratio in dB

        % Estimation 
        J;        % Cost function value
        H;        % Basis matrix (N x Nc)
        Hhat;     % An important intermediate term
        P;        % Signal projection matrix
        Po;       % Noise projection matrix
        dH_phi;   % Gradient of H wrt phi
        dH_beta;  % Gradient of H wrt beta
        dH_gamma; % Gradient of H wrt gamma
        
        % MCMC algorithm
        lmc = classLangevinMonteCarlo;
    end

    methods
        function obj = classChirpParamEst(filePath)
            %CLASSCHIRPPARAMEST Construct an instance of this class
            %   Detailed explanation goes here
            
            % Load chirp param file along with tuning
            load(filePath);
            
            % At the end
            obj.switchMode('comp');
        end

        % Function declarations are provided below. Definitions are in
        % separate .m files in the class folder. 

        % Init related
        obj = genChirpSignal(obj);
        obj = addGaussianNoise(obj);
        obj = initChirpSignals(obj);
        obj = switchMode(obj, select);

        % Computational steps
        obj = compCostFunc(obj);
        obj = compScalarGains(obj);
        obj = compBasisMatrix(obj);
        obj = compProjMatrix(obj);
        obj = compAllGrads(obj); % tough one
        obj = compIndividualGrads(obj); % very tough one

        % The process function
        obj = runEstProcess(obj);

    end
end