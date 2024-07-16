classdef classChirpParamEst
    %CLASSCHIRPPARAMEST Summary of this class goes here
    %   Detailed explanation goes here

    properties

        % System properties
        fs; % Sampling rate (Hz)
        Td; % Duration (sec)
        Nc; % Number of Chirps
        N;  % Total number of samples (fs x Td)
        x;  % Clean mixture of chirps;
        y;  % Noisy mixture of chirps
        w;  % White Gaussian noise at specified snr

        % Chirp properties
        alpha; % Ssalar gains
        beta;  % Amplitude envelope parameter
        gamma; % The other amplitude envelope paramater
        phi;   % Phase parameter
        snr;

        % Estimation 
        J;        % Cost function value
        H;        % Basis matrix
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
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            
            % Load chirp param file along with tuning
            %...
        end

        % Function declaration is provided below. Definitions are in
        % separate .m files. 

        obj = genChirpSignal(obj);
        obj = compCostFunc(obj);
        obj = compScalarGains(obj);
        obj = compBasisMatrix(obj);
        obj = compProjMatrix(obj);
        obj = compAllGrads(obj); % tough one
        obj = compIndividualGrads(obj); % very tough one

    end
end