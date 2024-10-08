classdef classChirpParamEst < handle
    %CLASSCHIRPPARAMEST Summary of this class goes here
    %   Detailed explanation goes here

    properties

        % System properties
        fs = 1000; % Sampling rate (Hz)
        Td = 1;    % Duration (sec)
        Nc = 2;    % Number of Chirps
        N  = 1000; % Total number of samples (fs * Td)

        % Tuning
        gamma;

        % Indexing
        n;        % Sample indexing vector from 0:N-1
        c;        % Chirp index
        p;        % Polynomial phase coeff index
        k;        % Index of the current polynomial phase coeff (in the entire set)
        Kp;       % Total number of polynomial phase coeff
        Ka;       % Total number of polynomial amplitude coeff
        Pc;       % Array containing number of polynomial phase coeff in each chirp
        Ac;       % Array containing number of polynomial amplitude coeff in each chirp

        % Signals for iteration steps
        a;        % Chirp amplitude envelope (N x Nc)
        e;        % Chirps without amplitude envelope (N x Ka)
        u;        % Chirps with amplitude envelope (N x Nc)
        x;        % Clean chirp (N x 1)
        xg;       % Basis for gamma (Nx1)
        w;        % White Gaussian noise at specified snr (N x 1)

        % Generated signals
        am;       % Chirp amplitude envelopes (N x Nc)
        em;       % Chirps without amplitude envelope (N x Nc)
        um;       % Chirps with amplitude envelepe (N x Nc)
        xm;       % Clean mixture of chirps (N x 1)
        ym;       % Noisy mixture of chirps (N x 1)
        wm;       % White Gaussian noise at specified snr (N x 1)
        fim;      % Instantaneous frequency

        % Generated signals
        amRecon;       % Chirp amplitude envelopes (N x Nc)
        emRecon;       % Chirps without amplitude envelope (N x Nc)
        umRecon;       % Chirps with amplitude envelepe (N x Nc)
        xmRecon;       % Clean mixture of chirps (N x 1)
        ymRecon;       % Noisy mixture of chirps (N x 1)
        fimRecon;      % Instantaneous frequency

        % Chirp properties (settings)
        phi;      % Phase polynomial parameters (1 x Nc but cell array)
        rho;      % Amplitude polynomial parameters (1 x Nc but cell array)
        snr;      % Signal-to-Noise ratio in dB

        % Terms used in iterative estimation
        J;          % Objective function value (scalar)
        rhoEst;     
        rhoEstCell;
        phiEst;     % Phase polynomial parameters array (1 x K)
        phiEstCell; % Phase polynomial parameters in cell (1 x Nc)
        phi0Est;
        H;          % Basis matrix (N x Ka)
        Hhat;       % An important intermediate term (N x Nc)
        P;          % Signal projection matrix
        Po;         % Noise projection matrix
        Idn;         % Identity matrix used to compute Po = I - P;
        Idk;
        bvec;       % Vector of amplitude & phase offsets

        % Gradients
        dH_phi;   % Gradient of H wrt phi   (N x Nc x K)
        dJ_phi;   % Gradient of J wrt phi   (1 x K)
        d2J_phi;  % Hessian of J wrt phi    (K x K)

        % Booleans
        bMinFound = false;
        bApplyWin = false;

        % Tolerances
        minObjTol;

        % Errors
        sqrPhiError;
        sqrRhoError;
        crbnd;      % Cramer-Rao Bound
        fimat;      % Fisher Information Matrix
        sigma2;      % Noise var
        
    end

    methods
        function obj = classChirpParamEst(cpeSetting)
            %CLASSCHIRPPARAMEST Construct an instance of this class
            %   Detailed explanation goes here
            
            % Load chirp param file along with tuning
            obj.fs    = cpeSetting.fs;
            obj.Td    = cpeSetting.Td;
            obj.Nc    = cpeSetting.Nc;
            obj.phi   = cpeSetting.phi;
            obj.rho   = cpeSetting.rho;
            obj.snr   = cpeSetting.snr;
            obj.minObjTol = cpeSetting.minObjTol;
            obj.gamma     = cpeSetting.gamma;
            obj.bApplyWin = cpeSetting.bApplyWin;
            
            obj.N = obj.fs * obj.Td;
            obj.n = (0:obj.N-1).';
            obj.c = 1;
            obj.p = 1;
            obj.k = 1;

            % Fill out Pc and Ac arrays and compute total phase params K
            obj.Pc  = zeros(obj.Nc, 1);
            for c = 1:obj.Nc
                obj.Pc(c,1) = size(obj.phi{1,c}, 1);
            end
            obj.Kp = sum(obj.Pc);

            obj.Ac  = zeros(obj.Nc, 1);
            for c = 1:obj.Nc
                obj.Ac(c,1) = size(obj.rho{1,c}, 1);
            end
            obj.Ka = sum(obj.Ac);

            % Init chirp signals based on the settings file
            obj = obj.resetArrays();
            obj = obj.initChirpSignals();
            obj = compCramerRaoBounds(obj);
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
        obj = compScalarGains(obj, params);
        obj = reconChirpSignals(obj)

        % The process function
        obj = runCpeCore(obj, params);

        % Error Analysis
        obj = evalParamErrors(obj);

    end
end