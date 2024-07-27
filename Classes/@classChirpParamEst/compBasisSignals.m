function obj = compBasisSignals(obj)
%GENCHIRPSIGNAL Summary of this function goes here
%   Detailed explanation goes here

% Cache some variables for speed
N  = obj.N;
fs = obj.fs;
c  = obj.c;

phiEstCell = obj.phiEstCell;
betaEst    = obj.betaEst;
gammaEst   = obj.gammaEst;

% max polynomial degree
P = size(phiEstCell{1,c},1);
pvec = (0:P-1).';

% Avoiding divides
oneOverFs = 1/fs;

if ~obj.bAmpGamma
    for n = 1:N % -- loop over number of samples
    
        % Get the amplitude envelope
        obj.A(n,1) = exp(-betaEst(1,c) * (n-1) * oneOverFs) * (1 - exp(-gammaEst(1,c) * (n-1)  * oneOverFs));
    
        % Get the exponential polynomial phase sinusoid
        npvec = ((n-1) * oneOverFs).^pvec; % vectors of powers of n/fs
        obj.e(n,1) = exp(2*pi*1j .* (phiEstCell{1,c}.' * npvec));
    
    end
    % Multiply the two
    if obj.bAmpEnv
        obj.x = obj.A .* obj.e;
    else
        obj.x = obj.e;
    end
else
    for n = 1:N % -- loop over number of samples
    
        % Get the amplitude envelope
        obj.A(n,1) = exp(-(betaEst(1,c) + gammaEst(1,c)) * (n-1) * oneOverFs); % Unnecessary sum of beta and gamma. Optimize. 
    
        % Get the exponential polynomial phase sinusoid
        npvec = ((n-1) * oneOverFs).^pvec; % vectors of powers of n/fs
        obj.e(n,1) = exp(2*pi*1j .* (phiEstCell{1,c}.' * npvec));
    
    end
    % Multiply the two
    if obj.bAmpEnv
        obj.xg = obj.A .* obj.e;
    else
        obj.xg = obj.e;
    end
end

end


