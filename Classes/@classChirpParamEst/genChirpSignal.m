function obj = genChirpSignal(obj)
%GENCHIRPSIGNAL Summary of this function goes here
%   Detailed explanation goes here

Nc = obj.Nc;
N  = obj.N;
fs = obj.fs;

% Calculate the polynomial chirp
for c = 1:Nc  % -- loop over number of chirps
    % max polynomial degree
    P = size(obj.phi{1,c},1);

    for nind = 1:N % -- loop over number of samples

        % Get the amplitude envelope
        obj.Am(nind,c) = obj.alpha(1,c) * exp(-obj.beta(1,c) * (nind-1)/fs) * (1 - exp(-obj.gamma(1,c) * (nind-1)/fs));

        % Get the exponential polynomial phase sinusoid
        npvec = ((nind-1) / fs).^(0:P-1).'; % vectors of powers of n/fs
        obj.em(nind,c) = exp(2*pi*1j .* (obj.phi{:,c}.' * npvec));

    end % -- end loop over number of chirps
end % -- end loop over number of samples

if obj.bAmpEnv
    obj.xm = obj.Am .* obj.em;
else
    obj.xm = obj.em;
end

% Combined to form multicomponent signal
obj.ym = sum(obj.xm, 2);

end


