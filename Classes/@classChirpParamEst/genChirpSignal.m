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

        % Get the exponential polynomial phase sinusoid
        npvec = ((nind-1) / fs).^(0:P-1).'; % vectors of powers of n/fs
        obj.em(nind,c) = exp(2*pi*1j .* (obj.phi{:,c}.' * npvec));

    end % -- end loop over number of chirps
end % -- end loop over number of samples

obj.xm = obj.em;

% Combined to form multicomponent signal
obj.ym = sum(obj.xm, 2);

end


