function obj = compBasisSignals(obj)
%GENCHIRPSIGNAL Summary of this function goes here
%   Detailed explanation goes here

% Cache some variables for speed
N  = obj.N;
fs = obj.fs;
c  = obj.c;

phiEstCell = obj.phiEstCell;

% max polynomial degree
P = size(phiEstCell{1,c},1);
pvec = (0:P-1).';

% Avoiding divides
oneOverFs = 1/fs;

for n = 1:N % -- loop over number of samples

    % Get the exponential polynomial phase sinusoid
    npvec = ((n-1) * oneOverFs).^pvec; % vectors of powers of n/fs
    obj.e(n,1) = exp(2*pi*1j .* (phiEstCell{1,c}.' * npvec));

end

obj.x = obj.e;

end


