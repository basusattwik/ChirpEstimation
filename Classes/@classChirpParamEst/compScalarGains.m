function obj = compScalarGains(obj, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Cache for speed
fs = obj.fs;
N  = obj.N;
Nc = obj.Nc;
Pc = obj.Pc;
ym = obj.ym;

% Preallocate for speed
phi   = cell(1, Nc);
e     = zeros(N, 1);
H     = zeros(N, Nc);

% Avoiding divides
oneOverFs = 1/fs;

% Extract current params into data structs
startInd = 1;
for c = 1:Nc
    phi{1,c} = params(startInd:startInd+Pc(c)-1);
    startInd = startInd + Pc(c);
end

% Compute basis matrix
for c = 1:Nc
    P = size(phi{1,c},1);
    pvec = (0:P-1).';
    for n = 1:N % -- loop over number of samples

        % Get the exponential polynomial phase sinusoid
        npvec  = ((n-1) * oneOverFs).^pvec; % vectors of powers of n/fs
        e(n,1) = exp(2*pi*1j .* (phi{1,c}.' * npvec));
    end
    x = e;
    H(:,c) = x;
end

% Get the projection matrix and the orthogonal projection matrix
obj.alpha  = (H' * H)^(-1) * H' * obj.ym; 

end