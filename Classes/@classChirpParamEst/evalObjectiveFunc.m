function J = evalObjectiveFunc(obj, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Cache for speed
fs = obj.fs;
N  = obj.N;
Nc = obj.Nc;
Pc = obj.Pc;
ym = obj.ym;

% Preallocate for speed
beta  = zeros(1, Nc); 
gamma = zeros(1, Nc);
phi   = cell(1, Nc);
A     = zeros(N, 1);
e     = zeros(N, 1);
H     = zeros(N, Nc);

% Avoiding divides
oneOverFs = 1/fs;

% Extract current params into data structs
for c = 1:Nc
    beta(1,c)  = params(c,1);
    gamma(1,c) = params(c + Nc,1);
end
startInd = 2*Nc + 1;
for c = 1:Nc
    phi{1,c} = params(startInd:startInd+Pc(c)-1);
    startInd = startInd + Pc(c);
end

% Compute basis matrix
for c = 1:Nc
    P = size(phi{1,c},1);
    pvec = (0:P-1).';
    for n = 1:N % -- loop over number of samples
        % Get the amplitude envelope
        A(n,1) = exp(-beta(1,c) * (n-1) * oneOverFs) * (1 - exp(-gamma(1,c) * (n-1) * oneOverFs));
    
        % Get the exponential polynomial phase sinusoid
        npvec  = ((n-1) * oneOverFs).^pvec; % vectors of powers of n/fs
        e(n,1) = exp(2*pi*1j .* (phi{1,c}.' * npvec));
    end
    x = A .* e;
    H(:,c) = x;
end

% Get the projection matrix and the orthogonal projection matrix
P  = H * (H' * H)^(-1) * H'; 
Po = eye(size(P)) - P;

% Objective function value
J = real(ym' * Po * ym);

end