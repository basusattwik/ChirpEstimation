function J = evalCostFunc(obj, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

fs = obj.fs;
N  = obj.N;
Nc = obj.Nc;
Pc = obj.Pc;

beta  = zeros(1, obj.Nc); 
gamma = zeros(1, obj.Nc);
phi = cell(1, obj.Nc);
A = zeros(obj.N, 1);
e = zeros(obj.N, 1);
H = zeros(2*obj.N, Nc);


for c = 1:Nc
    beta(1,c)  = params(c,1);
    gamma(1,c) = params(c + Nc,1);
end
startInd = 2*Nc + 1;
for c = 1:Nc
    phi{1,c} = params(startInd:startInd+Pc(c)-1);
    startInd = startInd + Pc(c);
end

for c = 1:Nc
    P = size(phi{1,c},1);
    for n = 1:N % -- loop over number of samples
        % Get the amplitude envelope
        A(n,1) = exp(-beta(1,c) * n/fs) * (1 - exp(-gamma(1,c) * n/fs));
    
        % Get the exponential polynomial phase sinusoid
        npvec = ((n-1) / fs).^(0:P-1).'; % vectors of powers of n/fs
        e(n,1) = exp(2*pi*1j .* (phi{1,c}.' * npvec));
    end
    x = A .* e;
    H(:,c) = [real(x) ; imag(x)];
end

% Get the projection matrix and the orthogonal projection matrix
P  = H * pinv(H.' * H) * H.'; 
Po = eye(size(P)) - P;

% Cost function value
J = obj.ymvec.' * Po * obj.ymvec;

end