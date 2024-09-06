function obj = compScalarGains(obj, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Cache for speed
fs = obj.fs;
N  = obj.N;
Nc = obj.Nc;
Pc = obj.Pc;
Ac = obj.Ac;
ym = obj.ym;
Ka = obj.Ka;

% Preallocate for speed
phiCell = cell(1, Nc);
H = zeros(N, Ka);

nvecOverFs = obj.n / fs;
twoPij = 2*pi*1j;

startInd = 1;
for c = 1:Nc
    phiCell{1,c} = params(startInd:startInd+Pc(c)-2);
    startInd = startInd + Pc(c) - 1;
end

startInd = 1;
for c = 1:Nc
    phi  = phiCell{1,c};
    pvec = (1:Pc(c)-1);
    avec = (0:Ac(c)-1);

    endInd = startInd+Ac(c)-1;
    H(:, startInd:endInd) = (nvecOverFs.^avec) .* exp(twoPij .* ((nvecOverFs.^pvec) * phi));
    startInd = startInd + Ac(c);
end

% Get the projection matrix and the orthogonal projection matrix
obj.bvec  = ((H' * H + obj.gamma * obj.Idk) \ H') * ym; 

% Get absolute vals (rho)
obj.rhoEst = abs(obj.bvec);

% Get phase (phi_0)
phi0 = angle(obj.bvec);

startInd = 1;
for c = 1:Nc
    endInd   = startInd + Ac(c) - 1;
    obj.rhoEstCell{1,c} = obj.rhoEst(startInd:endInd);
    obj.phi0Est(1,c) = mean(phi0(startInd:endInd));
    startInd = startInd + Ac(c);
end

end