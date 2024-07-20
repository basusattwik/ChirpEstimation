function obj = convertParamArray2Cell(obj, params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Nc = obj.Nc;
Pc = obj.Pc;

for c = 1:Nc
    obj.betaEst(1,c)  = params(c,1);
    obj.gammaEst(1,c) = params(c + Nc,1);
end

startInd = 2*Nc + 1;
for c = 1:Nc
    obj.phiEstCell{1,c} = params(startInd:startInd+Pc(c)-1);
    startInd = startInd + Pc(c);
end

end