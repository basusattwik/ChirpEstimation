function obj = convertParamArray2Cell(obj, params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Nc = obj.Nc;
Pc = obj.Pc;

startInd = 1;
for c = 1:Nc
    obj.phiEstCell{1,c} = params(startInd:startInd+Pc(c)-1);
    startInd = startInd + Pc(c);
end

end