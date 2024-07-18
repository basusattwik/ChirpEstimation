function obj = convertGradJArray2Cell(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Nc = obj.Nc;
Pc = obj.Pc;

startInd = 1;
for c = 1:Nc
    obj.phiEstCell{1,c} = phiEst(startInd:Pc(c));
    startInd = Pc(c) + 1;
end

end