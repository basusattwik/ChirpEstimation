function obj = compBasisMatrix(obj)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Nc = obj.Nc;

% Fill columns with the basis signal 
for c = 1:Nc
    obj.c = c;
    obj = obj.compBasisSignals(); 
    obj.H(:,c) = obj.x;
end

end