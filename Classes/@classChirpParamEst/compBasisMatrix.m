function obj = compBasisMatrix(obj)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Nc = obj.Nc;

obj.H = zeros(2*obj.N, Nc);
for c = 1:Nc
    obj.c = c;
    obj = obj.compBasisSignals(); 
    obj.H(:,c) = [real(obj.x) ; imag(obj.x)];
end

end