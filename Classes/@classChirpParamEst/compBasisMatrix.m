function obj = compBasisMatrix(obj)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Nc = obj.Nc;

obj.H = zeros(2*obj.N, Nc);
for c = 1:Nc
    yc = obj.compBasisSignals(); % Nc and alpha arguments are set to 1 here. Make this modification
    obj.H(:,c) = [real(yc) ; imag(yc)];
end

end