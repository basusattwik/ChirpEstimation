function obj = compBasisMatrix(obj)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Nc = obj.Nc;

if size(obj.phi,2) ~= Nc || size(obj.alpha,2) ~= Nc || size(obj.beta,2) ~= Nc || size(obj.gamma,2) ~= Nc
    error('Number of chirps does not match number of parameters provided'); % This should move outside
end

obj.H = zeros(2*obj.N, Nc);
for c = 1:Nc
    yc = obj.genChirpSignal(obj); % Nc and alpha arguments are set to 1 here. Make this modification
    obj.H(:,c) = [real(yc) ; imag(yc)];
end

end