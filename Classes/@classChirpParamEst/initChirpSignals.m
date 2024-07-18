function obj = initChirpSignals(obj)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

obj = genChirpSignal(obj);   % Generate the chirp according to the settings file
obj = addGaussianNoise(obj); % Add noise and store signals in class

% Fill out Kc array
for c = 1:obj.Nc
    obj.Pc(c,1) = size(obj.phi{1,c}, 1);
end

obj.K = sum(obj.Pc);

end

