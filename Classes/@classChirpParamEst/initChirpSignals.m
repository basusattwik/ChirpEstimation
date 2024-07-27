function obj = initChirpSignals(obj)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

obj = genChirpSignal(obj);   % Generate the chirp according to the settings file
obj = addGaussianNoise(obj); % Add noise and store signals in class

end

