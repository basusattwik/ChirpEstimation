function obj = runEstProcess(obj)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Init chirp signals based on the settings file
obj = initChirpSignals(obj);
obj = resetArrays(obj);

% These things should go into a loop inside LMC process()

% Get the basis matrix H
obj = compBasisMatrix(obj);

% Get the projection matrix and Hhat
obj = compProjMatrix(obj);

% Get gradients
obj = compGradsForH(obj);
obj = compGradsforJ(obj);

% Run Langevin Monte Carlo with Simulated Tempering
% obj = obj.lmc.runLmcUpdates();


end