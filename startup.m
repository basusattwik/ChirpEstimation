% Startup script for the Chirp Estimation simulation framework

close all
clearvars
clc

subfolders = {'Input', 'Output'};
for i = 1:numel(subfolders)
    name = ['Data/', subfolders{i}];  % todo: this is not working
    if ~exist(name, 'dir')
       mkdir(name);
       addpath(genpath(name));
    end
end

% Add folders to MATLAB path
folders = {'Generate', 'Simulations', 'Helpers', 'Examples', 'Experiments', 'Classes', 'Analysis', 'BasicPlots', 'Baselines', 'Data'};
for i = 1:numel(folders)
    addpath(genpath(folders{i}));
end

disp("ChirpEstimation: Ready to run!");