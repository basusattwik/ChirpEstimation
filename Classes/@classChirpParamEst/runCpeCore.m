function obj = runCpeCore(obj, params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Reformat arrays to cells
obj = obj.convertParamArray2Cell(params);

% Main process
obj = obj.compBasisMatrix();
obj = obj.compProjMatrix();
obj = obj.compAllGradients();

end