function obj = runCpeCore(obj, params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Reformat arrays to cells
obj = obj.convertParamArray2Cell(params);

% Update the basis matrix
obj = obj.compBasisMatrix();

% Update the projection matrix
obj = obj.compProjMatrix();

% Compute new gradients based on the new parameters
obj = obj.compAllGradients();

% Compute the new cost function value
obj = obj.compObjectiveFunc();

end

% General comments:
% compBasisSignals() is being run more times than necessary. Can be done
% once and used in Basis, Project and Grads. But need to be careful about
% the gamma part in grads.