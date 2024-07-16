function obj = compProjMatrix(obj)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Compute Hhat
obj.Hhat = pinv(obj.H.' * obj.H) * obj.H.'; 

% Get the projection matrix and the orthogonal projection matrix
obj.P  = obj.H * obj.Hhat;
obj.Po = eye(size(obj.P)) - obj.P;

end