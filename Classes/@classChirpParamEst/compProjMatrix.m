function obj = compProjMatrix(obj)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Compute Hhat
obj.Hhat = (obj.H' * obj.H)^(-1) * obj.H'; 

% Get the projection matrix and the orthogonal projection matrix
obj.P  = obj.H * obj.Hhat;
obj.Po = eye(obj.N, obj.N) - obj.P;

end