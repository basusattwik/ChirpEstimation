function obj = compCostFunc(obj)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Cost function value
obj.J = obj.y.' * obj.P * obj.y;

end