function obj = compObjectiveFunc(obj)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Cost function value
obj.J = real(obj.ym' * obj.Po * obj.ym); % Force it to be real to prevent tiny imaginary values ~ e-16

end