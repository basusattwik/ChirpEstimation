function obj = compScalarGains(obj)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Compute alpha 
obj.alpha = obj.Hhat * obj.ymvec; 

end