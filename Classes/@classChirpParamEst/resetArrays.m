function obj = resetArrays(obj)
%RESETARRAYS Summary of this function goes here
%   Detailed explanation goes here

obj.a   = zeros(obj.N, 1); % not used
obj.e   = zeros(obj.N, obj.Ka);
obj.x   = zeros(obj.N, 1);
obj.xg  = zeros(obj.N, 1);
obj.w   = zeros(obj.N, 1);

obj.am  = zeros(obj.N, obj.Nc); % not used
obj.em  = zeros(obj.N, obj.Nc);
obj.xm  = zeros(obj.N, 1);
obj.ym  = zeros(obj.N, 1);
obj.wm  = zeros(obj.N, 1);

obj.alphaEst   = zeros(1, obj.Ka);
obj.phiEst     = zeros(1, obj.Kp);
obj.rhoEst     = zeros(1, obj.Ka);
obj.phiEstCell = cell(1,  obj.Nc);
obj.rhoEstCell = cell(1,  obj.Nc);

obj.J        = 0;
obj.H        = zeros(obj.N,  obj.Ka);
obj.Hhat     = zeros(obj.Ka, obj.N);
obj.P        = zeros(obj.N,  obj.N);
obj.Po       = zeros(obj.N,  obj.N);
obj.dH_phi   = zeros(obj.N,  obj.Ka, obj.Kp - obj.Nc);
obj.dJ_phi   = zeros(1, obj.Kp - obj.Nc);
obj.Id       = zeros(obj.N, obj.N);


end