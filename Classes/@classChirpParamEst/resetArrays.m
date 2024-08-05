function obj = resetArrays(obj)
%RESETARRAYS Summary of this function goes here
%   Detailed explanation goes here

obj.e   = zeros(obj.N, 1);
obj.u   = zeros(obj.N, 1);
obj.x   = zeros(obj.N, 1);
obj.xg  = zeros(obj.N, 1);
obj.w   = zeros(obj.N, 1);

obj.em  = zeros(obj.N, obj.Nc);
obj.um  = zeros(obj.N, obj.Nc);
obj.xm  = zeros(obj.N, 1);
obj.ym  = zeros(obj.N, 1);
obj.wm  = zeros(obj.N, 1);

obj.alphaEst = zeros(1, obj.Nc);
obj.phiEst   = zeros(1, obj.K);
obj.phiEstCell = cell(1, obj.Nc);

obj.J        = 0;
obj.H        = zeros(obj.N, obj.Nc);
obj.Hhat     = zeros(obj.Nc, obj.N);
obj.P        = zeros(obj.N, obj.N);
obj.Po       = zeros(obj.N, obj.N);
obj.dH_phi   = zeros(obj.N, obj.Nc, obj.K);
obj.dJ_phi   = zeros(1, obj.K);


end