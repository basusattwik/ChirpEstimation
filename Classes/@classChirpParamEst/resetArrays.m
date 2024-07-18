function obj = resetArrays(obj)
%RESETARRAYS Summary of this function goes here
%   Detailed explanation goes here

obj.Pc  = zeros(obj.Nc, 1);
obj.A   = zeros(2*obj.N, obj.Nc);
obj.e   = zeros(2*obj.N, obj.Nc);
obj.u   = zeros(2*obj.N, obj.Nc);
obj.x   = zeros(2*obj.N, 1);
obj.w   = zeros(2*obj.N, 1);

obj.Am  = zeros(2*obj.N, obj.Nc);
obj.em  = zeros(2*obj.N, obj.Nc);
obj.um  = zeros(2*obj.N, obj.Nc);
obj.xm  = zeros(2*obj.N, 1);
obj.ym  = zeros(2*obj.N, 1);
obj.wm  = zeros(2*obj.N, 1);

obj.alphaEst = zeros(1, obj.Nc);
obj.betaEst  = zeros(1, obj.Nc); 
obj.gammaEst = zeros(1, obj.Nc);
obj.phiEst   = zeros(1, obj.K);
obj.phiEstCell = cell(1, obj.Nc);

obj.J        = zeros(2*obj.N, 1);
obj.H        = zeros(2*obj.N, obj.Nc);
obj.Hhat     = zeros(obj.Nc,  2*obj.N);
obj.P        = zeros(2*obj.N, 2*obj.N);
obj.Po       = zeros(2*obj.N, 2*obj.N);
obj.dH_phi   = zeros(obj.N, obj.Nc, obj.K);
obj.dH_beta  = zeros(obj.N, obj.Nc);
obj.dH_gamma = zeros(obj.N, obj.Nc);
obj.dJ_phi   = zeros(1, obj.K);
obj.dJ_beta  = zeros(1, obj.Nc);
obj.dJ_gamma = zeros(1, obj.Nc);


end