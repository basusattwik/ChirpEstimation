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
obj.fim = zeros(obj.N, obj.Nc);
obj.xm  = zeros(obj.N, 1);
obj.ym  = zeros(obj.N, 1);
obj.wm  = zeros(obj.N, 1);

obj.amRecon  = zeros(obj.N, obj.Nc); % not used
obj.emRecon  = zeros(obj.N, obj.Nc);
obj.fimRecon = zeros(obj.N, obj.Nc);
obj.xmRecon  = zeros(obj.N, 1);
obj.ymRecon  = zeros(obj.N, 1);

obj.bvec       = zeros(obj.Ka, 1);
% obj.phiEst     = zeros(1, obj.Kp);
obj.rhoEst     = zeros(1, obj.Ka);
obj.phi0Est    = zeros(1, obj.Nc);
obj.phiEstCell = cell(1,  obj.Nc);
obj.rhoEstCell = cell(1,  obj.Nc);

obj.J       = 0;
obj.H       = zeros(obj.N,  obj.Ka);
obj.Hhat    = zeros(obj.Ka, obj.N);
obj.P       = zeros(obj.N,  obj.N);
obj.Po      = zeros(obj.N,  obj.N);
obj.dH_phi  = zeros(obj.N,  obj.Ka, obj.Kp - obj.Nc);
obj.dJ_phi  = zeros(1, obj.Kp - obj.Nc);
obj.d2J_phi = zeros(obj.Kp - obj.Nc, obj.Kp - obj.Nc);
obj.Idn     = eye(obj.N, obj.N);
obj.Idk     = eye(obj.Ka, obj.Ka);

obj.fimat = zeros(obj.Kp - obj.Nc, obj.Kp - obj.Nc);
obj.crbnd = zeros(1, obj.Kp - obj.Nc);


end