function obj = compGradients(obj)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Get gradients of H with respect to phi, beta and gamma
[dH_phi, dH_beta, dH_gamma] = genAllGradH(fs, beta, gamma, phi, c, p);

% Get gradient of H col wrt phi, beta and gamma (This is too expensive!)
dJ_phi   = 2 * pagemtimes(y.' * Po, pagemtimes(dH_phi,   Hhat * y));
dJ_beta  = 2 * pagemtimes(y.' * Po, pagemtimes(dH_beta,  Hhat * y));
dJ_gamma = 2 * pagemtimes(y.' * Po, pagemtimes(dH_gamma, Hhat * y));

end