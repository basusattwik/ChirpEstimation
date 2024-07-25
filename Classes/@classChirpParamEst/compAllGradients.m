function obj = compAllGradients(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Cache some variables for speed
fs = obj.fs;
Nc = obj.Nc;
Pc = obj.Pc;
ym = obj.ym;
Po = obj.Po;
H  = obj.H;
Hhat = obj.Hhat;

% Avoiding divides
oneOverFs = 1/fs;

k = 1;
for c = 1:Nc

    obj.c = c;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Gradient of J wrt beta ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Get the gradient of H wrt beta first
    obj.dH_beta(:,c) = -(obj.n .* oneOverFs) .* H(:,c);
    obj.dJ_beta(1,c) = real(-2 * ym' * Po * obj.dH_beta(:,c) * Hhat(c,:) * ym);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Gradient of J wrt phi ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for p = 1:Pc(c)

        obj.p = p;
        obj.k = k;

        % Get the gradient of H wrt phi
        obj.dH_phi(:,c,k) = 2*pi*1j * (obj.n .* oneOverFs).^p .* H(:,c);
        obj.dJ_phi(1,k)   = real(-2 * ym' * Po * obj.dH_phi(:,c,k) * Hhat(c,:) * ym);

        k = k + 1;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Gradient of J wrt gamma ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    obj.bAmpGamma = true;
    obj = obj.compBasisSignals(); % Now we need a different basis signals for gamma
    obj.bAmpGamma = false;

    % Get the gradient of H wrt gamma 
    obj.dH_gamma(:,c) = (obj.n .* oneOverFs) .* obj.xg;
    obj.dJ_gamma(1,c) = real(-2 * ym' * Po * obj.dH_gamma(:,c) * Hhat(c,:) * ym);
end

end