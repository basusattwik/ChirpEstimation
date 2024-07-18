function obj = compGradsForJ(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Nc = obj.Nc;
Pc = obj.Pc;

y    = obj.y;
Po   = obj.Po;
Hhat = obj.Hhat;

k = 1;
for c = 1:Nc

    obj.c = c;

    % Compute the basis signal once for beta and phi
    temp = obj.compBasisSignals();

    % --- Gradient of H wrt beta ---
    % obj.compGradsForH_beta();
    obj.dH_beta(:,obj.c) = [-(obj.n / obj.fs) .* real(temp) ; ...
                            -(obj.n / obj.fs) .* imag(temp)];
    obj.dJ_beta(1,c)  = 2 * y.' * Po * obj.dH_beta(:,c)  * Hhat(c,:) * y;

    % --- Gradient of H wrt phi ---
    for p = 1:Pc

        obj.p = p;
        obj.k = k;
        % obj.compGradsForH_phi();
        % temp = obj.compBasisSignals();
        obj.dH_phi(:, obj.c, obj.k) = [-(obj.n / obj.fs)^obj.p .* imag(temp) ; ...
                                       +(obj.n / obj.fs)^obj.p .* real(temp)];
        obj.dJ_phi(1,k) = 2 * y.' * Po * squeeze(obj.dH_phi(:,c,k)) * Hhat(c,:) * y;
        k = k+1;
    end

    % --- Gradient of H wrt gamma ---
    % obj.compGradsForH_gamma();
    obj.bAmpGamma = true;
    temp = obj.compBasisSignals(); % Now we need a different basis signals for gamma
    obj.dH_gamma(:,obj.c) = [-(obj.n / obj.fs) .* real(temp) ; ...
                             -(obj.n / obj.fs) .* imag(temp)];
    obj.bAmpGamma = false;
    obj.dJ_gamma(1,c) = 2 * y.' * Po * obj.dH_gamma(:,c) * Hhat(c,:) * y;

end

end