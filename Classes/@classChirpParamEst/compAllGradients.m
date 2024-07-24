function obj = compAllGradients(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Cache some variables for speed
Nc    = obj.Nc;
Pc    = obj.Pc;
ymvec = obj.ymvec;
Po    = obj.Po;
Hhat  = obj.Hhat;

k = 1;
for c = 1:Nc

    obj.c = c;

    % Compute the basis signal once for beta and phi
    obj  = obj.compBasisSignals(); % I think we can reuse H here. Optimize!
    temp = obj.x;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Gradient of J wrt beta ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Get the gradient of H wrt beta first
    obj.dH_beta(:,c) = [-(obj.n / obj.fs) .* real(temp) ; ...
                        -(obj.n / obj.fs) .* imag(temp)];

    obj.dJ_beta(1,c) = -2 * ymvec.' * Po * obj.dH_beta(:,c) * Hhat(c,:) * ymvec;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Gradient of J wrt phi ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for p = 1:Pc(c)

        obj.p = p;
        obj.k = k;

        % Get the gradient of H wrt phi
        obj.dH_phi(:,c,k) = [-(obj.n ./ obj.fs).^p .* imag(temp) ; ...
                              (obj.n ./ obj.fs).^p .* real(temp)];

        obj.dJ_phi(1,k) = -2 * ymvec.' * Po * obj.dH_phi(:,c,k) * Hhat(c,:) * ymvec;

        k = k + 1;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Gradient of J wrt gamma ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    obj.bAmpGamma = true;
    obj  = obj.compBasisSignals(); % Now we need a different basis signals for gamma
    temp = obj.x;
    obj.bAmpGamma = false;

    % Get the gradient of H wrt gamma 
    obj.dH_gamma(:,c) = [(obj.n / obj.fs) .* real(temp) ; ...
                         (obj.n / obj.fs) .* imag(temp)];

    obj.dJ_gamma(1,c) = -2 * ymvec.' * Po * obj.dH_gamma(:,c) * Hhat(c,:) * ymvec;
end

end