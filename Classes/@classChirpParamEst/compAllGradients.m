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

k = 0;
for c = 1:Nc

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Gradient of J wrt phi ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for p = 0:Pc(c)-1
       
        k = k + 1;

        % Get the gradient of H wrt phi
        obj.dH_phi(:,c,k) = 2*pi*1j * (obj.n .* oneOverFs).^p .* H(:,c);
        obj.dJ_phi(1,k)   = -2 * real(ym' * Po * obj.dH_phi(:,c,k) * Hhat(c,:) * ym);     
    end

end

end