function obj = runCpeCore(obj, params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Nc = obj.Nc;
Ac = obj.Ac;
Pc = obj.Pc;
fs = obj.fs;

% Avoiding divides
nvecOverFs = obj.n / fs;
twoPij  = 2*pi*1j;

startInd = 1;
for c = 1:Nc
    obj.phiEstCell{1,c} = params(startInd:startInd+Pc(c)-2); %  minus 2 because we are not including phase offset phi0
    startInd = startInd + Pc(c)-1;
end

%%%%%%%%%%%%%%%%%%%%%%
%                    %
% --- Basis Matrix ---
%                    %
%%%%%%%%%%%%%%%%%%%%%%

% Create H matrix (basis vectors)
startInd = 1;
for c = 1:Nc

    phi  = obj.phiEstCell{1,c};
    pvec = (1:Pc(c)-1);
    avec = (0:Ac(c)-1);

    endInd = startInd+Ac(c)-1;
    obj.H(:, startInd:endInd) = (nvecOverFs.^avec) .* exp(twoPij .* ((nvecOverFs.^pvec) * phi));
    startInd = startInd + Ac(c);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
% --- Projection & Annhialator Matrices ---
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the projection matrix and the orthogonal projection matrix
obj.Hhat = (obj.H' * obj.H)^(-1) * obj.H';
obj.P    = obj.H * obj.Hhat;
obj.Po   = obj.Id - obj.P; % Can optimize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             %
% --- Gradient of J wrt phi ---
%                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = 0;
startInd = 1;
for c = 1:Nc

    endInd = startInd + Ac(c) - 1;
    for p = 1:Pc(c)-1
        k = k+1;

        % Gradient of H wrt phi
        obj.dH_phi(:,startInd:endInd,k) = twoPij .* nvecOverFs.^p .* obj.H(:,startInd:endInd);

        % Gradient of J wrt phi
        obj.dJ_phi(1,k) = -2 * real(obj.ym' * (obj.Po * (obj.dH_phi(:,startInd:endInd,k) * (obj.Hhat(startInd:endInd,:) * obj.ym))));

    end
    startInd = startInd + Ac(c);
end

% Objective function value: want to minimize this
obj.J = real(obj.ym' * (obj.Po * obj.ym)); % Force it to be real to prevent tiny imaginary values ~ e-16

end
