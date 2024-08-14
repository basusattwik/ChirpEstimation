function obj = runCpeCore(obj, params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N  = obj.N;
Nc = obj.Nc;
Ac = obj.Ac;
Pc = obj.Pc;
fs = obj.fs;

% Avoiding divides
oneOverFs = 1/fs;

startInd = 1;
for c = 1:Nc
    obj.phiEstCell{1,c} = params(startInd:startInd+Pc(c)-2); %  minus 2 because we are not including phase offset phi0
    startInd = startInd + Pc(c) - 1;
end

startInd = 1;
for c = 1:Nc

    phi  = obj.phiEstCell{1,c};
    pvec = (1:Pc(c)-1).';
    avec = (1:Ac(c)-1);

    endInd = startInd + Ac(c) - 1;
    for n = 1:N % -- loop over number of samples

        % Get the exponential polynomial phase sinusoid
        nind  = (n-1) * oneOverFs;
        npvec = nind.^pvec; % vectors of powers of n/fs
        navec = nind.^avec;
        obj.H(n, startInd:endInd) = [1, navec] .* exp(2*pi*1j .* (phi.' * npvec));
    end
    startInd = startInd + Ac(c);
end

% Compute Hhat
% Get the projection matrix and the orthogonal projection matrix
obj.Hhat = (obj.H' * obj.H)^(-1) * obj.H';
obj.P    = obj.H * obj.Hhat;
obj.Po   = eye(obj.N, obj.N) - obj.P;

k = 0;
startInd = 1;
for c = 1:Nc

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % --- Gradient of J wrt phi ---
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    endInd = startInd + Ac(c) - 1;
    for p = 1:Pc(c)-Nc

        k = k+1;

        % Get the gradient of H wrt phi
        obj.dH_phi(:,startInd:endInd,k) = 2*pi*1j .* (obj.n .* oneOverFs) .* obj.H(:,startInd:endInd);
        obj.dJ_phi(1,k) = -2 * real(obj.ym' * (obj.Po * (obj.dH_phi(:,startInd:endInd,k) * (obj.Hhat(startInd:endInd,:) * obj.ym))));     
    end
    startInd = startInd + Ac(c);
end

% Objective function value: want to minimize this
obj.J = real(obj.ym' * (obj.Po * obj.ym)); % Force it to be real to prevent tiny imaginary values ~ e-16


end
