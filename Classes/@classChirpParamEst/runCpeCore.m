function obj = runCpeCore(obj, params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% % Reformat arrays to cells
% obj = obj.convertParamArray2Cell(params);
% 
% % Update the basis matrix
% obj = obj.compBasisMatrix();
% 
% % Update the projection matrix
% obj = obj.compProjMatrix();
% 
% % Compute new gradients based on the new parameters
% obj = obj.compAllGradients();
% 
% % Compute the new cost function value
% obj = obj.compObjectiveFunc();

% -----------
N  = obj.N;
Nc = obj.Nc;
Pc = obj.Pc;
fs = obj.fs;

% Avoiding divides
oneOverFs = 1/fs;

startInd = 1;
for c = 1:Nc
    obj.phiEstCell{1,c} = params(startInd:startInd+Pc(c)-1);
    startInd = startInd + Pc(c);
end

% Fill columns with the basis signal 
for c = 1:Nc

    P = size(obj.phiEstCell{1,c},1);
    pvec = (0:P-1).';
    for n = 1:N % -- loop over number of samples

        % Get the exponential polynomial phase sinusoid
        npvec = ((n-1) * oneOverFs).^pvec; % vectors of powers of n/fs
        obj.e(n,1) = exp(2*pi*1j .* (obj.phiEstCell{1,c}.' * npvec));

    end

    obj.x = obj.e;
    obj.H(:,c) = obj.x;
end

% Compute Hhat
% Get the projection matrix and the orthogonal projection matrix
obj.P  = obj.H * (obj.H' * obj.H)^(-1) * obj.H';
obj.Po = eye(obj.N, obj.N) - obj.P;

k = 0;
for c = 1:Nc

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % --- Gradient of J wrt phi ---
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for p = 0:Pc(c)-1

        k = k+1;

        % Get the gradient of H wrt phi
        obj.dH_phi(:,c,k) = 2*pi*1j * (obj.n .* oneOverFs).^p .* obj.H(:,c);
        obj.dJ_phi(1,k)   = -2 * real(obj.ym' * obj.Po * obj.dH_phi(:,:,k) * ((obj.H' * obj.H)^(-1) * obj.H') * obj.ym);     
    end

end

% Objective function value: want to minimize this
obj.J = real(obj.ym' * obj.Po * obj.ym); % Force it to be real to prevent tiny imaginary values ~ e-16


end
