function J = evalObjectiveFunc(obj, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Cache for speed
fs = obj.fs;
N  = obj.N;
Nc = obj.Nc;
Pc = obj.Pc;
Ac = obj.Ac;
ym = obj.ym;
Ka = obj.Ka;

% Preallocate for speed
phiCell = cell(1, Nc);
H = zeros(N, Ka);

% Avoiding divides
oneOverFs = 1/fs;

% Extract current params into data structs
startInd = 1;
for c = 1:Nc
    phiCell{1,c} = params(startInd:startInd+Pc(c)-2);
    startInd = startInd + Pc(c) - 1;
end

% % Compute basis matrix
% for c = 1:Nc
%     P = size(phi{1,c},1);
%     pvec = (1:P-1).';
%     for n = 1:N % -- loop over number of samples
% 
%         % Get the exponential polynomial phase sinusoid
%         npvec  = ((n-1) * oneOverFs).^pvec; % vectors of powers of n/fs
%         e(n,1) = exp(2*pi*1j .* (phi{1,c}.' * npvec));
%     end
%     x = e;
%     H(:,c) = x;
% end

startInd = 1;
for c = 1:Nc

    phi  = phiCell{1,c};
    pvec = (1:Pc(c)-1).';
    avec = (1:Ac(c)-1);

    endInd = startInd + Ac(c) - 1;
    for n = 1:N % -- loop over number of samples

        % Get the exponential polynomial phase sinusoid
        nind  = (n-1) * oneOverFs;
        npvec = nind.^pvec; % vectors of powers of n/fs
        navec = nind.^avec;
        H(n, startInd:endInd) = [1, navec] .* exp(2*pi*1j .* (phi.' * npvec));
    end
    startInd = startInd + Ac(c);
end

% Get the projection matrix and the orthogonal projection matrix
P  = H * (H' * H)^(-1) * H'; 
Po = eye(size(P)) - P;

% Objective function value
J = real(ym' * (Po * ym));

end