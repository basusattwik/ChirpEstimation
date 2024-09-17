function obj = evalParamErrors(obj)
%EVALERRORS
% Details

Nc = obj.Nc;
fs = obj.fs;
Pc = obj.Pc;
% Ac = cpe.Ac;

% Get the parameters
phiEstVec = []; %zeros(obj.numParams, 1);
phiActVec = []; %zeros(obj.numParams, 1);
rhoEstVec = [];
rhoActVec = [];

% Convert cells to array
for c = 1:Nc
    phiActVec = [phiActVec ; obj.phi{1,c}(2:end,1)];
    rhoActVec = [rhoActVec ; obj.rho{1,c}];
end

q = zeros(Nc,1);
for c1 = 1:Nc
    for c2 = 1:Nc
        q(c2,1) = (obj.phi{1,c1}(2,1) - obj.phiEstCell{1,c2}(1,1)).^2;
    end
    [~, m] = min(q);
    phiEstVec = [phiEstVec; obj.phiEstCell{1,m}];
    rhoEstVec = [rhoEstVec; obj.rhoEstCell{1,m}];
    q(:) = 0;
end

% q = zeros(Nc,1);
% for c1 = 1:Nc
%     for c2 = 1:Nc
%         q(c1,1) = sum((cpe.phi{1,c1}(2:end,1) - cpe.phiEstCell{1,c2}).^2);
%     end
%     [~, m] = min(q);
%     phiEstVec = [phiEstVec; cpe.phiEstCell{1,m}];
%     q(:) = 0;
% end

normPhiError = (phiActVec - phiEstVec) .* repmat((1/fs).^(1:Pc(1,1)-1).', Nc, 1);
obj.sqrPhiError = normPhiError.^2;

normRhoError = (rhoActVec - rhoEstVec);
obj.sqrRhoError = normRhoError.^2;

end
