function obj = evalErrors(obj)
%EVALERRORS
% Details

bp  = obj.bestParticleInd;
cpe = obj.cpe{1,bp};

Nc = cpe.Nc;
fs = cpe.fs;
% Pc = cpe.Pc;
% Ac = cpe.Ac;

% Get the parameters
phiEstVec = []; %zeros(obj.numParams, 1);
phiActVec = []; %zeros(obj.numParams, 1);

% Convert cells to array
for c = 1:Nc
    phiActVec = [phiActVec ; cpe.phi{1,c}(2:end,1)];
end

q = zeros(Nc,1);
for c1 = 1:Nc
    for c2 = 1:Nc
        q(c1,1) = sum((cpe.phi{1,c1}(2:end,1) - cpe.phiEstCell{1,c2}).^2);
    end
    [~, m] = min(q);
    phiEstVec = [phiEstVec; cpe.phiEstCell{1,m}];
    q(:) = 0;
end

normError    = (phiActVec - phiEstVec);% .* repmat((1/fs).^(1:obj.numParams/2).', obj.numParams/2, 1);
obj.sqrError = normError.^2;

end
