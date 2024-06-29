function [normPseudoPdf, pseudoPdf] = genPseudoPdf(x, rho, M)
%GENPSEUDOPDF Summary of this function goes here
%   Detailed explanation goes here

P = 1; % one chirp at a time
alphaGrid = linspace(0, 1/2, M);
betaGrid  = linspace(0, 2/2, M);
propFunc  = zeros(M, M); %zeros(M * ones(1,nPhi));
for i = 1:M
    for j = 1:M
        propFunc(i,j) = genProposalFunc(x, P, alphaGrid(i), betaGrid(j), rho);
    end
end
pseudoPdf = exp(real(propFunc));

% Normalizing factor
normFactor    = trapz(alphaGrid, trapz(betaGrid, pseudoPdf, 2));
normPseudoPdf = pseudoPdf ./ (normFactor);

end