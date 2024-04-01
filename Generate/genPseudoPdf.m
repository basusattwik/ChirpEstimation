function [normPseudoPdf, pseudoPdf] = genPseudoPdf(x, nPhi, rho, M)
%GENPSEUDOPDF Summary of this function goes here
%   Detailed explanation goes here

P = 1; % one chirp at a time
alphaGrid = linspace(0, 1, M);
betaGrid  = linspace(0, 2, M);
propFunc  = zeros(M * ones(1,nPhi));
for i = 1:M
    for j = 1:M
        propFunc(i,j) = genProposalFunc(x, P, alphaGrid(i), betaGrid(j), rho);
        % ch = genExpPolyChirp3(1, N, 1, -[paramGrid(1,i), paramGrid(1,j)]);
        % propFunc(i,j) = (1/N) * abs(x' * ch).^2;
    end
end
pseudoPdf = exp(real(propFunc));

% Normalizing factor
normFactor    = trapz(alphaGrid, trapz(betaGrid, pseudoPdf, 2));
normPseudoPdf = pseudoPdf ./ (normFactor);

end