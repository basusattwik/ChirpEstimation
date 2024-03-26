function [normPseudoPdf, pseudoPdf] = genPseudoPdf(x, phi, rho, M, Fs)
%GENPSEUDOPDF Summary of this function goes here
%   Detailed explanation goes here

lenSignal = length(x);
numParams = size(phi,2);

paramGrid = repmat(linspace(0, 1, M), [numParams, 1]);
pseudoPdf = zeros(M * ones(1,numParams));
for i = 1:M
    for j = 1:M
        c = genExpPolyChirp2(Fs, lenSignal / Fs, 1, -[paramGrid(1,i), paramGrid(1,j)]);
        pseudoPdf(i,j) = exp(rho * (1/lenSignal) * abs(x.' * c)^2);
    end
end

normFactor    = (1/M^(numParams)) .* sum(pseudoPdf, 'all');
normPseudoPdf = pseudoPdf ./ normFactor;

% Maybe think of using numerical integration cumtrapz twice

end