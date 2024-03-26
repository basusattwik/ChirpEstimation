function I = genProposalDist(x, Fs, alpha, beta, rho, P, M)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

N = length(x);

gdash = 1;
for p = 1:P
    for i = 1:M
        for j = 1:M
            ch = genPolynomialChirp(Fs, N ,1, [0; alpha(i); beta(i)]);
            I(p,i,j) = exp(rho * (1/N) * sum(x .* ch));
        end
    end
    gdash = gdash * I;
end

end