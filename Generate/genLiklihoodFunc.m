function Lc = genLiklihoodFunc(x, P, alpha, beta, lambda)
%GENLIKLIHOODFUNC Summary of this function goes here
%   Detailed explanation goes here

N = length(x);

H = zeros(N, P);
A = ones(P,1);
for p = 1:P
    H(:,p) = genExpPolyChirp3(1, N, A, [alpha(p), beta(p)]);
end

% Likelihood Function
Lc = lambda * x' * (H * pinv(H' * H) * H') * x;

end