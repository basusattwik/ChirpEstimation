function g = genProposalFunc(x, P, alpha, beta, rho1)
%GENPROPOSALFUNC Summary of this function goes here
%   Detailed explanation goes here

N = length(x);
H = zeros(N,P);
A = ones(P,1);
Fs = 1; % Ignore Fs here
for p = 1:P
    H(:,p) = genExpPolyChirp3(Fs, N, A, [alpha(p), beta(p)]);
end

% Likelihood Function
g = real(rho1 * (1/N) * x' * (H * H') * x);

end