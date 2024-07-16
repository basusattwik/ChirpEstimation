function y = genChirpSignal(fs, Nc, Td, alpha, beta, gamma, phi)
%GENCHIRPSIGNAL Summary of this function goes here
%   Detailed explanation goes here

if isempty(Nc)
    Nc = 1;
end
if isempty(alpha)
    alpha = 1;
end
if size(phi,2) ~= Nc || size(alpha,2) ~= Nc || size(beta,2) ~= Nc || size(gamma,2) ~= Nc
    error('Number of chirps does not match number of parameters provided');
end

N = Td * fs;     % Number of samples

% Preallocate for speed
A = zeros(N, Nc); % amplitude envelop for each chirp
e = zeros(N, Nc); % individual chirps without amplitude envelop
x = zeros(N, Nc); % individual chirps with amplitude envelop

% Calculate the polynomial chirp
for c = 1:Nc  % -- loop over number of chirps

    % max polynomial degree
    P = size(phi{1,c},1);

    for n = 1:N % -- loop over number of samples

        % Get the amplitude envelop
        A(n,c) = alpha(1,c) * exp(-beta(1,c) * n/fs) * (1 - exp(-gamma(1,c) * n/fs));

        % Get the exponential polynomial phase sinusoid
        npvec  = ((n-1) / fs).^(0:P-1).'; % vectors of powers of n/fs
        e(n,c) = exp(2*pi*1j .* (phi{:,c}.' * npvec));

    end % -- end loop over number of chirps

    % Multiply the two
    x(:,c) = A(:,c) .* e(:,c);

end % -- end loop over number of samples

% Combined multicomponent signal
y = sum(x,2);

end


