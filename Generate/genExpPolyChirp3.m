function [xtot, n] = genExpPolyChirp3(Fs, Td, A, phi)
%GENEXPPOLYCHIRP Creates a exponential polynomial phase sinusoid (no phase
%offset)

    D = size(phi,2); % polynomial degree
    N = Td * Fs;     % Number of samples
    P = size(phi,1); % Number of parameter sets.. i.e., number of chirps

    % Calculate the polynomial chirp
    n = (0:N-1).';
    x = zeros(N,1);    % individual chirps
    xtot = zeros(N,1); % combined chirps
    for p = 1:P
        for nind = 1:N        
            x(nind,1) = exp(2*pi*1j .* (phi(p,:) .* [1, 0.5]) * (n(nind) / Fs).^(1:D).');
        end
        xtot = xtot + A(p) * x;
        x(:) = 0;
    end 
end