function [xtot, t] = genExpPolyChirp2(Fs, Td, A, phi)
%GENEXPPOLYCHIRP Creates a exponential polynomial phase sinusoid (no phase
%offset)

    t = (0:1/Fs:Td-1/Fs).';
    D = size(phi,2); % polynomial degree
    L = Td * Fs;
    P = size(phi,1);

    % Calculate the polynomial chirp
    x    = zeros(L,1); % individual chirps
    xtot = zeros(L,1); % combined chirps
    for p = 1:P
        for tind = 1:L        
            x(tind,1) = exp(2*pi*1j .* phi(p,:) * t(tind).^(1:D).');
        end
        xtot = xtot + A(p) * x;
        x(:) = 0;
    end
end