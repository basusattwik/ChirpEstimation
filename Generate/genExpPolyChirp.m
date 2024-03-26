function [xpps, iphs, t] = genExpPolyChirp(Fs, Td, A, phi)
%GENEXPPOLYCHIRP Creates a exponential polynomial phase sinusoid

    t = (0:1/Fs:Td-1/Fs).';
    N = numel(phi);
    L = Td*Fs;

    % Calculate the polynomial phase sinusoid
    xpps = zeros(L,1);
    iphs = zeros(L,1);
    for tind = 1:L        
        iphs(tind,1) = phi.' * t(tind).^(0:N-1).';
        xpps(tind,1) = exp(2*pi*1j * iphs(tind,1));
    end
    xpps = A * xpps;
end