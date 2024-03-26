function [xpps, iphs, t] = genPolyChirp(Fs, Td, A, phi)
%GENPOLYCHIRP Creates a polynomial phase sinusoid

    t = (0:1/Fs:Td-1/Fs).';
    N = numel(phi);
    L = Td*Fs;

    % Calculate the polynomial phase sinusoid
    xpps = zeros(L,1);
    iphs = zeros(L,1);
    for tind = 1:L        
        iphs(tind,1) = phi.' * t(tind).^(1:N).';
        xpps(tind,1) = sin(2*pi * iphs(tind,1));
    end
    xpps = A * xpps;
end