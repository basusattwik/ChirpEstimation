function obj = genChirpSignal(obj)
%GENCHIRPSIGNAL Summary of this function goes here
%   Detailed explanation goes here

Nc = obj.Nc;
N  = obj.N;
fs = obj.fs;
Pc = obj.Pc;
Ac = obj.Ac;

% Avoiding divides
oneOverFs = 1 / fs;

% Calculate the polynomial chirp (NOT vectorized... but does not run in loops)
for c = 1:Nc  % -- loop over number of chirps

    rho  = obj.rho{1,c};
    phi0 = obj.phi{1,c}(1,1);
    phi  = obj.phi{1,c}(2:end,1);

    for nind = 1:N % -- loop over number of samples

        % Get the exponential polynomial phase sinusoid with polynomial
        % amplitude
        npvec = ((nind-1) * oneOverFs).^(1:Pc(c)-1); % vectors of powers of n/fs
        nfvec = ((nind-1) * oneOverFs).^(0:Pc(c)-2);
        navec = ((nind-1) * oneOverFs).^(0:Ac(c)-1);

        obj.am(nind, c) = navec  * rho;
        obj.fim(nind,c) = nfvec  * ((1:Pc(c)-1).' .* phi);
        obj.em(nind, c) = exp(1j * (phi0 + 2*pi .* (npvec * phi))); % exp(j * phi_0 + 2pij * phi(t))

    end % -- end loop over number of chirps
end % -- end loop over number of samples

obj.xm = obj.am .* obj.em;

% Combined to form multicomponent signal
obj.ym = sum(obj.xm, 2);

end


