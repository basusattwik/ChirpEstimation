function obj = genChirpSignal(obj)
%GENCHIRPSIGNAL Summary of this function goes here
%   Detailed explanation goes here

Nc = obj.Nc;
N  = obj.N;
fs = obj.fs;

if strcmpi(obj.mode, 'gen')

    % Calculate the polynomial chirp
    for c = 1:Nc  % -- loop over number of chirps
    
        % max polynomial degree
        P = size(obj.phi{1,c},1);
    
        for n = 1:N % -- loop over number of samples
    
            % Get the amplitude envelope
            obj.Am(n,c) = obj.alpha(1,c) * exp(-obj.beta(1,c) * n/fs) * (1 - exp(-obj.gamma(1,c) * n/fs));
    
            % Get the exponential polynomial phase sinusoid
            npvec = ((n-1) / fs).^(0:P-1).'; % vectors of powers of n/fs
            obj.em(n,c) = exp(2*pi*1j .* (obj.phi{:,c}.' * npvec));
    
        end % -- end loop over number of chirps
    
        % Multiply the two
        obj.xm(:,c) = obj.Am(:,c) .* obj.em(:,c);
    
    end % -- end loop over number of samples
    
    % Combined to form multicomponent signal
    obj.ym = sum(obj.xm, 2);

elseif strcmpi(obj.mode, 'comp')

    % max polynomial degree
    P = size(obj.phi, 1);

    for n = 1:N % -- loop over number of samples

        % Get the amplitude envelop
        obj.A(n,1) = exp(-obj.beta * n/fs) * (1 - exp(-obj.gamma * n/fs));

        % Get the exponential polynomial phase sinusoid
        npvec = ((n-1) / fs).^(0:P-1).'; % vectors of powers of n/fs
        obj.e(n,1) = exp(2*pi*1j .* (obj.phi.' * npvec));

    end

    % Multiply the two
    obj.x = obj.A .* obj.e;

end


