function y = addWhiteGaussianNoise(x, snr)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if ~isreal(x)
    w = randn(size(x)) + 1j * randn(size(x)); %  random noise signal
else
    w = randn(size(x)); %  random noise signal
end

xpow = sum(abs(x.^2)) ;
wpow = sum(abs(w.^2));

pn = xpow / (10^(snr/10)); % required SNR
wscaled = (w ./ sqrt(wpow)) .* sqrt(pn); % normalizing the noisy signal to get unity power then multiplying it by the new power to achieve the required SNR

y = x + wscaled; % add noise to signal

% wnewpow = sum(abs(wscaled.^2));
% snrtrue = 10*log10(xpow / wnewpow);

end