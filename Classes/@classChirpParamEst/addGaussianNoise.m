function obj = addGaussianNoise(obj)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if ~isreal(obj.x)
    obj.w = randn(size(obj.x)) + 1j * randn(size(obj.x)); %  random noise signal
else
    obj.w = randn(size(obj.x)); %  random noise signal
end

xpow = sum(abs(obj.x.^2)) ;
wpow = sum(abs(obj.w.^2));

pn = xpow / (10^(obj.snr/10)); % required SNR
wscaled = (obj.w ./ sqrt(wpow)) .* sqrt(pn); % normalizing the noisy signal to get unity power then multiplying it by the new power to achieve the required SNR

obj.y = obj.x + wscaled; % add noise to signal

% wnewpow = sum(abs(wscaled.^2));
% snrtrue = 10*log10(xpow / wnewpow);

end