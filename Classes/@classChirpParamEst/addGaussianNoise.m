function obj = addGaussianNoise(obj)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if ~isreal(obj.xm)
    obj.wm = randn(size(obj.xm)) + 1j * randn(size(obj.xm)); %  random noise signal
else
    obj.wm = randn(size(obj.xm)); %  random noise signal
end

xpow = sum(abs(obj.xm.^2)) ;
wpow = sum(abs(obj.wm.^2));

pn = xpow / (10^(obj.snr/10)); % required SNR
wscaled = (obj.wm ./ sqrt(wpow)) .* sqrt(pn); % normalizing the noisy signal to get unity power then multiplying it by the new power to achieve the required SNR

obj.ym = obj.xm + wscaled; % add noise to signal

% wnewpow = sum(abs(wscaled.^2));
% snrtrue = 10*log10(xpow / wnewpow);

end