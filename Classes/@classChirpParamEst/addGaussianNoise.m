function obj = addGaussianNoise(obj)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if ~isreal(obj.ym)
    obj.wm = randn(size(obj.ym)) + 1j * randn(size(obj.ym)); %  random noise signal
else
    obj.wm = randn(size(obj.ym)); %  random noise signal
end

ypow = sum(abs(obj.ym.^2)) ;
wpow = sum(abs(obj.wm.^2));

pn  = ypow / (10^(obj.snr/10)); % required SNR
wsc = (obj.wm ./ sqrt(wpow)) .* sqrt(pn); % normalizing the noisy signal to get unity power then multiplying it by the new power to achieve the required SNR

obj.ym = obj.ym + wsc; % add noise to signal

% wnewpow = sum(abs(wscaled.^2));
% snrtrue = 10*log10(ypow / wnewpow);

end