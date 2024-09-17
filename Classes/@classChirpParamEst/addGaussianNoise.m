function obj = addGaussianNoise(obj)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% if ~isreal(obj.ym)
%     obj.wm = randn(size(obj.ym)) + 1j * randn(size(obj.ym)); %  random noise signal
% else
%     obj.wm = randn(size(obj.ym)); %  random noise signal
% end
% 
% ypow = sum(abs(obj.ym.^2)) ;
% wpow = sum(abs(obj.wm.^2));
% 
% pn  = ypow / (10^(obj.snr/10)); % required SNR
% 
% obj.sigma = sqrt(pn) / sqrt(wpow);
% 
% wsc = obj.wm .* obj.sigma; % normalizing the noisy signal to get unity power then multiplying it by the new power to achieve the required SNR
% 
% obj.ym = obj.ym + wsc; % add noise to signal

% signal: Input complex exponential signal
% snr_db: Desired SNR in dB
% noisy_signal: Output signal with added complex Gaussian noise
% noise_variance: Variance of the added noise

% Step 1: Compute signal power
ypow = mean(abs(obj.ym).^2);

% Step 2: Convert SNR from dB to linear scale
snr_linear = 10^(obj.snr / 10);

% Step 3: Compute noise power based on SNR
wpow = ypow / snr_linear;

% Step 4: Compute noise variance (since it's complex, the variance is split between real and imaginary parts)
wvar = wpow / 2;

% Step 5: Generate complex Gaussian noise
wreal = sqrt(wvar) * randn(size(obj.ym));
wimag = sqrt(wvar) * randn(size(obj.ym));
obj.wm = wreal + 1i * wimag;

% Step 6: Add noise to the signal
obj.ym = obj.ym + obj.wm;

% Output noise variance
obj.sigma2 = wvar;

% Check
wnewpow = mean(abs(obj.wm).^2);
snrtrue = 10*log10(ypow / wnewpow);

end

