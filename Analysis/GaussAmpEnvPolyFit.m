
clc
close all
clearvars

% Parameters
sample_rate = 3000;   % Sampling rate in Hz
duration = 1;         % Duration of the signal in seconds
variance = 0.035;      % Variance of the Gaussian function

% Generate Gaussian envelope
[gaussian_envelope, t] = generate_gaussian_envelope(sample_rate, duration, variance);

ord   = 6;
coeff = polyfit(t, gaussian_envelope, ord);
poly  = polyval(coeff, t);

figure('WindowStyle','docked')

plot(t, gaussian_envelope); hold on;
plot(t, (poly));
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

function [gaussian_envelope, t] = generate_gaussian_envelope(sample_rate, duration, variance)
    % Generate Gaussian envelope
    % sample_rate: sampling rate in Hz
    % duration: duration of the signal in seconds
    % variance: variance of the Gaussian function
   
    
    % Time vector starting from 0
    t = 0:1/sample_rate:duration-1/sample_rate; %linspace(0, duration, num_samples);
    
    % Center of the Gaussian function
    center = duration / 2;
    
    % Gaussian envelope
    gaussian_envelope = exp(-(t - center).^2 / (2 * variance));
    
    % Plotting the Gaussian envelope
    % figure;
    % plot(t, gaussian_envelope);
    % title('Gaussian Envelope');
    % xlabel('Time (s)');
    % ylabel('Amplitude');
    % grid on;
end
