function [energy_decay, log_f] = calcEDR(s,f,num_freq_to_plot)
% calculate energy decay relief based on input spectrogram
% requires audio toolbox for hz2mel, mel2hz

arguments
    s = []; % magnitude spectrogram, nfreq x ntime
    f = []; % frequency values, nfreq
    num_freq_to_plot = 50;
end

% s = clip(s,1e-15, 1);

% Backwards integration
M = size(s,2);
energy_decay = zeros(num_freq_to_plot,size(s,2));

% Remap frequencies in log space
log_f=logspace(log10(min(f)), log10(max(f)), num_freq_to_plot);
log_y=interp1(f,s,log_f);

for m = 1:M
    energy_decay(:,m) = 10*log10(sum(log_y(:,m:end),2));
end

energy_decay = energy_decay - max(energy_decay,[],"all");