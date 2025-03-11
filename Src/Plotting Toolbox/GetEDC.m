% 
% This function plots an Energy Decay Curve (aka Schroeder curve) of an IR
% using backwards integration.
%
% Optional Parameters:
%
% line_style - defaults to solid ("-")
% x_max - Maximum time value to plot in seconds
% y_min - Minimum energy value to plot in dB (defaults to -70 dB)
% integration_limit_sec - integration limit in seconds (defaults to the length of the IR)

function [integrated_values_dB, time_values_seconds] = GetEDC(ir, sample_rate, octave_band_centre, integration_limit_sec)
    if octave_band_centre ~= false
        oct_filt = octaveFilter(octave_band_centre,"SampleRate",sample_rate);
        ir = oct_filt(ir);
    end
    
    if exist('integration_limit_sec','var')
        integration_limit_samples = integration_limit_sec * sample_rate;
    else
        integration_limit_samples = size(ir,1);
    end

    % Calculate Schroeder decay
    integrated_values_dB(integration_limit_samples:-1:1)=10*log10(cumsum(ir(integration_limit_samples:-1:1).^2)/sum(ir(1:integration_limit_samples).^2));

    time_values_samples = 1:length(integrated_values_dB);
    time_values_seconds = time_values_samples / sample_rate;
end