function edc_derivative = GetEDCDerivative(ir,fs,octave_band_centre,x_max)
    if ~exist("x_max", "var")
        x_max = 4;
    end
    
    [edc, ~] = GetEDC(ir,fs,octave_band_centre);
    skip = 20;
    
    edc_derivative = zeros(ceil(length(edc) / skip), 1);
    derivative_pos = 1;
    
    for sample_index = 1:skip:length(edc)-skip
        % Get gradient of line from index to index + 10
        period = skip / fs;
        gradient = (edc(sample_index + skip) - edc(sample_index)) / period;
        edc_derivative(derivative_pos) = gradient;
        derivative_pos = derivative_pos + 1;
    end

    edc_derivative = smooth(edc_derivative, 900);
    
    time_values = ((1:length(edc_derivative)) * skip) / fs;
end