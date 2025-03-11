function PlotEDCDerivative(ir,fs,octave_band_centre,x_max)
    if ~exist("x_max", "var")
        x_max = 4;
    end
    if ~exist("octave_band_centre", "var")
        octave_band_centre = false;
    end
    
    edc = GetEDC(ir,fs,octave_band_centre);
    edc_derivative = zeros(length(edc), 1);
    window_size = 4000;

    for sample_index = 1:length(edc)-window_size
        period = 1 / fs;
        % gradient = (edc(sample_index + window_size) - edc(sample_index)) / period;
        sample_range = sample_index:sample_index + window_size;
        lin_reg = polyfit(sample_range * period, edc(sample_range), 1);
        gradient = lin_reg(1);
        edc_derivative(sample_index) = gradient;
    end
    
    time_values = (1:length(edc_derivative)) / fs;
    plot(time_values, edc_derivative);
    title("Derivative");
    ylabel("EDC Decay Rate (dB/s)");
    xlabel("Time (s)");
    xlim([0 x_max]);
    ylim([-100 0]);
    set(gcf,'position',[400,2000,700,250]);
    grid on;
end