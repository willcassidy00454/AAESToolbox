% Returns T_30 of an IR in the specified octave band
% This is based on the gradient calculated between -5 dB and -35 dB
function t30 = FindT30(ir, fs, band_centre, bandwidth_mode)
    if ~exist("bandwidth_mode", "var")
        bandwidth_mode = "1 octave";
    end

    if ~exist("band_centre", "var")
        edc_dB = GetEDC(ir, fs);
    else
        edc_dB = GetEDC(ir, fs, band_centre, bandwidth_mode);
    end

    minus_5_index = find(edc_dB <= -5, 1);
    minus_35_index = find(edc_dB <= -35, 1);
    sampling_period = 1 / fs;
    minus_5_time = minus_5_index * sampling_period;
    minus_35_time = minus_35_index * sampling_period;
    gradient = -30 / (minus_35_time - minus_5_time);
    t30 = -60 / gradient;
end