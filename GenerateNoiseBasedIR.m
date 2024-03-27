% This script creates a set of Gaussian pink noise IRs which decay
% exponentially given a desired RT. This will be saved at the specified
% sample rate as a mono waveform file.
% These are arranged diagonally, assuming a 1-to-1 reverberator matrix.

% This matrix maps reverb times to their specific cases
% room num x alpha set
reverb_time_map = [1.08];%[0.25 0.51 1.14
                   % 0.54 1.08 2.11
                   % 1.19 2.30 3.90];

rt_factors = [1.5 2.5 3 3.5];

sample_rate = 48000;
num_channels = 16;
bit_depth = 24;

cutoff = 20;
[zhi,phi,khi] = butter(8,2*cutoff/sample_rate,"high");
sos_hpf = zp2sos(zhi,phi,khi);

for room_num = 1:size(reverb_time_map, 1)
    for alpha_set = 1:size(reverb_time_map, 2)
        for rt_factor_index = 1:size(rt_factors, 2)
            rt_factor = rt_factors(rt_factor_index);
            output_dir = "/Users/willcassidy/Documents/PhD/MATLAB/AAES Model/Pink Reverberator IRs/Decaying Noise Ch["+num_channels+"] Room[2] AlphaSet[2] RTFactor["+rt_factor+"]/";
            % output_dir = "/Users/willcassidy/Documents/PhD/MATLAB/AAES Model/Pink Reverberator IRs/Decaying Noise Ch["+num_channels+"] Room["+room_num+"] AlphaSet["+alpha_set+"] RTFactor["+rt_factor+"]/";
            mkdir(output_dir);

            GenerateAndSaveIR(reverb_time_map(room_num, alpha_set) * rt_factor, sample_rate, output_dir, num_channels, bit_depth, sos_hpf);
        end
    end
end

function GenerateAndSaveIR(reverb_time, sample_rate, file_directory, num_channels, bit_depth, sos_hpf)
    % reverb_time = 1;
    tail_length_seconds = reverb_time; % Time after 60dB of decay
    % bit_depth = 24;
    % file_directory = "/Users/willcassidy/Documents/PhD/MATLAB/Reverberator IRs/Decaying Noise Ch[8] Room[1] AlphaSet[3] SampleRate[48000]/";
    % num_channels = 8;
    
    % Saves a diagonal of uncorrelated noise (x_1_1, X_2_2 etc.)
    for channel = 1:num_channels
        ir = GenerateIR(reverb_time, tail_length_seconds, sample_rate, sos_hpf);
        ir = ir / max(abs(ir)); % Normalise
        
        SaveIR(ir, file_directory + "X_" + channel + "_" + channel + ".wav", sample_rate, bit_depth);
    end
end

function y = GetExpDecayCoeff(time_seconds, reverb_time)
    y = exp((time_seconds * log(power(10, -3))) / reverb_time);
end

function ir = GenerateIR(reverb_time, tail_length_seconds, sample_rate, sos_hpf)
    if reverb_time ~= 0
        num_samples = floor((reverb_time + tail_length_seconds) * sample_rate);
        ir = pinknoise(num_samples, 'double'); %randn(1, num_samples, 'double');
        ir = ir / max(abs(ir));
        ir = ir';
    
        decay_timesteps = 0:1/sample_rate:(num_samples-1)/sample_rate;
        decay_coefficients = GetExpDecayCoeff(decay_timesteps, reverb_time);

        ir = ir .* decay_coefficients;

        ir = sosfilt(sos_hpf, ir);
    else
        ir = 1; % If reverb time is 0, output an impulse
    end
end

function SaveIR(ir, filename, sample_rate, bit_depth)
    audiowrite(filename, ir, sample_rate, 'BitsPerSample', bit_depth);
end