
octave_cutoff = 1000;

hold on

[ir, fs] = audioread("../../../../../Documents/PhD/Kentish Town Predicted AAES RIRs/ReverberatorRTFactor[1]_LoopGain[-5].wav");
oct_filt = octaveFilter(octave_cutoff,"SampleRate",fs);
ir = oct_filt(ir);

ir = ir / max(abs(ir));

PlotEDC(ir, fs, "-", 3, -50);

function PlotAllSpectrograms()
    num_channels_set = [8 12 16];
    room_nums = [1 2 3];
    alpha_sets = [1 2 3];
    rt_factors = [0 0.5 1 2 4];
    loop_gain_biases_dB = [0 -2 -4 -6];

    % 3 parameters by 3^3 combinations
    combined_param_map = GenerateCombinedParamMap(num_channels_set, room_nums, alpha_sets);
    
    for combined_index = 1:size(combined_param_map, 2)
        num_channels = combined_param_map(1, combined_index);
        room_num = combined_param_map(2, combined_index);
        alpha_set = combined_param_map(3, combined_index);
    
        for rt_factor_index = 1:size(rt_factors, 2)
            rt_factor = rt_factors(rt_factor_index);

            fig = figure("Visible","off");

            loop_gain_bias_dB = 0;

            [ir, fs] = audioread("AAES Modelled IRs KT/ReverberatorRTFactor["+rt_factor+"]_LoopGain["+loop_gain_bias_dB+"].wav");
            PlotSpectrogram(ir, fs, 6);
            title("Modelled");

            if rt_factor == 0.5
                rt_factor = "0_5";
            end

            [ir, fs] = audioread("../../Kentish Town Lab RIRs 20240215/Best Measurements/RT"+rt_factor+"LG0.wav");
            PlotSpectrogram(ir, fs, 6);
            title("Measured");

            saveas(fig, "Plots/KT Spectrograms/Kentish Town AAES Modelled vs Measured RTFactor["+rt_factor+"]","png");
        end
    end
end

function PlotAllEDCs(octave_centre_f)
    oct_filt = octaveFilter;

    if exist('octave_centre_f','var')
        oct_filt = octaveFilter(octave_centre_f);
    end

    num_channels_set = [8 12 16];
    room_nums = [1 2 3];
    alpha_sets = [1 2 3];
    rt_factors = [0 0.5 1 2 4];
    loop_gain_biases_dB = [0 -2 -4 -6];
    
    % 3 parameters by 3^3 combinations
    combined_param_map = GenerateCombinedParamMap(num_channels_set, room_nums, alpha_sets);
    
    for combined_index = 1:size(combined_param_map, 2)
        num_channels = combined_param_map(1, combined_index);
        room_num = combined_param_map(2, combined_index);
        alpha_set = combined_param_map(3, combined_index);
    
        [ir_example, fs] = audioread("AAES Pink Model Data/AAES Modelled IRs/Ch["+num_channels+"] Room["+room_num+"] AlphaSet["+alpha_set+"]/ReverberatorRTFactor[4]_LoopGain[0].wav");

        if exist('octave_centre_f','var')
            ir_example = oct_filt(ir_example);
        end

        room_case_max_ir_length_sec = size(ir_example,1) / fs;

        for rt_factor_index = 1:size(rt_factors, 2)
            rt_factor = rt_factors(rt_factor_index);

            fig = figure("Visible","off");
            hold on

            legend_cells = cell(size(loop_gain_biases_dB,2));

            for loop_gain_index = 1:size(loop_gain_biases_dB, 2)
                loop_gain_bias_dB = loop_gain_biases_dB(loop_gain_index);

                [ir, fs] = audioread("AAES Pink Model Data/AAES Modelled IRs/Ch["+num_channels+"] Room["+room_num+"] AlphaSet["+alpha_set+"]/ReverberatorRTFactor["+rt_factor+"]_LoopGain["+loop_gain_bias_dB+"].wav");

                if exist('octave_centre_f','var')
                    ir = oct_filt(ir);
                end

                PlotEDC(ir, fs, "-", room_case_max_ir_length_sec);

                legend_cells{loop_gain_index} = "LoopGain["+loop_gain_bias_dB+"]";
            end

            [ir, fs] = audioread("AAES Pink Model Data/AAES Modelled IRs/Ch["+num_channels+"] Room["+room_num+"] AlphaSet["+alpha_set+"]/E_1_1.wav");

            if exist('octave_centre_f','var')
                ir = oct_filt(ir);
            end

            PlotEDC(ir, fs, "--", room_case_max_ir_length_sec);

            legend_cells = resize(legend_cells,[size(loop_gain_biases_dB,2) + 1, 1]);
            legend_cells{size(loop_gain_biases_dB,2) + 1} = "Passive Room";
            legend(legend_cells);

            if rt_factor == 0.5
                rt_factor = "0_5";
            end

            saveas(fig, "AAES Pink Model Data/Plots/1kHz EDCs/AAES Model EDCs Ch["+num_channels+"] Room["+room_num+"] AlphaSet["+alpha_set+"] RTFactor["+rt_factor+"]","png");
        end
    end
end