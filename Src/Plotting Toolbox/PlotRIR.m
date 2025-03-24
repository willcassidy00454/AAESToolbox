
% octave_cutoff = 1000;

% close all

% hold on

% [ir, fs] = audioread("/Users/willcassidy/Documents/PhD/AAESinMATLAB/Outputs/AAES_RIR_TI.wav");

% [ir, fs] = audioread("/Users/willcassidy/Development/GitHub/AAESToolbox/Active Acoustics Review/AAES Receiver RIRs/AAES Condition 9/ReceiverRIR.wav");
% 
% edc_deriv = zeros(7, ceil(length(ir) / 20));
% edc_deriv(1,:) = GetEDCDerivative(ir, fs, 200, 5);
% edc_deriv(2,:) = GetEDCDerivative(ir, fs, 400, 5);
% edc_deriv(3,:) = GetEDCDerivative(ir, fs, 800, 5);
% edc_deriv(4,:) = GetEDCDerivative(ir, fs, 1600, 5);
% edc_deriv(5,:) = GetEDCDerivative(ir, fs, 3200, 5);
% edc_deriv(6,:) = GetEDCDerivative(ir, fs, 6400, 5);
% edc_deriv(7,:) = GetEDCDerivative(ir, fs, 12800, 5);
% 
% plt = waterfall(edc_deriv);
% zlim([-100 0]);
% set(plt, 'FaceColor', 'k');

% PlotEnergyDecayRelief(ir, fs);

% PlotSpectrogram(ir, fs);

% PlotSpectrogramsForLoopGains("Active Acoustics Review/AAES Receiver RIRs/Condition 1/", [-50 -4 -2 -0.5 0], 4);

% PlotSpectAndEDCForConditions("Active Acoustics Review/AAES Receiver RIRs/", 1:16, 3);
% tiledlayout(1,4);
% nexttile;
% for i = 14
%     figure;
%     PlotEnergyDecayReliefForCondition("Active Acoustics Review/AAES Receiver RIRs/", i);
% end
% nexttile;
% PlotEnergyDecayReliefForCondition("Active Acoustics Review/AAES Receiver RIRs/", 7);
% nexttile;
% PlotEnergyDecayReliefForCondition("Active Acoustics Review/AAES Receiver RIRs/", 8);
% nexttile;
% PlotEnergyDecayReliefForCondition("Active Acoustics Review/AAES Receiver RIRs/", 9);


% figure
% nexttile
% for i=1:4
% hold on
% [ir, fs] = audioread("AbsorbCoeffsTest/ConcertHall.wav");
% [ir, fs] = audioread("Reverberators/Reverberator 3/X_R1_S1.wav");
% [irs, fs] = audioread("Audio Data/AAES Receiver RIRs/AAES Room 1 Absorption 1 RT 1 Loop Gain -1 Filter 1 Routing 1/ReceiverRIR.wav");
% ir = irs(:,4);
% [ir, fs] = audioread("AAESinMATLAB/Outputs/AAES_Inline_TV_+2.wav");
% [ir, fs] = audioread("Active Acoustics Review/Generated AAES RIRs/Room Condition 1/G_R1_S1.wav");
% tiledlayout(3,1);
hold on
[ir, fs] = audioread("Active Acoustics Review/AAES Receiver RIRs/AAES Condition 3/ReceiverRIR.wav");
% ir1 = zeros(48000 * 4, 1);
% [ir, fs] = audioread("Active Acoustics Review/Generated AAES RIRs/Room Condition 1/E_R1_S1.wav");
% [ir, fs] = audioread("Active Acoustics Review/AAES Receiver RIRs/Pilsen.wav");

% delay_matrix = readmatrix("Active Acoustics Review/Directivities/delay_matrix_room_10.dat");
% PlotMatrixDRR("Active Acoustics Review/Generated AAES RIRs/Room Condition 9/", "H", 8, 8, "b. Omni Mics", delay_matrix);
% plot(ir);
% ir1(1:length(ir)) = ir;
% disp(mean(ir));

PlotEDC(ir, fs, false, "-.", 2.5);
% PlotSpectrogram(ir, fs, 1.2, true);

% PlotEDCDerivative(ir, fs, 3000, 1.8);
% disp(FindT30(ir, fs, 125));
% disp(FindT30(ir, fs, 250));
% disp(FindT30(ir, fs, 500));
% disp(FindT30(ir, fs, 1000));
% disp(FindT30(ir, fs, 2000));
% disp(FindT30(ir, fs, 4000));
% disp(FindT30(ir, fs, 8000));
% end

% For single EDC plots:
set(gcf, "position", [300 300 600 500]);


% For triple-stacked figures:
% set(gcf, "position", [300 000 550 900]);


% For two side-by-side figures:
% set(gcf,'position',[300,300,1000,400]);


% tiledlayout(1,3);
% nexttile;
% PlotMatrixDRR("Active Acoustics Review/Generated AAES RIRs/Room Condition 1/","H",16,16,"DRR: All Mics, All LS");
% nexttile;
% PlotMatrixDRR("Active Acoustics Review/Generated AAES RIRs/Room Condition 6/","H",8,8,"DRR: Central Omni Mics, Wall LS");
% nexttile;
% PlotMatrixDRR("Active Acoustics Review/Generated AAES RIRs/Room Condition 7/","H",8,8,"DRR: Wall Cardioid Mics, Wall LS");

% [ir, fs] = audioread("Active Acoustics Review/AAES Receiver RIRs/Pilsen.wav");

% for rec = 10:19
% [ir, fs] = audioread("/Users/willcassidy/Documents/PhD/Pilsen/Full_TP_Matrix3/E013_R008_M01.wav");
% [ir, fs] = audioread("Active Acoustics Review/Generated AAES RIRs/Room Condition 1/E_1_1.wav");
% ir = ir / max(abs(ir));
% PlotSpectrogram(ir, fs, 1);
% % end
% title("Estimated Absorption");
% figure
% PlotHeatmapForRouting("AAESDatasetGenerator/Routings/routing_4.dat", 16, 16, true);

% [sine, ~] = audioread("Active Acoustics Review/Sine.wav");
% [tvrev, ~] = audioread("Active Acoustics Review/TVReverb.wav");
% [ir, fs] = audioread("Active Acoustics Review/Reverberators/Reverberator 6/X_1_1.wav");
% PlotSpectrogram(conv(ir,sine(:,1)),fs,10);
% title("TI");

% sound(conv(ir,sine(:,1)), 48000);

% [ir, fs] = audioread("Active Acoustics Review/Reverberators/Reverberator 9/X_1_1.wav");
% PlotSpectrogram(conv(ir,sine(:,1)),fs,10);
% title("TV");

% sound(conv(ir,sine(:,1)), 48000);

% PlotSpectrogram(conv(tvrev,sine(:,1)),fs,10);
% title("TV");

% sound(conv(tvrev,sine(:,1)), 48000);

% % % Use this for the spherical harmonics test:
% for plot_idx = 1:9
%     % figure
%     % PlotIRs("Active Acoustics Review/Generated AAES RIRs SH Test/Azimuth/", 8, plot_idx);
%     % set(gcf,'position',[300 * plot_idx, 400, 300, 800]);
%     PlotSHDirectivity("Active Acoustics Review/Generated AAES RIRs SH Test/Azimuth/", 8, plot_idx, 400);
% end

% % % Use this for directivity test:
% for plot_idx = 1:9
%     nexttile
%     set(gcf,'position',[300, 400, 300, 800]);
%     [ir, fs] = audioread("Active Acoustics Review/Generated AAES RIRs/Src Directivity Test/E_R1_S"+plot_idx+".wav");
%     % [ir, fs] = audioread("Active Acoustics Review/Generated AAES RIRs/Rec Directivity Test/E_R"+plot_idx+"_S1.wav");
%     plot(ir);
%     xlim([0 1000]);
%     ylim([-0.025 0.025]);
% end

% oct_filt = octaveFilter(octave_cutoff,"SampleRate",fs);
% ir = oct_filt(ir);

% ir = ir / max(abs(ir));

% PlotEDC(ir, fs, "-", 10, -50);
% PlotSpectrogram(ir,fs,5);

function PlotSHDirectivity(read_dir, num_indices, channel_to_plot, trunc_length_samples)
    irs = zeros(num_indices,trunc_length_samples);
    maxima = zeros(num_indices);
    thetas = zeros(num_indices);

    for index = 1:num_indices
        [ir, ~] = audioread(read_dir + "Src Position "+index+"_R1_S1.wav");
        irs(index,:) = ir(1:trunc_length_samples, channel_to_plot);

        maxima(index) = max(abs(irs(index,:)),[],"all");
        thetas(index) = (index - 1) * pi / 4;
        % plot(irs(index,:));
        % ylim([-0.5 0.5]);
    end

    maxima(7) = 2 * maxima(7) / 3;
    maxima = 2 * maxima;

    nexttile;
    polarplot(thetas,maxima);
    rlim([0 1.0])
    title("Azimuth - Spherical Harmonic "+channel_to_plot)
end

function PlotIRs(read_dir, num_indexes_to_append, channel_to_plot)
    tiledlayout(num_indexes_to_append,1);

    for index = 1:num_indexes_to_append
        nexttile;
        [ir, ~] = audioread(read_dir + "Src Position " + index + "_R1_S1.wav");
        % plot (1:length(ir), 20*log10(abs(ir)));
        plot(ir(:,channel_to_plot));
        % PlotEDC(ir, 48000);
        title("Index " + index);
        xlim([0, 1000]);
        ylim([-0.5 0.5]);
    end
end

function PlotSpectrogramsForLoopGains(read_dir, loop_gains, length_sec)
    for loop_gain = loop_gains
        [ir, fs] = audioread(read_dir + "LoopGain[" + loop_gain + "].wav");
        PlotSpectrogram(ir,fs,length_sec);
        title("Loop Gain = "+loop_gain+"dB");
    end
end

function PlotEnergyDecayReliefForCondition(read_dir, condition_index)
    [ir, fs] = audioread(read_dir + "AAES Condition "+condition_index+"/ReceiverRIR.wav");
    PlotEnergyDecayRelief(ir,fs);
    title("Condition "+condition_index);
end

function PlotEnergyDecayRelief(ir,fs)
    ir_trunc_length = 5;
    plot_time_length_ratio = ir_trunc_length / (length(ir) / fs);
    % ir = zeros(ir_trunc_length * fs,1);
    % 
    % num_nonzero_samples = min(length(ir_raw), length(ir));
    % ir(1:num_nonzero_samples) = ir_raw(1:num_nonzero_samples);

    z_min = -220;
    [p,f] = pspectrum(ir,fs,"spectrogram","FrequencyLimits",[20 20000],"OverlapPercent",80,"FrequencyResolution", 40,"MinThreshold",z_min);
    
    num_freqs_to_plot = 10;
    [energy_decay, log_f] = calcEDR(p,f,num_freqs_to_plot);

    % surf(energy_decay);
    plt = waterfall(energy_decay);
    set(plt, 'FaceColor', 'k');

    freq_step = num_freqs_to_plot / 10;
    yticks(1:freq_step:num_freqs_to_plot);
    yticklabels(round(log_f(1:freq_step:num_freqs_to_plot)));
    ylabel("Frequency / Hz");

    max_time_bins = size(p,2) * plot_time_length_ratio;
    xticks(1:(max_time_bins/8):max_time_bins);
    xticklabels(0:(ir_trunc_length/8):ir_trunc_length);
    xlabel("Time / s");
    xlim([0 max_time_bins])

    % zticklabels([]);
    zlabel("Energy / dB");

    view(121, 38); %view (45, 40);% view(114, 8); %
end

function PlotSpectAndEDCForConditions(read_dir, condition_indices, length_sec)
    [passive_ir_base, fs] = audioread("Active Acoustics Review/Generated AAES RIRs/Room Condition 1/E_1_1.wav");
    [passive_ir_cond11, ~] = audioread("Active Acoustics Review/Generated AAES RIRs/Room Condition 4/E_1_1.wav");
    [passive_ir_cond12, ~] = audioread("Active Acoustics Review/Generated AAES RIRs/Room Condition 5/E_1_1.wav");

    for condition = condition_indices
        [ir, ~] = audioread(read_dir + "AAES Condition " + condition + "/ReceiverRIR.wav");
        PlotSpectrogram(ir,fs,length_sec);
        title("Condition " + condition);

        nexttile
        if condition == 11
            PlotEDC(passive_ir_cond11, fs, false, "--", length_sec, -70);
        elseif condition == 12
            PlotEDC(passive_ir_cond12, fs, false, "--", length_sec, -70);
        else
            PlotEDC(passive_ir_base, fs, false, "--", length_sec, -70);
        end
        hold on
        PlotEDC(ir, fs, false, "-", length_sec, -70);
        hold off
        title("Condition " + condition);
    end

    PlotSpectrogram(passive_ir_base,fs,length_sec);
    title("Passive Room (base absorption)");
end

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

function PlotHeatmapForRouting(routing_dir, num_rows, num_cols, convert_to_dB)
    nexttile

    routing = readmatrix(routing_dir);% zeros(num_rows, num_cols);

    % for row = 1:num_rows
    %     for col = 1:num_cols
    %         [routing(row, col), ~] = audioread(routing_dir + "X_R"+row+"_S"+col+".wav");
    %     end
    % end

    if (convert_to_dB)
        if ~isempty(find(routing == 0))
            routing(find(routing == 0)) = 0.001;
        end

        routing = 20 * log10(abs(routing));
    end

    heatmap(routing, "Colormap", parula, "CellLabelColor", "none");

    if (convert_to_dB)
        title("Routing Matrix Magnitude / dB");
    else
        title("Routing Matrix Gain");
    end

    xlabel("Microphones");
    ylabel("Loudspeakers");

    if (convert_to_dB)
        clim([-60 0]);
    else
        clim([-1 1]);
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

function PlotMatrixDRR(read_dir, matrix_prefix, num_rows, num_cols, plot_title, src_rec_delay_matrix)
    nexttile
    
    if ~exist("src_rec_delay_matrix", "var")
        [matrix_drr, delay_matrix] = GetMatrixDRR(read_dir, ...
                                matrix_prefix, ...
                                num_rows, ...
                                num_cols);
    else
        [matrix_drr, ~] = GetMatrixDRR(read_dir, ...
                                matrix_prefix, ...
                                num_rows, ...
                                num_cols, ...
                                [], ...
                                src_rec_delay_matrix);
    end

    heatmap(matrix_drr, "Colormap", parula, "CellLabelColor", "none");

    clim([-60, 0]);

    xlabel("Microphones");
    ylabel("Loudspeakers");

    title(plot_title);
end