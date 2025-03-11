% Saves num_src of multichannel .wav files each with channel count =
% num_rec
function SaveIRsMultichannelReceivers(irs, sample_rate, bit_depth, output_dir, group_label)
    for src_pos = 1:size(irs,2)
        filename = group_label + "_R1_S" + src_pos + ".wav";
        output = squeeze(irs(:, src_pos, :))';

        audiowrite(output_dir + filename, output, sample_rate, 'BitsPerSample', bit_depth);
    end
end