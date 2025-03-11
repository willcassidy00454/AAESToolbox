% Saves .wav files for a tensor of IRs with dimensions (rec pos, src pos, sample pos)
function SaveIRs(irs, sample_rate, bit_depth, output_dir, group_label)
    for rec_pos = 1:size(irs,1)
        for src_pos = 1:size(irs,2)
            filename = group_label + "_R" + rec_pos + "_S" + src_pos + ".wav";

            audiowrite(output_dir + filename, squeeze(irs(rec_pos, src_pos, :)), sample_rate, 'BitsPerSample', bit_depth);
        end
    end
end