function ConvertMatToWavs(mat_file_dir, wav_output_dir)
    audio_matrix = load(mat_file_dir).IRs;
    mkdir(wav_output_dir);

    audio_matrix = audio_matrix / max(abs(audio_matrix),[],"all");

    for row = 1:size(audio_matrix, 2)
        for col = 1:size(audio_matrix, 3)
            audiowrite(wav_output_dir + "X_" + row + "_" + col + ".wav", audio_matrix(:,row,col), 48000, "BitsPerSample", 24);
        end
    end
end