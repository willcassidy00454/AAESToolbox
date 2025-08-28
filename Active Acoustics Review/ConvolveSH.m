programme_item_read_dir = "Active Acoustics Review/Speech.wav";
[ir_speech, ~] = audioread(programme_item_read_dir,[1 (5 * 48000)]);

for pos_index = 1:8
    sh_read_dir = "Active Acoustics Review/Generated AAES RIRs SH Test/Azimuth/";
    [ir_sh, fs] = audioread(sh_read_dir + "Src Position " + pos_index + "_1.wav");
    output = zeros(16, length(ir_sh) + length(ir_speech) - 1);
    
    for channel = 1:16
        output(channel,:) = conv(ir_speech, ir_sh(:,channel));
    end
    
    output = output / max(abs(output),[],"all");

    audiowrite(sh_read_dir + "Speech_"+pos_index+".wav", output', fs);
end