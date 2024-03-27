num_channels_set = [8];%[8 12 16];
room_nums = [1];%[1 2 3];
alpha_sets = [1];% 3];
loop_gain_biases_dB = [0];%-2 -4 -6];

output_directory = "AAES Reversal Test/";

% 3 parameters by 3^3 combinations
combined_param_map = GenerateCombinedParamMap(num_channels_set, room_nums, alpha_sets);

% if isempty(gcp('nocreate'))
%     parpool('Processes');
% end

%parfor
for combined_index = 1:size(combined_param_map, 2)
    num_channels = combined_param_map(1, combined_index);
    room_num = combined_param_map(2, combined_index);
    alpha_set = combined_param_map(3, combined_index);
    % num_channels = 16;
    % room_num = 2;
    % alpha_set = 2;

    GenerateAAESIRs(num_channels, room_num, alpha_set, loop_gain_biases_dB, output_directory);
end

% delete(gcp('nocreate'));