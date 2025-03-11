%% User Parameters

rir_dir = "Active Acoustics Review/Generated AAES RIRs/AAES IRs Ch[16x16] Room[35 18 8]/";
reverberator_dir = "Active Acoustics Review/Indentity Reverberator 16ch/";
output_dir = "Active Acoustics Review/AAES Receiver RIRs/";

% num_channels_set = [16];%[8 12 16];
% room_nums = [2];%[1 2 3];
% alpha_sets = [1 2];% 3];
loop_gain_biases_dB = [-20];%-2 -4 -6];

uses_parallel_processing = true; % Change this to false if you don't have access to the Parallel Computing Toolbox

%% Generation

% 3 parameters by 3^3 combinations
% combined_param_map = GenerateCombinedParamMap(num_channels_set, room_nums, alpha_sets);

if uses_parallel_processing
    % if isempty(gcp('nocreate'))
    %     parpool('Processes');
    % end
    % 
    % parfor combined_index = 1:size(combined_param_map, 2)
    %     num_channels = combined_param_map(1, combined_index);
    %     room_num = combined_param_map(2, combined_index);
    %     alpha_set = combined_param_map(3, combined_index);
        
        % Currently generates a square AAES
        GenerateAAESIRs(rir_dir, reverberator_dir, output_dir, loop_gain_biases_dB, 16, 16);
    % end
    % 
    % delete(gcp('nocreate'));
else
    % for combined_index = 1:size(combined_param_map, 2)
    %     num_channels = combined_param_map(1, combined_index);
    %     room_num = combined_param_map(2, combined_index);
    %     alpha_set = combined_param_map(3, combined_index);

        % Currently generates a square AAES
        GenerateAAESIRs(rir_dir, reverberator_dir, loop_gain_biases_dB, loop_gain_biases_dB, 16, 16);
    % end
end