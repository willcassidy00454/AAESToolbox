% This script generates AAES IRs for a set of rooms with varying
% dimensions, absorption coefficients and AAES channel counts.

%% User Parameters

% General audio parameters
sample_rate = 48000;
bit_depth = 24;

% Room dimensions in metres [x y z; x y z; etc.]
room_dims = [5.7    7.35    2.5;
             8.74   17      5.5;
             19.52  30.83   15];

absorptions_dir = "Example Absorption Coefficients/";
transducer_coords_dir = "Example Transducer Coordinates/";
output_parent_dir = "Automated RIRs/";

uses_parallel_processing = true; % Change this to false if you don't have access to the Parallel Computing Toolbox

%% Generation

% Wall absorption coefficients
% Read "Example Absorption Coefficients/Absorption Coefficients Info.txt" for
% more info
alphas = zeros(7, 6, 3);

for alpha_set = 1:3
    alphas(:, :, alpha_set) = readmatrix(absorptions_dir + "alpha_set_" + alpha_set + ".dat");
end

% Currently assumes there is only one source and receiver per room
% These use the dimensions (room number, axis)
% [x y z; x y z; x y z]
% The sources are 2/3 along the x and y axes and 1.2m high
src_sets = [2*room_dims(1,1)/3 2*room_dims(1,2)/3 1.2
            2*room_dims(2,1)/3 2*room_dims(2,2)/3 1.2
            2*room_dims(3,1)/3 2*room_dims(3,2)/3 1.2];
% The receivers are 1/3 along the x and y axes and 1.2m high
rec_sets = [room_dims(1,1)/3 room_dims(1,2)/3 1.2
            room_dims(2,1)/3 room_dims(2,2)/3 1.2
            room_dims(3,1)/3 room_dims(3,2)/3 1.2];

% ls_sets_... are arranged with dimensions: (ls number, axis, room number)
% These aren't concatenated due to their dimensions
ls_sets_8ch = zeros(8,3,3);
ls_sets_12ch = zeros(12,3,3);
ls_sets_16ch = zeros(16,3,3);

% mic_sets_... are arranged with dimensions: (mic number, axis, room number)
mic_sets_8ch = zeros(8,3,3);
mic_sets_12ch = zeros(12,3,3);
mic_sets_16ch = zeros(16,3,3);

% Load coordinates
for room_num = 1:3
    ls_sets_8ch(:,:,room_num) = readmatrix(transducer_coords_dir + "ls_positions_room_" + room_num + "_8ch.dat");
    ls_sets_12ch(:,:,room_num) = readmatrix(transducer_coords_dir + "ls_positions_room_" + room_num + "_12ch.dat");
    ls_sets_16ch(:,:,room_num) = readmatrix(transducer_coords_dir + "ls_positions_room_" + room_num + "_16ch.dat");
    
    mic_sets_8ch(:,:,room_num) = readmatrix(transducer_coords_dir + "mic_positions_room_" + room_num + "_8ch.dat");
    mic_sets_12ch(:,:,room_num) = readmatrix(transducer_coords_dir + "mic_positions_room_" + room_num + "_12ch.dat");
    mic_sets_16ch(:,:,room_num) = readmatrix(transducer_coords_dir + "mic_positions_room_" + room_num + "_16ch.dat");
end

if uses_parallel_processing
    if isempty(gcp('nocreate'))
        parpool('Processes');
    end
    
    % This currently assumes channel counts of 8, 12 and 16
    parfor channel_count_index = 1:3
        disp(4 + 4 * channel_count_index + "-channel AAES...");
        
        if channel_count_index == 1
            GenerateRIRsParallel(room_dims, alphas, src_sets, rec_sets, ls_sets_8ch, mic_sets_8ch, sample_rate, bit_depth);
        elseif channel_count_index == 2
            GenerateRIRsParallel(room_dims, alphas, src_sets, rec_sets, ls_sets_12ch, mic_sets_12ch, sample_rate, bit_depth);
        else
            GenerateRIRsParallel(room_dims, alphas, src_sets, rec_sets, ls_sets_16ch, mic_sets_16ch, sample_rate, bit_depth);
        end
    end
    
    delete(gcp('nocreate'));
else
    % This currently assumes channel counts of 8, 12 and 16
    for channel_count_index = 1:3
        disp(4 + 4 * channel_count_index + "-channel AAES...");
        
        if channel_count_index == 1
            GenerateRIRs(room_dims, alphas, src_sets, rec_sets, ls_sets_8ch, mic_sets_8ch, sample_rate, bit_depth);
        elseif channel_count_index == 2
            GenerateRIRs(room_dims, alphas, src_sets, rec_sets, ls_sets_12ch, mic_sets_12ch, sample_rate, bit_depth);
        else
            GenerateRIRs(room_dims, alphas, src_sets, rec_sets, ls_sets_16ch, mic_sets_16ch, sample_rate, bit_depth);
        end
    end
end

function GenerateRIRs(output_parent_dir, rooms, alphas, src_sets, rec_sets, ls_sets, mic_sets, sample_rate, bit_depth)
    for room_num = 1:size(rooms, 1)
        disp("Room " + room_num + "...");
    
        current_config = RoomWithAAES(rooms(room_num,:), alphas, src_sets(room_num,:), rec_sets(room_num,:), ls_sets(:,:,room_num), mic_sets(:,:,room_num), sample_rate, bit_depth);
        current_config.GenerateSystemIRs(output_parent_dir);
    end
end

function GenerateRIRsParallel(output_parent_dir, rooms, alphas, src_sets, rec_sets, ls_sets, mic_sets, sample_rate, bit_depth)
    parfor room_num = 1:size(rooms, 1)
        disp("Room " + room_num + "...");
    
        current_config = RoomWithAAES(rooms(room_num,:), alphas, src_sets(room_num,:), rec_sets(room_num,:), ls_sets(:,:,room_num), mic_sets(:,:,room_num), sample_rate, bit_depth);
        current_config.GenerateSystemIRs(output_parent_dir);
    end
end
