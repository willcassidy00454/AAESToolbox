% This script generates AAES IRs for a set of rooms with varying
% dimensions, absorption coefficients and AAES channel counts.

close all
clear all

%% User Parameters

% General audio parameters
sample_rate = 48000;
bit_depth = 24;

absorptions_dir = "Active Acoustics Review/Absorption Coefficients/";
coords_dir = "Active Acoustics Review/Directivity Test/";
rotations_dir = "Active Acoustics Review/Directivity Test/";
directivities_dir = "Active Acoustics Review/Directivity Test/";
output_dir = "Active Acoustics Review/Generated AAES RIRs/";

num_absorption_sets = 1;

%% Generation

room_dims = readmatrix("Active Acoustics Review/Room Dimensions/room_dimensions.dat");
alphas = readmatrix(absorptions_dir + "absorption_coeffs_1.dat");

% Src Test:
src_coords = readmatrix(coords_dir + "test_coords.dat");
rec_coords = readmatrix(coords_dir + "static_coords.dat");

src_rotations = readmatrix(rotations_dir + "test_rotations.dat");
rec_rotations = readmatrix(rotations_dir + "static_rotations.dat");

src_directivities = string(readcell(directivities_dir + "test_directivities.csv"));
rec_directivities = string(readcell(directivities_dir + "static_directivities.csv"));

ir_dir = output_dir + "Src Directivity Test/";
source_is_rotating = true;

% Rec Test:
% src_coords = readmatrix(coords_dir + "static_coords.dat");
% rec_coords = readmatrix(coords_dir + "test_coords.dat");
% 
% src_rotations = readmatrix(rotations_dir + "static_rotations.dat");
% rec_rotations = readmatrix(rotations_dir + "test_rotations.dat");
% 
% src_directivities = string(readcell(directivities_dir + "static_directivities.csv"));
% rec_directivities = string(readcell(directivities_dir + "test_directivities.csv"));
% 
% ir_dir = output_dir + "Rec Directivity Test/";
% source_is_rotating = false;

current_config = RoomWithAAES(room_dims, ...
    alphas, ...
    src_coords, ...
    rec_coords, ...
    [], ...
    [], ...
    src_rotations, ...
    rec_rotations, ...
    [], ...
    [], ...
    src_directivities, ...
    rec_directivities, ...
    [], ...
    [], ...
    sample_rate, ...
    bit_depth);
current_config.GenerateSystemIRs(ir_dir, true);

PlotDirectivity(ir_dir, 1, 17, 300, source_is_rotating, "Azimuth");
PlotDirectivity(ir_dir, 18, 17, 300, source_is_rotating, "Elevation");

function PlotDirectivity(read_dir, start_index, num_indices, trunc_length_samples, source_is_rotating, label)
    irs = zeros(num_indices, trunc_length_samples);
    maxima = zeros(num_indices, 1);
    thetas = zeros(num_indices, 1);

    for index = 1:num_indices
        if (source_is_rotating)
            [ir, ~] = audioread(read_dir + "E_R1_S"+(start_index + index - 1)+".wav");
        else
            [ir, ~] = audioread(read_dir + "E_R"+(start_index + index - 1)+"_S1.wav");
        end

        irs(index,:) = ir(1:trunc_length_samples);

        maxima(index) = max(abs(irs(index,:)),[],"all");
        thetas(index) = (index-1) * 2 * pi / (num_indices - 1);
    end

    nexttile
    polarplot(thetas, maxima);
    rlim([0 1]);
    title(label);
end