% This script simulates multiple IRs relating to a specified rectangular
% room and saves them into a folder in the format required by
% "MultichannelAAESModel.m". This script does not generate the input
% matrix, U, or the reverberator matrix, X.

% Inputs                Dimensions
%------------------------------------------
% room_dim              1x3
% alphas                7x6
% src_positions         Nx3
% rec_positions         Mx3
% ls_positions          Kx3
% mic_positions         Lx3
% src_rotations         Nx2
% rec_rotations         Mx2
% ls_rotations          Kx2
% mic_rotations         Lx2
% src_directivities     Nx1
% rec_directivities     Mx1
% ls_directivities      Kx1
% mic_directivities     Lx1
% sample_rate           1x1
% output_dir            1x1

function GenerateAKToolsRIRs3rdOrderOutput(room_dim, alphas, src_positions, rec_positions, ls_positions, mic_positions, src_rotations, rec_rotations, ls_rotations, mic_rotations, src_directivities, rec_directivities, ls_directivities, mic_directivities, sample_rate, output_dir, bit_depth, should_high_pass, should_normalise)
    if ~exist('should_normalise','var')
        should_normalise = true;
    end
    
    %% Generate all IRs

    u = 1; % Source is a single dirac delta

    % This assumes one receiver, but with 3rd order spherical harmonics
    for third_order_sh_output_index = 1:16
        e(third_order_sh_output_index,1,:) = GenerateSrcToRecIRs(src_positions, rec_positions, src_rotations, rec_rotations, src_directivities, rec_directivities, room_dim, alphas, sample_rate, "E", should_high_pass, third_order_sh_output_index); % Sources to Receivers
        f(third_order_sh_output_index,:,:) = GenerateSrcToRecIRs(ls_positions, rec_positions, ls_rotations, rec_rotations, ls_directivities, rec_directivities, room_dim, alphas, sample_rate, "F", should_high_pass, third_order_sh_output_index); % AAES to Receivers
    end

    g = GenerateSrcToRecIRs(src_positions, mic_positions, src_rotations, mic_rotations, src_directivities, mic_directivities, room_dim, alphas, sample_rate, "G", should_high_pass); % Sources to AAES
    h = GenerateSrcToRecIRs(ls_positions, mic_positions, ls_rotations, mic_rotations, ls_directivities, mic_directivities, room_dim, alphas, sample_rate, "H", should_high_pass); % AAES to AAES

    %% Normalise all IRs as one group

    if should_normalise
        disp("Normalising IRs...");

        max_e = max(abs(e),[],"all");
        max_f = max(abs(f),[],"all");
        max_g = max(abs(g),[],"all");
        max_h = max(abs(h),[],"all");

        max_amplitude = max([max_e max_f max_g max_h],[],"all");

        e = e / max_amplitude;
        f = f / max_amplitude;
        g = g / max_amplitude;
        h = h / max_amplitude;
    end

    %% Save all IR files

    disp("Saving IRs...");

    SaveIRs(u, sample_rate, bit_depth, output_dir, "U"); % Untested - this should output a single dirac delta
    SaveIRsMultichannelReceivers(e, sample_rate, bit_depth, output_dir, "E");
    SaveIRsMultichannelReceivers(f, sample_rate, bit_depth, output_dir, "F");
    SaveIRs(g, sample_rate, bit_depth, output_dir, "G");
    SaveIRs(h, sample_rate, bit_depth, output_dir, "H");

    %% Write info to log txt

    writelines("Sample rate: " + sample_rate, output_dir + "README.txt", WriteMode="overwrite");
    writelines("Bit depth: " + bit_depth, output_dir + "README.txt", WriteMode="append");
    writelines("Room dimensions (m):", output_dir + "README.txt", WriteMode="append");
    writetable(array2table(room_dim), output_dir + "README.txt", WriteMode="append");
    writelines("Absorption Coefficients ([x1, x2, y1, y2, z1, z2] across frequency):", output_dir + "README.txt", WriteMode="append");
    writetable(array2table(alphas), output_dir + "README.txt", WriteMode="append");
    writelines("The following co-ordinates are arranged in columns, where each column is X, Y and Z, respectively:", output_dir + "README.txt", WriteMode="append");
    writelines("Source positions:", output_dir + "README.txt", WriteMode="append");
    writetable(array2table(src_positions), output_dir + "README.txt", WriteMode="append");
    writelines("Receiver positions:", output_dir + "README.txt", WriteMode="append");
    writetable(array2table(rec_positions), output_dir + "README.txt", WriteMode="append");
    writelines("AAES loudspeaker positions:", output_dir + "README.txt", WriteMode="append");
    writetable(array2table(ls_positions), output_dir + "README.txt", WriteMode="append");
    writelines("AAES mic positions:", output_dir + "README.txt", WriteMode="append");
    writetable(array2table(mic_positions), output_dir + "README.txt", WriteMode="append");
    
    disp("Finished folder: "+output_dir);
end