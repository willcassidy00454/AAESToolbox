
% This script is a standalone for calculating GBis

% Load transfer functions

num_channels = 16;
rt_factor = 4;

rir_directory = "Kentish Town Lab RIRs/"+num_channels+"-channel/";

[ir_example, fs] = audioread(rir_directory + "H_1_1.wav");

num_bins = size(ir_example, 1);

H = zeros(num_channels, num_channels, num_bins);
X = zeros(num_channels, num_channels, num_bins);

H = FillTransferFunctionMatrix(H, num_bins, "H", rir_directory);
X = FillReverberatorMatrix(X, num_bins, "Reverberator IRs/Decaying Noise Ch[16] Room[2] AlphaSet[2] RTFactor["+rt_factor+"]/");

% Calculate feedback loop per bin

feedback_loop = zeros(num_channels, num_channels, num_bins);

for bin = 1:num_bins
    feedback_loop(:,:,bin) = X(:,:,bin) * H(:,:,bin);
end

% Calculate GBI
gbi_dB = FindWorstCaseGBI(feedback_loop);

disp(gbi_dB);

function matrix_to_fill = FillTransferFunctionMatrix(matrix_to_fill, desired_ir_length, filename_base_id, ir_directory)
    num_rows = size(matrix_to_fill,1);
    num_cols = size(matrix_to_fill,2);

    % Load each IR, zero pad, take FFT and insert into transfer function matrix
    for row = 1:num_rows
        for col = 1:num_cols
            padded_ir = zeros(1, desired_ir_length);
    
            [raw_ir, ~] = audioread(ir_directory + filename_base_id + "_" + row + "_" + col + ".wav");
        
            nonzero_length = min(length(raw_ir), desired_ir_length); % Iterate up to the end of the audio, truncating if too long

            for sample_pos = 1:nonzero_length
                padded_ir(sample_pos) = raw_ir(sample_pos);
            end
        
            matrix_to_fill(row, col, :) = fft(padded_ir);
        end
    end
end

% This function assumes the reverberator matrix is diagonal
function matrix_to_fill = FillReverberatorMatrix(matrix_to_fill, desired_ir_length, ir_directory)
    num_rows = size(matrix_to_fill,1);
    num_cols = size(matrix_to_fill,2);

    % Load each IR, zero pad, take FFT and insert into transfer function matrix
    for row = 1:num_rows
        for col = 1:num_cols
            padded_ir = zeros(1, desired_ir_length);
    
            % Only read the diagonal IRs from file (since all others are zero)
            if (row == col)
                [raw_ir, ~] = audioread(ir_directory + "X_" + row + "_" + col + ".wav");
            
                nonzero_length = min(length(raw_ir), desired_ir_length); % Iterate up to the end of the audio, truncating if too long
    
                for sample_pos = 1:nonzero_length
                    padded_ir(sample_pos) = raw_ir(sample_pos);
                end
            end
        
            matrix_to_fill(row, col, :) = fft(padded_ir);
        end
    end
end