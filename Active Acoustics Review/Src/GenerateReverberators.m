
% GenerateIdentityReverberator(16, "Active Acoustics Review/Reverberators/Reverberator 1/", 48000, 24);

% hh = GenerateRandomHouseholder(4,4)

hadamard = GenerateNextHadamardIteration(GenerateNextHadamardIteration(GenerateNextHadamardIteration(GenerateNextHadamardIteration(1))));
% 
% heatmap(hadamard);
% 
SaveReverberatorMatrix(hadamard, "Active Acoustics Review/Reverberators/Reverberator 14/", 48000, 32);

% SaveIdentityReverberator(16, "Active Acoustics Review/Reverberators/Reverberator 1/", 48000, 32);

function SaveReverberatorMatrix(matrix, output_dir, sample_rate, bit_depth)
    mkdir(output_dir);

    for row = 1:size(matrix, 1)
        for col = 1:size(matrix, 2)
            audiowrite(output_dir + "X_R" + row + "_S" + col + ".wav", matrix(row,col), sample_rate, BitsPerSample=bit_depth);
        end
    end
end

function householder = GenerateRandomHouseholder(num_output_channels, num_input_channels)
    % Generate random vector with length equal to input dimension
    v = rand([1 max([num_input_channels num_output_channels])]);

    norm = 0;

    for value = v
        norm = norm + value * value;
    end

    norm = sqrt(norm);

    v = v / norm;

    householder = zeros([num_output_channels num_input_channels]);

    % Fill householder
    for row = 1:num_output_channels
        for col = 1:num_input_channels
            householder (row, col) = -2 * v(row) * v(col);

            if row == col
                householder(row, col) = householder(row, col) + 1;
            end
        end
    end
end

function output_matrix = GenerateNextHadamardIteration(input_matrix)
    num_rows = length(input_matrix);
    num_cols = num_rows;

    output_matrix = zeros(num_rows * 2, num_cols * 2);

    % Fill top left
    output_matrix(1:num_rows, 1:num_cols) = input_matrix;

    % Fill bottom left
    output_matrix(num_rows + 1:num_rows * 2, 1:num_cols) = input_matrix;

    % Fill top right
    output_matrix(1:num_cols, num_cols + 1:num_cols * 2) = input_matrix;

    % Fill bottom right
    output_matrix(num_rows + 1:num_rows * 2, num_cols + 1:num_cols * 2) = -input_matrix;
end

function SaveIdentityReverberator(num_channels, output_dir, sample_rate, bit_depth)
    mkdir(output_dir);

    for row = 1:num_channels
        for col = 1:num_channels
            output = 0;
            if row == col
                output = 1;
            end
            audiowrite(output_dir + "X_R" + row + "_S" + col + ".wav", output, sample_rate, BitsPerSample=bit_depth);
        end
    end
end