% clear all
% close all
% tiledlayout(3,1);

%% User inputs (currently, L must equal K since the reverberator is square)

% If the sources/receivers are part of the AAES, then simply use their IRs
% as if they were separate transducers. This duplication acts as a probing
% of the desired AAES transducers.

% num_aaes_channels = 8;
% 
% room_num = 1;
% alpha_set = 1;

% Bias is added to loop gain where 0 results in the limit of instability
% loop_gain_biases_dB = [-2 -4 -6 -200];

function GenerateAAESIRs(num_aaes_channels, room_num, alpha_set, loop_gain_biases_dB, output_directory, sample_rate)
    %% Internal parameters
    N = 1; % Number of source loudspeakers
    M = 1; % Number of receiver microphones
    K = num_aaes_channels; % Number of AAES loudspeakers
    L = K; % Number of AAES microphones

    reverberator_rt_factors = [1];% 1.5 2 2.5 3];%[0 0.5 1 2 4];

    room_dims = [5.7 7.35 2.5; 8.74 17 5.5; 19.52 30.83 15];
    rir_directory = "Automated RIRs/AAES IRs Ch["+L+"x"+K+"] Room"+mat2str(room_dims(room_num, :))+" AlphaSet["+alpha_set+"]/";
    % rir_directory = "Kentish Town Lab RIRs/16-channel/";
    
    [example_ir, sample_rate] = audioread(rir_directory + "E_1_1.wav");
    passive_ir_length = size(example_ir, 1); % Use the first IR of the set to determine IR lengths
    num_bins = round(passive_ir_length * 3); % Number of frequency bins (equals length of IR to save)
    bit_depth = 24;
    
    %% Initialise matrices
    
    U = zeros(N, 1, num_bins); % U = 1xN Inputs
    E = zeros(M, N, num_bins); % E = NxM matrix of transfer functions from each room source to each observer microphone
    F = zeros(M, K, num_bins); % F = KxM matrix of transfer functions from each AAES loudspeaker to each observer microphone
    G = zeros(L, N, num_bins); % G = NxL matrix of transfer functions from each room source to each AAES microphone
    H = zeros(L, K, num_bins); % H = KxL matrix of transfer functions from each AAES loudspeaker to each AAES microphone
    X = zeros(K, L, num_bins); % X = LxK matrix of transfer functions defining the reverberator
    
    %% Fill transfer function matrices by reading IR files and performing FFTs
    
    U = FillTransferFunctionMatrix(U, num_bins, "U", rir_directory, false);
    E = FillTransferFunctionMatrix(E, num_bins, "E", rir_directory, true);
    F = FillTransferFunctionMatrix(F, num_bins, "F", rir_directory, true);
    G = FillTransferFunctionMatrix(G, num_bins, "G", rir_directory, true);
    H = FillTransferFunctionMatrix(H, num_bins, "H", rir_directory, true);

    mkdir(output_directory);
    % writelines("Worst-case gains before instabilities:", output_directory + "GBIs.txt", WriteMode="overwrite");

    for reverberator_rt_factor = reverberator_rt_factors
        % The 16ch folders are always used for the reverberator, since
        % lower channel counts just read a subset

        X = FillReverberatorMatrix(X, num_bins, "Pink Reverberator IRs/Decaying Noise Ch[16] Room["+room_num+"] AlphaSet["+alpha_set+"] RTFactor["+reverberator_rt_factor+"]/");

        gbi_dB = 0.0;

        for loop_gain_bias_dB = loop_gain_biases_dB
            disp("Computing Ch["+K+"] Room["+room_num+"] AlphaSet["+alpha_set+"] RTFactor["+reverberator_rt_factor+"] LoopGain["+loop_gain_bias_dB+"]...");

            %% Isolate feedback loop and find GBI

            feedback_loop = zeros(L, L, num_bins);
            
            for bin = 1:num_bins
                feedback_loop(:,:,bin) = X(:,:,bin) * H(:,:,bin);
            end
            
            gbi_dB = FindWorstCaseGBI(feedback_loop);
            
            %% Set loop gain to maximum before instability (minus a bias)
            
            % mu = AAES feedback loop gain
            mu = power(10, (gbi_dB + loop_gain_bias_dB) / 20);
            
            %% Compute output
            
            % V = 1xM Outputs
            V = zeros(M, 1, num_bins);
            
            % Compute output one frequency bin at a time (element-wise)
            for bin = 1:num_bins
                % V = E U
                %   + mu F (I - mu X H)^-1 X G U
                V(:, :, bin) = E(:, :, bin) * U(:, :, bin) ...
                             + mu .* F(:, :, bin) * inv(eye(K) - mu .* X(:, :, bin) * H(:, :, bin)) * X(:, :, bin) * G(:, :, bin) * U(:, :, bin);
            end
            
            % Convert receiver transfer function back to the time domain
            output_signal = ifft(squeeze(V(1, 1, :)));
            
            %% Save output
            
            audiowrite(output_directory + "ReverberatorRTFactor["+reverberator_rt_factor+"]_LoopGain["+loop_gain_bias_dB+"].wav", output_signal, sample_rate, 'BitsPerSample', bit_depth); 
        end

        % Write GBI to .txt
        % writelines("Ch["+K+"] Room["+room_num+"] AlphaSet["+alpha_set+"] ReverberatorRTFactor["+reverberator_rt_factor+"]: " + gbi_dB + " dB", output_directory + "GBIs.txt", WriteMode="append");
    end

    disp("Finished folder: " + output_directory);
end

%% Functions

function matrix_to_fill = FillTransferFunctionMatrix(matrix_to_fill, desired_ir_length, filename_base_id, ir_directory)
    num_rows = size(matrix_to_fill,1);
    num_cols = size(matrix_to_fill,2);

    % Load each IR, zero pad, take FFT and insert into transfer function matrix
    for row = 1:num_rows
        for col = 1:num_cols
            padded_ir = zeros(1, desired_ir_length);
    
            [raw_ir, ~] = audioread(ir_directory + filename_base_id + "_" + col + "_" + row + ".wav");
        
            nonzero_length = min(length(raw_ir), desired_ir_length); % Iterate up to the end of the audio, truncating if too long

            for sample_pos = 1:nonzero_length
                padded_ir(sample_pos) = raw_ir(sample_pos);
            end
        
            matrix_to_fill(row, col, :) = fft(padded_ir);
        end
    end
end

% Assumes the reverberator matrix is diagonal
function matrix_to_fill = FillReverberatorMatrix(matrix_to_fill, desired_ir_length, ir_directory)
    num_rows = size(matrix_to_fill,1);
    num_cols = size(matrix_to_fill,2);

    % Load each IR, zero pad, take FFT and insert into transfer function matrix
    for row = 1:num_rows
        for col = 1:num_cols
            padded_ir = zeros(1, desired_ir_length);
    
            % Only read the diagonal IRs from file (since all others are zero)
            if (row == col)
                [raw_ir, ~] = audioread(ir_directory + "X_" + col + "_" + row + ".wav");
            
                nonzero_length = min(length(raw_ir), desired_ir_length); % Iterate up to the end of the audio, truncating if too long
    
                for sample_pos = 1:nonzero_length
                    padded_ir(sample_pos) = raw_ir(sample_pos);
                end
            end
        
            matrix_to_fill(row, col, :) = fft(padded_ir);
        end
    end
end