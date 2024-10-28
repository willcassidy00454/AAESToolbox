
% If the sources/receivers are part of the AAES, then simply use their IRs
% as if they were separate transducers. This duplication acts as a probing
% of the desired AAES transducers.

function GenerateAAESIRs(rir_parent_dir, num_aaes_loudspeakers, num_aaes_mics, room_num, alpha_set, loop_gain_biases_dB, output_directory)
    %% Internal parameters
    N = 1; % Number of source loudspeakers
    M = 1; % Number of receiver microphones
    K = num_aaes_loudspeakers; % Number of AAES loudspeakers
    L = num_aaes_mics; % Number of AAES microphones

    reverberator_rt_factors = [1 1.5 2 2.5 3];

    room_dims = [5.7 7.35 2.5; 8.74 17 5.5; 19.52 30.83 15];
    rir_directory = rir_parent_dir + "AAES IRs Ch["+L+"] Room"+mat2str(room_dims(room_num, :))+" AlphaSet["+alpha_set+"] SampleRate[48000]/";
    
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
    
    U = FillTransferFunctionMatrix(U, num_bins, "U", rir_directory);
    E = FillTransferFunctionMatrix(E, num_bins, "E", rir_directory);
    F = FillTransferFunctionMatrix(F, num_bins, "F", rir_directory);
    G = FillTransferFunctionMatrix(G, num_bins, "G", rir_directory);
    H = FillTransferFunctionMatrix(H, num_bins, "H", rir_directory);

    mkdir(output_directory);

    for reverberator_rt_factor = reverberator_rt_factors
        % The 16ch folders are always used for the reverberator, since
        % lower channel counts just read a subset

        reverberator_dir = rir_parent_dir + "../Pink Reverberator IRs/Decaying Noise Ch[16] Room["+room_num+"] AlphaSet["+alpha_set+"] RTFactor["+reverberator_rt_factor+"]/";
        X = FillReverberatorMatrix(X, num_bins, reverberator_dir);

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
            
            audiowrite(output_directory + "Ch["+K+"]_Room["+room_num+"]_AlphaSet["+alpha_set+"]_RTFactor["+reverberator_rt_factor+"]_LoopGain["+loop_gain_bias_dB+"].wav", output_signal, sample_rate, 'BitsPerSample', bit_depth); 
        end
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