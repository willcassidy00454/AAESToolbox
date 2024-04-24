%
% This function returns the available gain before instability of the feedback
% matrix A in decibels.
%
% Negative output indicates input is unstable
%
% This may be conservative since it only c0nsiders the magnitude of the
% eigenvalues
%
% Input:
% A is a 3D tensor with dimensions (num_inputs x num_outputs x num_bins)
% where num_inputs = num_outputs
%
function gbi_dB = FindWorstCaseGBI(A)
    num_bins = size(A,3);

    all_channel_eigenvalues = zeros(round(num_bins/2),1);
    
    % num_bins / 2 used to discard symmetry
    for bin = 1:round(num_bins/2)
        all_channel_eigenvalues(bin) = max(abs(eig(A(:,:,bin))), [], "all");
    end

    % figure
    % hold on
    % semilogx(all_channel_eigenvalues);

    gbi_dB = -20 * log10(max(all_channel_eigenvalues,[],"all"));
end