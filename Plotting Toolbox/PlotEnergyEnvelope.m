%
% This function plots the energy over time for an IR using forward
% integration. This uses rectangular windows with no overlap.
%
% Optional Parameters:
%
% integration_length_samples - the number of samples each window spans (defaults to 500)
% normalise - whether the energy should be normalised to 0 dB or not (defaults to false)
% line_style - defaults to solid ("-")

function PlotEnergyEnvelope(ir, sample_rate, integration_length_samples, normalise, line_style)
    if ~exist('integration_length_samples','var')
        integration_length_samples = 500;
    end    

    if ~exist('normalise','var')
        normalise = false;
    end
    
    if ~exist('line_style','var')
        line_style = "-";
    end    

    % Square IR
    ir = power(ir, 2);
    
    % For each integration window along the length of the IR, calculate
    % integration and store in vector of energy values
    
    num_integration_steps = ceil(length(ir) / integration_length_samples);
    integrated_values = zeros(1, num_integration_steps);
    
    current_integration_step = 1;
    
    for integration_start_sample = 1:integration_length_samples:length(ir)
        % Extract current integration samples from IR
        current_integration_samples = ir(integration_start_sample:min(integration_start_sample + integration_length_samples,length(ir)));
    
        % Calculate integration
        current_integration_value = trapz(current_integration_samples);
    
        % Store in integration vector
        integrated_values(current_integration_step) = current_integration_value;
        current_integration_step = current_integration_step + 1;
    end
    
    if normalise
        integrated_values = integrated_values / max(abs(integrated_values));
    end

    integrated_values_dB = 10*log10(integrated_values);

    time_values_samples = 0:integration_length_samples:num_integration_steps * integration_length_samples - integration_length_samples;
    time_values_seconds = time_values_samples / sample_rate;
    
    % Plot energy values over time
    plot(time_values_seconds,integrated_values_dB,'LineStyle',line_style);
    title("IR Energy Envelope","FontSize",15);
    xlabel("Time (s)");
    ylabel("Energy (dB)");
    set(gcf,'position',[400,2000,700,500]);
    grid on
end