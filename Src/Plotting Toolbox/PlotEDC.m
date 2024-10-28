% 
% This function plots an Energy Decay Curve (aka Schroeder curve) of an IR
% using backwards integration.
%
% Optional Parameters:
%
% line_style - defaults to solid ("-")
% x_max - Maximum time value to plot in seconds
% y_min - Minimum energy value to plot in dB (defaults to -70 dB)
% integration_limit_sec - integration limit in seconds (defaults to the length of the IR)

function PlotEDC(ir, sample_rate, line_style, x_max, y_min, integration_limit_sec)
    integration_limit_samples = 0;

    if ~exist('line_style','var')
        line_style = "-";
    end

    if ~exist('y_min','var')
        y_min = -50;
    end

    if exist('integration_limit_sec','var')
        integration_limit_samples = integration_limit_sec * sample_rate;
    else
        integration_limit_samples = size(ir,1);
    end

    % Calculate Schroeder decay
    integrated_values_dB(integration_limit_samples:-1:1)=10*log10(cumsum(ir(integration_limit_samples:-1:1).^2)/sum(ir(1:integration_limit_samples).^2));

    time_values_samples = 1:length(integrated_values_dB);
    time_values_seconds = time_values_samples / sample_rate;
    
    % Plot energy values over time
    plot(time_values_seconds,integrated_values_dB,'LineStyle',line_style);
    title("Energy Decay Curve","FontSize",12);
    xlabel("Time (s)");
    ylabel("Energy (dB)");
    set(gcf,'position',[400,2000,350,250]);
    ylim([y_min 0]);
    grid on

    if exist('x_max','var')
        xlim([0 x_max]);
    end
end