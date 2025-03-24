%
% This function plots a spectrogram of the IR
%
% Optional Parameters:
%
% x_limit - maximum time to plot in seconds (defaults to length of IR)

function PlotSpectrogram(ir,sample_rate,x_limit,should_normalise)
    if ~exist('x_limit','var')
        x_limit = size(ir,1) / sample_rate;
    end

    if ~exist('should_normalise','var')
        should_normalise = false;
    end

    nexttile

    if (should_normalise)
        [p, ~, ~] = pspectrum(ir,sample_rate,"spectrogram","FrequencyLimits",[20 20000],"OverlapPercent",90,"FrequencyResolution", 100);
        max_power_dB = 10 * log10(max(p, [], "all"));
        gain = power(10, abs(max_power_dB) / 20);
        ir = ir * gain;
    end

    min_threshold = -60;

    pspectrum(ir,sample_rate,"spectrogram","FrequencyLimits",[20 20000],"OverlapPercent",90,"FrequencyResolution",100,"MinThreshold",min_threshold);
    % set(gcf,'position',[600, 400, 500, 400]);
    xlabel("Time (s)");
    ylabel("Frequency / Hz");
    clim([min_threshold, 0]);

    set(gca,'Yscale','log');
    ylim([0.02, 20]);
    yticks([0.02, 0.2, 2, 20]);
    yticklabels(["20", "200", "2k", "20k"]);

    xlim([0, x_limit]);

    set(gca,'fontsize', 15);
    title("Title","FontSize",16,"FontWeight","normal");
end