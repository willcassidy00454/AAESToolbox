%
% This function plots a spectrogram of the IR
%
% Optional Parameters:
%
% x_limit - maximum time to plot in seconds (defaults to length of IR)

function PlotSpectrogram(ir,sample_rate,x_limit)
%     M = floor(length(ir)/32);
%     g =  hamming(M);
%     L = floor(M/2);
%     Ndft = max(128,2^nextpow2(M));
%     
%     [sx,fx,tx] = spectrogram(ir);
%     
%     [st,ft,tt] = stft(ir,Window=g,OverlapLength=L, ...
%         FFTLength=Ndft,FrequencyRange="onesided");
%     
%     dff = max(max(sx-st));
% 
%     figure
%     waterplot(st,ft/pi,tt/sample_rate);

    if ~exist('x_limit','var')
        x_limit = size(ir,1) / sample_rate;
    end

    nexttile
    %[S,F,T] = 
    pspectrum(ir,sample_rate,"spectrogram","FrequencyLimits",[20 20000],"OverlapPercent",80,"FrequencyResolution",20,"MinThreshold",-105);
    % time_values = T;
    % pcolor(time_values,F,log10(abs(S)));
    % shading flat
    % colorbar
    set(gca,'Yscale','log');
    set(gcf,'position',[600, 400, 500, 400]);
    xlim([0, x_limit]);
    ylim([0.02, 20]);
    xlabel("Time (s)");
    ylabel("Frequency / Hz");
    yticks([0.02, 0.2, 2, 20]);
    yticklabels(["20", "200", "2k", "20k"]);
end