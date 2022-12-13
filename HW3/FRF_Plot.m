function FRF_Plot(name, H, frq_axis, maxFreq, fig, titleString)
    figure(fig)
    if exist('titleString','var')
        title(titleString)
    end

    ax1 = subplot(2, 1, 1);
    ax2 = subplot(2, 1, 2);
%     subfigure 211
    hold on
    axes(ax1)

    plot(frq_axis, db(abs(H)), LineWidth=1.2);
    ax = gca; 
    ax.XMinorGrid = 'on'; ax.YMinorGrid = 'on';
    xlabel('freq [Hz]'); ylabel(['|',name,'| [dB]'])
    xlim([1, maxFreq])
    if exist('titleString','var')
        title(strcat(titleString)," Amplitude")
    end
%     subfigure 212
    hold on
    axes(ax2)

    plot(frq_axis, angle(H), LineWidth=1.2);
    ax = gca; 
    ax.XMinorGrid = 'on'; ax.YMinorGrid = 'on';
    xlabel('freq [Hz]'); ylabel(['\angle',name,' [rad]'])
    xlim([1, maxFreq])
    if exist('titleString','var')
        title(strcat(titleString)," Phase")
    end
end