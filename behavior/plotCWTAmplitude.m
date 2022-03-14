function plotCWTAmplitude(time,amplitudes,freq,labels,trix)


f = figure(102);
sgtitle(['Trial ' num2str(trix)])
for i = 1:numel(amplitudes)
    subplot(2,3,i)
    imagesc('XData',time,'YData',freq,'CData',amplitudes{i})
    
    title(labels{i})
    xlabel('Time')
    ylabel('Frequency')
    ax = f.CurrentAxes;
%     ax.YTickLabel = freq;
    ax.FontSize = 20;
    cbar = colorbar;
    cbar.Label.String = 'CWT Amplitude';
    
end

end
