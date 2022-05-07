function plotAllModes_nullPotent(latents,trials_by_type,plt)


switch plt.plot_mean
    case 0
        plotSingleTrials(latents,trials_by_type,plt);
    case 1
        plotMeanAcrossTrials(latents,trials_by_type,plt);
end







end


%% Helper Functions

function plotSingleTrials(latents,trials_by_type,plt)


plotTimeIx = latents.plotTimeIx;

trialfns = fieldnames(trials_by_type);
trialfns = trialfns(plt.trial_types);

latentfns = fieldnames(latents);
[~,mask] = patternMatchCellArray(latentfns,{'proj'},'all');
latentfns = latentfns(mask);


for latentix = 1:numel(latentfns)
    figure(210 + latentix); clf
    for i = 1:size(latents.(latentfns{latentix}),2) % number of modes
        subplot(3,2,i); hold on;
        for j = 1:numel(trialfns)
            plot(latents.time(plotTimeIx),squeeze(latents.(latentfns{latentix})(plotTimeIx,i,trials_by_type.(trialfns{j}))),'Color',plt.colors{j},'LineWidth',0.5)
        end
        hold off
        title(latents.names{i})
        ax = gca;
        ax.FontSize = 20;
    end
    sgtitle(latentfns{latentix},'Interpreter','none')
end



end % plotSingleTrials


function plotMeanAcrossTrials(latents,trials_by_type,plt)

plotTimeIx = latents.plotTimeIx;

trialfns = fieldnames(trials_by_type);
trialfns = trialfns(plt.trial_types);

latentfns = fieldnames(latents);
[~,mask] = patternMatchCellArray(latentfns,{'proj'},'all');
latentfns = latentfns(mask);


for latentix = 1:numel(latentfns)
    figure(220 + latentix); clf
    for i = 1:size(latents.(latentfns{latentix}),2) % number of modes
        subplot(3,2,i); hold on;
        for j = 1:numel(trialfns)
            plot(latents.time(plotTimeIx),mean(squeeze(latents.(latentfns{latentix})(plotTimeIx,i,trials_by_type.(trialfns{j}))),2),'Color',plt.colors{j},'LineWidth',3)
        end
        hold off
        title(latents.names{i})
        ax = gca;
        ax.FontSize = 20;
    end
    sgtitle(latentfns{latentix},'Interpreter','none')
end

end % plotMeanAcrossTrials


















