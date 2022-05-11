function plotAllModes_nullPotent(obj,params,latents,trials_by_type,plt)


switch plt.plot_mean
    case 0
        plotSingleTrials(obj,params,latents,trials_by_type,plt);
    case 1
        plotMeanAcrossTrials(obj,params,latents,trials_by_type,plt);
end



end


%% Helper Functions

function plotSingleTrials(obj,params,latents,trials_by_type,plt)

align = mode(obj.bp.ev.(params.alignEvent));
sample = mode(obj.bp.ev.sample) - align;
delay = mode(obj.bp.ev.delay) - align;


plotTimeIx = latents.plotTimeIx;

trialfns = fieldnames(trials_by_type);
trialfns = trialfns(plt.trial_types);

latentfns = fieldnames(latents);
[~,mask] = patternMatchCellArray(latentfns,{'proj'},'all');
latentfns = latentfns(mask);


for latentix = 1:numel(latentfns)
    f = figure(210 + latentix); clf
    f.Position = [-1120         281         355         514];
    for i = 1:size(latents.(latentfns{latentix}),2) % number of modes
        ax = subplot(2,1,i); hold on;
        for j = 1:numel(trialfns)
            plot(latents.time(plotTimeIx),squeeze(latents.(latentfns{latentix})(plotTimeIx,i,trials_by_type.(trialfns{j}))),'Color',plt.colors{j},'LineWidth',0.5)
        end
        
        xline(sample,'k--','LineWidth',2)
        xline(delay,'k--','LineWidth',2)
        xline(0,'k--','LineWidth',2)
        
        
        if i~=2
            ax.XTick = [];
        else
            ax.XLabel.String = 'Time (s) from go cue';
        end
                
        hold off
        title(latents.names{i})
        ax = gca;
        ax.FontSize = 20;
        
        xlim([obj.time(latents.plotTimeIx(1)), obj.time(latents.plotTimeIx(end))]);
    end
    sgtitle(latentfns{latentix},'Interpreter','none')
end



end % plotSingleTrials


function plotMeanAcrossTrials(obj,params,latents,trials_by_type,plt)

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


















