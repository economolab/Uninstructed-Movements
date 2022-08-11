function plotAvgJawVelocityDuringStim(meta,obj,dfparams,params,kin,kinfeats,feats2plot,cond2plot,sav)

dfparams.plt.ms = {'.','.','.','.','.','.'};


[~,mask] = patternMatchCellArray(kin(1).featLeg,feats2plot,'any');
featix = find(mask);

times = [0 0.8]; % relative to alignEv (should be delay)
for i = 1:numel(times)
    [~,ix(i)] = min(abs(dfparams.time - times(i)));
end

% get avg feature value across time and trials for each condition
for i = 1:numel(obj)
    for k = 1:numel(featix)
        for j = 1:numel(cond2plot)
            trials = params(i).trialid{cond2plot(j)};
            temp = kinfeats{i}(ix(1):ix(2),trials,featix(k));
%             temp = normalizeInRange(temp,[0 1]);
            vel{i}{k}(j) = nanvar(nanvar(temp));
%             vel{i}{k}(j) = nanmean(nanmean(temp)); % vel{session}{feat}(cond)
        end
    end
end


%% plot
f = figure; hold on;
f.Position = [680   205   477   773];
t = tiledlayout('flow');
for k = 1:numel(featix)
    ax = nexttile; hold on;

    for i = 1:numel(obj)

        for j = 1:2:numel(cond2plot)
            s = plot(vel{i}{k}(j),vel{i}{k}(j+1),dfparams.plt.ms{j},'MarkerSize',30,'Color',dfparams.plt.color{cond2plot(j)});
        end
        title(feats2plot{k}, 'Interpreter','none');
        %         xlim([dfparams.times(1) dfparams.times(2)]);
        %         ylim([0 size(temp,2)]);
    end
    %             refline;

    xlabel('Avg No Stim')
    ylabel('Avg Stim')
    ax.FontSize = 20;
    axis(ax,'equal')

    mins = min([ax.XLim(1) ax.YLim(1)]);
    maxs = max([ax.XLim(2) ax.YLim(2)]);
    ax.XLim = [mins maxs];
    ax.YLim = [mins maxs];

    plot(ax.XLim,ax.YLim,'k--','LineWidth',2)


end



if sav
    pth = [ 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements\mc_stim\figs\'];
    fn = [ 'AvgJawVelStim_NoStim' ];
    mysavefig(f,pth,fn)
end


end