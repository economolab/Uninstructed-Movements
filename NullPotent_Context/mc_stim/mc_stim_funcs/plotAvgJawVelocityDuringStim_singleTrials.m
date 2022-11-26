function plotAvgJawVelocityDuringStim_singleTrials(meta,obj,dfparams,params,kin,kinfeats,feats2plot,cond2plot,sav)

dfparams.plt.ms = {'.','.','.','.','.','.'};


[~,mask] = patternMatchCellArray(kin(1).featLeg,feats2plot,'any');
featix = find(mask);

times = [0 0.8]; % relative to alignEv (should be delay)
for i = 1:numel(times)
    [~,ix(i)] = min(abs(dfparams.time - times(i)));
end

% get avg feature value across time and trials for each condition
toplot = cell(numel(feats2plot),1);
for k = 1:numel(featix)
    for j = 1:numel(cond2plot)
        fracmov{j} = [];
        toplot{k}{j} = [];
        for i = 1:numel(obj)
            trials = params(i).trialid{cond2plot(j)};
            temp = kinfeats{i}(ix(1):ix(2),trials,featix(k));

            if strcmpi(feats2plot{k},'motion_energy')
                % how much time spent moving
                mask = temp > 15;
                fracmov{j} = [fracmov{j} ; (sum(mask,1) ./ size(mask,1))'];
            end

            temp = nanmean(temp,1);
            toplot{k}{j} = [toplot{k}{j} ; temp'];
        end
    end
end




%% plot

f = figure; hold on;
t = tiledlayout('flow');
for k = 1:numel(featix)
    for j = 1:numel(cond2plot)
        ax(j) = nexttile; hold on;
        h{j} = histfit(toplot{k}{j},10,'kernel');
        h{j}(1).EdgeColor = 'none';
        h{j}(1).FaceColor = dfparams.plt.color{cond2plot(j)};
        h{j}(2).Color = 'k';
        %         histogram(ax(j),toplot{k}{j},50,'FaceColor',dfparams.plt.color{cond2plot(j)},'EdgeColor','none');
        ax(j).FontSize = 12;
        ax(j).XLim = [0 50];
    end
end
sgtitle('Motion Energy')
for i = 1:numel(h)
    ys(:,i) = h{1}(1).YData;
end
mn = 0;
mx = max(ys(:));
for i = 1:numel(ax)
    ax(i).YLim = [mn mx];
end

f = figure; hold on;
t = tiledlayout('flow');
for j = 1:numel(cond2plot)
    ax(j) = nexttile; hold on;
    h{j} = histfit(fracmov{j},10,'kernel');
    h{j}(1).EdgeColor = 'none';
    h{j}(1).FaceColor = dfparams.plt.color{cond2plot(j)};
    h{j}(2).Color = 'k';
    %     histogram(ax1,fracmov{j},50,'FaceColor',dfparams.plt.color{cond2plot(j)},'EdgeColor','none');
    ax(j).FontSize = 12;
    ax(j).XLim = [0 1];
end
sgtitle('Frac. time spent moving during delay')
for i = 1:numel(h)
    ys(:,i) = h{1}(1).YData;
end
mn = 0;
mx = max(ys(:));
for i = 1:numel(ax)
    ax(i).YLim = [mn mx];
end


if sav
    pth = [ 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements\mc_stim\figs\'];
    fn = [ 'AvgJawVelStim_NoStim' ];
    mysavefig(f,pth,fn)
end








end