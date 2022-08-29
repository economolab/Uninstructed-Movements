function plotPerformanceAllMice(meta,obj,rez,dfparams,params,cond2use,connect)

perf = rez(1).perf(cond2use);
for i = 2:numel(rez)
    perf = cat(1,perf,rez(i).perf(cond2use));
end
perf = perf .* 100;

if numel(meta) == 1
    f = figure;
    f.Position = [316          73        1232         905];
    ax = axes(f); hold on;
    for i = 1:numel(perf)
        b(i) = bar(categorical(dfparams.cond(cond2use(i))),perf(i));
        b(i).FaceColor = dfparams.plt.color{cond2use(i)};
        b(i).EdgeColor = 'none';
    end
    ylabel('Performance (%)')
    ylim([0,100])
    title('Performance on all sessions, all mice')
    ax.FontSize = 13;
    hold off;
    return
end


f = figure;
f.Position = [316          73        1232         905];
ax = axes(f);
violincols = reshape(cell2mat(dfparams.plt.color(cond2use)),3,numel(cond2use))';
vs = violinplot(perf,dfparams.cond(cond2use),...
    'EdgeColor',[1 1 1], 'ViolinAlpha',{0.2,1}, 'ViolinColor', violincols,'Width',0.2);
ylabel('Performance (%)')
ylim([0,100])
title('Performance on all sessions, all mice')
ax.FontSize = 13;


if connect
    % for each pair of conditions, plot connecting lines
    hold on;
    for i = 1:2:(numel(cond2use)-1)
        conds = [i i+1];
        xs = [vs(conds(1)).ScatterPlot.XData ; vs(conds(2)).ScatterPlot.XData];
        ys = [vs(conds(1)).ScatterPlot.YData ; vs(conds(2)).ScatterPlot.YData];
        patchline(xs, ys,'EdgeColor','k','EdgeAlpha',0.3)

%         plot(xs,ys,'k-')
    end
end

% t-test b/w groups
ct = 1;
for i = 1:2:(numel(cond2use)-1)
    conds = [i i+1];
    xs = [vs(conds(1)).ScatterPlot.XData ; vs(conds(2)).ScatterPlot.XData];
    ys = [vs(conds(1)).ScatterPlot.YData ; vs(conds(2)).ScatterPlot.YData];
    %     [h(ct),p(ct)] = ttest2(ys(1,:),ys(2,:));
    [p(ct),h(ct)] = ranksum(ys(1,:),ys(2,:)); % if unsure both datasets have equal variance, use mann u whtiney instead of t test

    if h(ct)
        lw = 0.4;
        x = [i i+1];
        y = [max(ys(1,:)) max(ys(2,:))] + 3;
        plot([x(1) x(1)],[y(1) max(y)+3],'k-','LineWidth',lw)
        plot([x(2) x(2)],[y(2) max(y)+3],'k-','LineWidth',lw)
        plot([x(1) x(2)],[max(y)+3 max(y)+3],'k-','LineWidth',lw)

        plot(mean(x),max(y)+4,'k*')


    end

    ct = ct + 1;


end



end






