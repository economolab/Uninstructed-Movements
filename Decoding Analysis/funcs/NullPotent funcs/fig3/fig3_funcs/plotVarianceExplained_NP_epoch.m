function plotVarianceExplained_NP_epoch(rez)


% curate data
ve.null_prep = zeros(numel(rez),1);
ve.null_move = zeros(numel(rez),1);
ve.potent_move = zeros(numel(rez),1);
ve.potent_prep = zeros(numel(rez),1);
for sessix = 1:numel(rez)
    ve.null_prep(sessix) = rez(sessix).ve.norm.null_prep;
    ve.null_move(sessix) = rez(sessix).ve.norm.null_move;
    ve.potent_move(sessix) = rez(sessix).ve.norm.potent_move;
    ve.potent_prep(sessix) = rez(sessix).ve.norm.potent_prep;
end


%% plot

cols{1} = [0.7 0.7 0.7];
cols{2} = [0.2 0.2 0.2];
cols{3} = cols{2};
cols{4} = cols{1};

fns = fieldnames(ve);

f=figure; hold on;
ax = gca;
div = 1;
xs = [1 2 5 6];
for i = 1:numel(fns)
    temp = ve.(fns{i});
    h(i) = bar(xs(i),mean(temp)); % mean across dims and sessions
%     if i<3
%         cix = i;
%     else
%         cix = i + 2;
%     end
    cix = i;
    h(i).FaceColor = cols{i};
    h(i).EdgeColor = 'none';
    h(i).FaceAlpha = 0.9;
    scatter(xs(i)*ones(size(temp)),temp,60,'MarkerFaceColor',cols{i}./div, ...
        'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',0.25, ...
        'MarkerFaceAlpha',0.7)
    errorbar(h(i).XEndPoints,mean(temp),std(temp)./sqrt(numel(temp)),'LineStyle','none','Color','k','LineWidth',1);
end
% ylim([-0.001 ax.YLim(2)])
ax.XTick = xs;
xlabels  = strrep(fns,'_','-');
xticklabels(xlabels);

ylabel('Fraction of VE')
ax.FontSize = 12;



end