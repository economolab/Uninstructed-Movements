function plotVarianceExplained_NP_overTime(rez)

dt = 20;

for sessix = 1:numel(rez)
    full = rez(sessix).N.full_reshape;
    tic
    ct = 1;
    for t = 1:dt:size(full,1)
        fullcov = cov(full(t:t+dt,:));
        e = sum(eig(fullcov));
        ve.null{sessix}(ct) = var_proj(rez(sessix).Qnull,fullcov,e);
        ve.potent{sessix}(ct) = var_proj(rez(sessix).Qpotent,fullcov,e);
        ct = ct + 1;
    end
    toc
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% curate data
ve.null_total = zeros(numel(rez),1);
ve.potent_total = zeros(numel(rez),1);
ve.null_prep = zeros(numel(rez),1);
ve.null_move = zeros(numel(rez),1);
ve.potent_move = zeros(numel(rez),1);
ve.potent_prep = zeros(numel(rez),1);
for sessix = 1:numel(rez)
    ve.null_total(sessix) = rez(sessix).ve.null_total;
    ve.potent_total(sessix) = rez(sessix).ve.potent_total;
    ve.null_prep(sessix) = rez(sessix).ve.null_prep;
    ve.null_move(sessix) = rez(sessix).ve.null_move;
    ve.potent_move(sessix) = rez(sessix).ve.potent_move;
    ve.potent_prep(sessix) = rez(sessix).ve.potent_prep;
end


%% plot

cols(1,:) = [50 50 50];
cols(2,:) = [100 100 100];
cols(3,:) = [128, 0, 0];
cols(4,:) = [4, 92, 1];
cols(5,:) = cols(4,:);
cols(6,:) = cols(3,:);
cols = cols ./ 255;

% cols = linspecer(10);

fns = fieldnames(ve);

f=figure; hold on;
ax = gca;
div = 1;
xs = [1 2 4 5 7 8];
for i = 1:numel(fns)
    temp = ve.(fns{i});
    h(i) = bar(xs(i),mean(temp)); % mean across dims and sessions
%     if i<3
%         cix = i;
%     else
%         cix = i + 2;
%     end
    cix = i;
    h(i).FaceColor = cols(cix,:);
    h(i).EdgeColor = 'none';
    h(i).FaceAlpha = 0.5;
    scatter(xs(i)*ones(size(temp)),temp,60,'MarkerFaceColor',cols(cix,:)./div, ...
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