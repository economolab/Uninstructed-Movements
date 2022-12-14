function plotSubspaceAlignment(ai)

dat = ai;

%% plot

cols = [0.1 0.1 0.1; 0.4 0.4 0.4];

fns = fieldnames(dat);

f=figure; hold on;
ax = gca;
div = 1;
xs = [1, 2];

% i = 1; % null afc_aw (frac of var explained of aw data using afc null space)
% temp = dat.null.afc_aw;
% h(i) = bar(xs(i), mean(temp));
% h(i).FaceColor = cols(i,:);
% h(i).EdgeColor = 'none';
% h(i).FaceAlpha = 0.5;
% scatter(xs(i)*ones(size(temp)),temp,60,'MarkerFaceColor',cols(i,:)./div, ...
%         'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',0.25, ...
%         'MarkerFaceAlpha',0.7)
% errorbar(h(i).XEndPoints,mean(temp),std(temp)./sqrt(numel(temp)),'LineStyle','none','Color','k','LineWidth',1);
% 
% i = 2; % null aw_afc (frac of var explained of afc data using aw null space)
% temp = dat.null.aw_afc;
% h(i) = bar(xs(i), mean(temp));
% h(i).FaceColor = cols(i,:);
% h(i).EdgeColor = 'none';
% h(i).FaceAlpha = 0.5;
% scatter(xs(i)*ones(size(temp)),temp,60,'MarkerFaceColor',cols(i,:)./div, ...
%         'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',0.25, ...
%         'MarkerFaceAlpha',0.7)
% errorbar(h(i).XEndPoints,mean(temp),std(temp)./sqrt(numel(temp)),'LineStyle','none','Color','k','LineWidth',1);


i = 1; % potent afc_aw 
temp = dat.potent.afc_aw;
h(i) = bar(xs(i), mean(temp));
h(i).FaceColor = cols(i,:);
h(i).EdgeColor = 'none';
h(i).FaceAlpha = 0.5;
scatter(xs(i)*ones(size(temp)),temp,60,'MarkerFaceColor',cols(i,:)./div, ...
        'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',0.25, ...
        'MarkerFaceAlpha',0.7)
errorbar(h(i).XEndPoints,mean(temp),std(temp)./sqrt(numel(temp)),'LineStyle','none','Color','k','LineWidth',1);

i = 2; % potent aw_afc 
temp = dat.potent.aw_afc;
h(i) = bar(xs(i), mean(temp));
h(i).FaceColor = cols(i,:);
h(i).EdgeColor = 'none';
h(i).FaceAlpha = 0.5;
scatter(xs(i)*ones(size(temp)),temp,60,'MarkerFaceColor',cols(i,:)./div, ...
        'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',0.25, ...
        'MarkerFaceAlpha',0.7)
errorbar(h(i).XEndPoints,mean(temp),std(temp)./sqrt(numel(temp)),'LineStyle','none','Color','k','LineWidth',1);

ax.XTick = xs;
xticklabels({'Qafc -> AW', 'Qaw -> AFC'})

ylabel('Potent Space Alignment Index')
ax.FontSize = 12;



end