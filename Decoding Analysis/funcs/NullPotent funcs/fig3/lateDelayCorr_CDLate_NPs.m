close all


cdix = 2;
% late delay
[~,tix(1)] = min(abs(obj(1).time - -0.4));
[~,tix(2)] = min(abs(obj(1).time - -0.02));

nullmeans = [];
potentmeans = [];

for sessix = 1:numel(rez)

    dat = rez(sessix).N_null;
    dims = size(dat);

    nullproj = reshape(dat,dims(1)*dims(2),dims(3)) * cd_null(sessix).cd_mode_orth(:,cdix);
    nullproj = reshape(nullproj,dims(1),dims(2));

    dat = rez(sessix).N_potent;
    dims = size(dat);

    potentproj = reshape(dat,dims(1)*dims(2),dims(3)) * cd_potent(sessix).cd_mode_orth(:,cdix);
    potentproj = reshape(potentproj,dims(1),dims(2));

    null = mean(nullproj(tix(1):tix(2),:),1)';
    potent = mean(potentproj(tix(1):tix(2),:),1)';

    nullmeans = [nullmeans ; null]; % late delay mean for each trial
    potentmeans = [potentmeans ; potent];

    figure; hold on;
    mdl = fitlm(null,potent);
    f = plot(mdl);
    f(1).Marker = 'o';
    f(1).MarkerSize = 5;
    f(1).MarkerFaceColor = [0.5 0.5 0.5];
    f(1).MarkerEdgeColor = 'w';
    % f(1).MarkerFaceAlpha = 0.3;

    f(2).Color = 'k';
    f(2).LineWidth = 2;

    for i = 3:4
        f(i).LineStyle = '-';
        f(i).Color = 'k';
    end
    r2(sessix) = mdl.Rsquared.Ordinary;
    title(['R^2 = ' num2str(r2(sessix))])
    xlabel('Late delay CD_Late_Null')
    ylabel('Late delay CD_Late_Potent')

    ax = gca;
    ax.Legend.Visible = 'off';
    ax.FontSize = 15;
    %     ylim([-12 ax.YLim(2)])
    hold off

end

%%

cols = [0.2 0.2 0.2];

figure;
ax = gca;
hold on;
xs = [1];
i = 1;
h(i) = bar(xs(i),mean(r2));
h(i).FaceColor = cols(i,:);
h(i).EdgeColor = 'none';
h(i).FaceAlpha = 0.5;
scatter(xs(i)*ones(size(r2)),r2,60,'MarkerFaceColor',cols(i,:), ...
    'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',.35, ...
    'MarkerFaceAlpha',0.7)
errorbar(h(i).XEndPoints,mean(r2),std(r2)./sqrt(size(r2,2)),'LineStyle','none','Color','k','LineWidth',1);

ax.XTick = [];
ylabel('R2')
ax.FontSize = 14;




















