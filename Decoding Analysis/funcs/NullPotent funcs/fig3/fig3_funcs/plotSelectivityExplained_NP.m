function plotSelectivityExplained_NP(allrez_potent,allrez_null,rez_potent,rez_null,sav,titlestring,plotmiss)


cols = getColors;
% cols(1,:) = [0.5 0.5 0.5];
% cols(2,:) = [0 0 0];

sample = mode(rez_potent(1).ev.sample) - mode(rez_potent(1).align);
delay  = mode(rez_potent(1).ev.delay) - mode(rez_potent(1).align);

lw = 2;
alph = 0.2;
f = figure; ax = axes(f); hold on;

sel = allrez_null.selexp;
% sel = allrez_null.selectivity_squared;

tempdat = squeeze(sel(:,4,:));
tempmean = nanmean(tempdat,2);
temperr = nanstd(tempdat,[],2) ./ sqrt(numel(rez_potent));
shadedErrorBar(rez_potent(1).time,tempmean,temperr,...
    {'Color',cols.null,'LineWidth',lw},alph,ax);



sel = allrez_potent.selexp;
% sel = allrez_potent.selectivity_squared;

tempdat = squeeze(sel(:,4,:));
tempmean = nanmean(tempdat,2);
temperr = nanstd(tempdat,[],2) ./ sqrt(numel(rez_potent));
shadedErrorBar(rez_potent(1).time,tempmean,temperr,...
    {'Color',cols.potent,'LineWidth',lw},alph,ax);


xline(sample,'k--','LineWidth',2)
xline(delay,'k--','LineWidth',2)
xline(0,'k--','LineWidth',2)

xlabel('Time (s) from go cue')
ylabel('Selectivity Ratio')
title(titlestring)
% legend('Total selectivity','early','late','go','early + late + go')
xlim([rez_potent(1).time(5),2])
% ylim([0 1])
ax = gca;
ax.FontSize = 14;


end