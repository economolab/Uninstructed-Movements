function plotCondAvgMEandCD(cond2plot, sessix, params, obj, meta, kin, regr, alph, cols, MEix,times)
sm = 21;

for c = 1:length(cond2plot)
    cond = cond2plot(c);
    condtrix = params(sessix).trialid{cond};
    condME = squeeze(kin(sessix).dat(:,condtrix,MEix));
    presampME = mean(condME(times.startix:times.stopix,:),1,'omitnan');
    presampME = mean(presampME);
    condME = condME-presampME;
    nTrials = size(condME,2);

    ax1 = subplot(2,1,1);
    toplot = mean(condME,2,'omitnan');
    err = 1.96*(std(condME,0,2,'omitnan')./sqrt(nTrials));
    ax = gca;
    shadedErrorBar(obj(sessix).time,mySmooth(toplot,sm),err,{'Color',cols{c},'LineWidth',2},alph,ax); hold on;

    ax2 = subplot(2,1,2);
    condRamp = regr(sessix).singleProj(:,condtrix);
    toplot = mean(condRamp,2,'omitnan');
    err = 1.96*(std(condRamp,0,2,'omitnan')./sqrt(nTrials));
    ax = gca;
    shadedErrorBar(obj(sessix).time,toplot,err,{'Color',cols{c},'LineWidth',2},alph,ax); hold on;
end
xlabel(ax1,'Time from go cue (s)')
ylabel(ax1,'Motion energy (a.u.)')
xline(ax1,0,'k--','LineWidth',1)
xline(ax1,-0.9,'k--','LineWidth',1)
xline(ax1,-2.2,'k--','LineWidth',1)
xlim(ax1,[-2.3 0])

xlabel(ax2,'Time from go cue (s)')
ylabel(ax2,'CDRamping (a.u.)')
xline(ax2,0,'k--','LineWidth',1)
xline(ax2,-0.9,'k--','LineWidth',1)
xline(ax2,-2.2,'k--','LineWidth',1)
xlim(ax2,[-2.3 0])
legend({'Right','Left'},'Location','best')

sgtitle([meta(sessix).anm ' ' meta(sessix).date '; Sesh ' num2str(sessix)])
end