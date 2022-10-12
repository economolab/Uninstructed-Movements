function plotSingleTrial_CDdelay_ME(obj,params,rez,me,meta,cond2use)



cols = getColors();
clrs{1} = cols.rhit;
clrs{2} = cols.lhit;

for sessix = 1:numel(rez)

    
    trix = cell2mat({params(sessix).trialid{cond2use}}');

    cddelayix = find(ismember( rez(sessix).cd_labels,'late'));
    plotcd = rez(sessix).cd_proj_trialdat(:,trix,cddelayix);

    plotme = me(sessix).data(:,trix);

    sample = mode(obj(sessix).bp.ev.sample) - 2.5;
    delay = mode(obj(sessix).bp.ev.delay) - 2.5;

%     f = figure;
%     imagesc(obj(sessix).time,1:numel(trix),plotcd')
% 
%     f = figure;
%     imagesc(obj(sessix).time,1:numel(trix),plotme')

    f = figure;
    subplot(2,1,1); hold on;
    plot(obj(sessix).time,-mySmooth(mean(plotcd(:,1:numel(params(sessix).trialid{2})),2),21),'Color',clrs{1},'LineWidth',2);
    plot(obj(sessix).time,-mySmooth(mean(plotcd(:,1+numel(params(sessix).trialid{2}):end),2),21),'Color',clrs{2},'LineWidth',2);
    xline(0,'k:','LineWidth',2)
    xline(sample,'k:','LineWidth',2)
    xline(delay,'k:','LineWidth',2)
    xlim([-1.7 0.05])
    xlabel('Time from go cue (s)')
    ylabel('Activity (a.u.)')
    
    subplot(2,1,2); hold on;
    plot(obj(sessix).time,mySmooth(mean(plotme(:,1:numel(params(sessix).trialid{2})),2),21),'Color',clrs{1},'LineWidth',2);
    plot(obj(sessix).time,mySmooth(mean(plotme(:,1+numel(params(sessix).trialid{2}):end),2),21),'Color',clrs{2},'LineWidth',2);
    xline(0,'k:','LineWidth',2)
    xline(sample,'k:','LineWidth',2)
    xline(delay,'k:','LineWidth',2)
    xlim([-1.7 0.05])
    xlabel('Time from go cue (s)')
    ylabel('Motion Energy')

    sgtitle([meta(sessix).anm '_' meta(sessix).date],'Interpreter','none')
end



end 