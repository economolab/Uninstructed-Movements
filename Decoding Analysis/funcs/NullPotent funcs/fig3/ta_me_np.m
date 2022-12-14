for sessix = 1:numel(meta)

    tempme = me(sessix);
    temprez = rez(sessix);
    tempobj = obj(sessix);
    tempparams = params(sessix);

    % trix = randsample(tempobj.bp.Ntrials,1,false); % 95
    trix = tempparams.trialid{1};
    temp = tempme.data(:,trix);
    plotme(:,sessix) = mean(temp,2);

    temp = squeeze(temprez.N_potent(:,trix,:));
    temp = squeeze(sum(temp.^2,2));
    potent(:,sessix) = normalize(mean(temp,2),'range',[0 1]);

end


%%
cols = getColors;
alph = 0.2;

f = figure;
ax = gca;
hold on

yyaxis(ax,'right')
means = nanmean(plotme,2);
errs = nanstd(plotme,[],2) ./ sqrt(numel(meta));
shadedErrorBar(tempobj.time,means,errs,{'LineWidth',2,'Color','k','LineStyle','-'},alph,ax);
xlabel('Time (s) from go cue')
ylabel('Motion Energy')
ax.YColor = [0 0 0];

yyaxis(ax,'left')
means = mean(potent,2);
errs = std(potent,[],2) ./ sqrt(numel(meta));
shadedErrorBar(tempobj.time,means,errs,{'LineWidth',2,'Color',cols.potent,'LineStyle','-'},alph,ax);
ylabel('Potent - Magnitude (a.u.)')
ax.YColor = [0.1500    0.1500    0.1500];


xlim([tempobj.time(5), 2])
xline(0,'k--')
sample = mode(tempobj.bp.ev.sample) - 2.5;
delay = mode(tempobj.bp.ev.delay) - 2.5;
xline(sample,'k--')
xline(delay,'k--')







