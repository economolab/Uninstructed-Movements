tempme = me(end-3);
temprez = rez(end-3);
tempobj = obj(end-3);
tempparams = params(end-3);

%%

trix = randsample(tempobj.bp.Ntrials,1,false); % 95

plotme = normalize(tempme.data(:,trix),'range',[0 1]);

potent = squeeze(temprez.N_potent(:,trix,:));
potent = normalize(sum(potent.^2,2),'range',[0 1]);

null = squeeze(temprez.N_null(:,trix,:));
null = normalize(sum(null.^2,2),'range',[0 1]);

%
close all
dy = 1;
cols = getColors;

f = figure;
ax = gca;
hold on
plot(tempobj.time,mySmooth(plotme,7),'k','LineWidth',2)
plot(tempobj.time,mySmooth(potent,7)+dy,'Color',cols.potent,'LineWidth',2)
plot(tempobj.time,mySmooth(null,7)+dy*2,'Color',cols.null,'LineWidth',2)

xlim([tempobj.time(5), 2.5])

lickL = tempobj.bp.ev.lickL{trix};
lickR = tempobj.bp.ev.lickR{trix};
lastLick = sort([lickL lickR]);
if ~isempty(lastLick)
    lastLick = lastLick(end);
else
    lastLick = nan;
end
lastLick-2.5
xline(lastLick-2.5,'k--')
xline(0,'k--')











