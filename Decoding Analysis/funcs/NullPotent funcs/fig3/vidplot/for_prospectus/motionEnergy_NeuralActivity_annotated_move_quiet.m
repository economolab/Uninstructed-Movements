% can run this script after running st_elsayed/pca/kaufman.m

% set session data to use
anm = 'JEB15';
date = '2022-07-27';

anms = {meta(:).anm};
dates = {meta(:).date};

anmix = ismember(anms,anm);
dateix = ismember(dates,date);
sessix = find(all([anmix;dateix],1));

dat.rez = rez(sessix);
dat.obj = obj(sessix);
dat.params = params(sessix);
dat.meta = meta(sessix);
dat.me = me(sessix);


%%

rng(5)

close all

k = 20;

trix = randsample(dat.obj.bp.Ntrials,k,false);
tempme = dat.me;


trialOffset = 0;
f = figure; hold on;
dy = 50;
for i = 1:numel(trix)
    plotme = mySmooth(tempme.data(:,trix(i)),21);

    ix = tempme.data(:,trix(i))>(tempme.moveThresh+5);
    z = plotme;
    z(~ix) = nan;
    patchline(dat.obj.time,trialOffset + plotme,'EdgeColor','k','EdgeAlpha',0.6,'LineWidth',2);
    plot(dat.obj.time,trialOffset + z,'r','LineWidth',1)
    trialOffset = trialOffset + dy;
end
xlabel('Time (s) from go cue')
ylabel('Trials')
title('Motion Energy')
xlim([dat.obj.time(15),dat.obj.time(end)]);
ax = gca;
ax.YTick = [];
ax.FontSize = 14;

align = 2.5;
sample = mode(dat.obj.bp.ev.sample) - align;
delay = mode(dat.obj.bp.ev.delay) - align;
xline(sample,'k--','LineWidth',1)
xline(delay,'k--','LineWidth',1)
xline(0,'k--','LineWidth',1)

hold off


%%


rng(5)

close all

k = 1;
trix = randsample(dat.obj.bp.Ntrials,k,false);
tempme = dat.me;


trialOffset = 0;
f = figure; hold on;
for i = 1:numel(trix)
    dy = 50;

    ix = tempme.data(:,trix(i))>(tempme.moveThresh+5);

    plotdat = dat.obj.trialdat(:,:,trix(i));
    cluix = randsample(size(plotdat,2),30,false);
    plotdat = plotdat(:,cluix);


    z = plotdat;
    z(~ix,:) = nan;

    for j = 1:size(plotdat,2)
        patchline(dat.obj.time,trialOffset + plotdat(:,j),'EdgeColor','k','EdgeAlpha',0.6,'LineWidth',2);
        plot(dat.obj.time,trialOffset + z(:,j),'r','LineWidth',1)
        trialOffset = trialOffset + dy;
    end


end
xlabel('Time (s) from go cue')
ylabel('Neurons')
title('Firing Rate')
xlim([dat.obj.time(15),dat.obj.time(end)]);
ax = gca;
ax.YTick = [];
ax.FontSize = 14;

align = 2.5;
sample = mode(dat.obj.bp.ev.sample) - align;
delay = mode(dat.obj.bp.ev.delay) - align;
xline(sample,'k--','LineWidth',1)
xline(delay,'k--','LineWidth',1)
xline(0,'k--','LineWidth',1)

hold off


%%
close all
temp = rez(sessix);

figure; imagesc(temp.covNull)
xlabel('Neurons')
ylabel('Neurons')
title('Quiet')
ax = gca;
ax.FontSize = 14;
% colormap(gray)
colormap(linspecer)
c1 = colorbar;


figure; imagesc(temp.covPotent)
xlabel('Neurons')
ylabel('Neurons')
title('Move')
ax = gca;
ax.FontSize = 14;
% colormap(hot)
colormap(linspecer)
c2 = colorbar;
c2.Limits = c1.Limits;

%%

close all

tempnull = dat.rez.N_null;
temppotent = dat.rez.N_potent;
tempme = dat.me.data;


% sumsqmag
ssmnull = sum(tempnull.^2,3,'omitnan');
ssmpotent = sum(temppotent.^2,3,'omitnan');

trix = 1:dat.obj.bp.Ntrials;
xlims = [-0.5 0.5];

[~,ix1] = min(abs(dat.obj.time + 0.3));
[~,ix2] = min(abs(dat.obj.time + 0.05));
[~,sorttrix] = sort(nanmean(tempme(ix1:ix2,:),1),'descend');
% sorttrix = trix;

cc = linspecer;
f = figure;
f.Position = [114   332   496   535];
imagesc(dat.obj.time,trix,tempme(:,sorttrix)','Interpolation','bilinear')
xlim(xlims)
xlabel('Time (s) from go cue')
ylabel('Trials')
colormap(cc);
c = colorbar;
c.Label.String = 'Motion Energy';

f = figure;
f.Position = [114   332   496   535];
imagesc(dat.obj.time,trix,ssmnull(:,sorttrix)','Interpolation','bilinear')
xlim(xlims)
xlabel('Time (s) from go cue')
ylabel('Trials')
colormap(cc);
c = colorbar;
% c.Limits = [0 200];
c.Label.String = 'Null Space Magnitude (a.u)';

f = figure;
f.Position = [114   332   496   535];
imagesc(dat.obj.time,trix,ssmpotent(:,sorttrix)','Interpolation','bilinear')
xlim(xlims)
xlabel('Time (s) from go cue')
ylabel('Trials')
colormap(cc);
c = colorbar;
% c.Limits = [0 200];
c.Label.String = 'Potent Space Magnitude (a.u)';

me_potent = corrcoef(tempme(:),ssmpotent(:));
me_null = corrcoef(tempme(:),ssmnull(:));

% make these plots pretty,
% calculate corrcoef for all sessions

%% for all sessions


close all

for sessix = 1:numel(rez)

    tempnull = rez(sessix).N_null;
    temppotent = rez(sessix).N_potent;
    tempme = me(sessix).data;


    % sumsqmag
    ssmnull = sum(tempnull.^2,3,'omitnan');
    ssmpotent = sum(temppotent.^2,3,'omitnan');
    me_potent = corrcoef(tempme(:),ssmpotent(:));
    me_null = corrcoef(tempme(:),ssmnull(:));

%     cc.null(sessix) = abs(me_null(1,2));
%     cc.potent(sessix) = abs(me_potent(1,2));

    cc.null(sessix) = me_null(1,2).^2;
    cc.potent(sessix) = me_potent(1,2).^2;


end


figure;
hold on;
scatter(cc.null,cc.potent,100,'MarkerFaceColor',[0 0 0] ./ 255,'MarkerEdgeColor','w')
xlim([0 0.25])
ylim([0 0.25])

ax = gca;
plot(ax.XLim,ax.YLim,'k--','LineWidth',1)
ax.FontSize = 14;

xlabel('$R^2 (ME,Null)$', 'Interpreter','latex','FontWeight','bold','FontSize',16)
ylabel('$R^2 (ME,Potent)$', 'Interpreter','latex','FontWeight','bold','FontSize',16);











