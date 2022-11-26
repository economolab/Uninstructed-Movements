function plotCorrelation_Cddelay_ME(obj,params,rez,me,meta,cond2use)

cols = linspecer(numel(obj));

r2 = [];
for sessix = 1:numel(rez)
    % only use 400 ms up to go cue for correlations
    [~,ix1] = min(abs(obj(sessix).time - -0.4));
    [~,ix2] = min(abs(obj(sessix).time - 0));

    trix = cell2mat({params(sessix).trialid{cond2use}}');

    cddelayix = find(ismember( rez(sessix).cd_labels,'late'));
    plotcd = rez(sessix).cd_proj_trialdat(ix1:ix2,trix,cddelayix);
    avgcd = mean(plotcd,1);

    plotme = me(sessix).data(ix1:ix2,trix);
    avgme = mean(plotme,1);

    r2 = [r2 ; corr(avgme',avgcd')];

end

f = figure; hold on
histogram(r2,10,'EdgeColor',[29, 50, 84]./255,'FaceColor',[62, 109, 184]./255, 'FaceAlpha',0.6)


% null distribution, shuffle avgcd and avgme and correlate
r2null = [];
for sessix = 1:numel(rez)
    % only use 400 ms up to go cue for correlations
    [~,ix1] = min(abs(obj(sessix).time - -0.4));
    [~,ix2] = min(abs(obj(sessix).time - 0));

    trix = cell2mat({params(sessix).trialid{cond2use}}');

    cddelayix = find(ismember( rez(sessix).cd_labels,'late'));
    plotcd = rez(sessix).cd_proj_trialdat(ix1:ix2,trix,cddelayix);
    avgcd = mean(plotcd,1);

    plotme = me(sessix).data(ix1:ix2,trix);
    avgme = mean(plotme,1);

    r2null = [r2null ; corr(randsample(avgme,numel(avgcd),false)',randsample(avgcd,numel(avgcd),false)')];

end

histogram(r2null,4,'EdgeColor',[0 0 0]./255,'FaceColor',[40 40 40]./255, 'FaceAlpha',0.6)
% h = kstest2(x1,x2)
xlabel('Correlation b/w CDdelay and ME on single trials (avg last 400 ms of delay)')


end