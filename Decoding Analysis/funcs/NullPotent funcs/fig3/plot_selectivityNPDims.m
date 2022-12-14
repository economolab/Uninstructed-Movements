close all

cols = getColors;

for sessix = 1:numel(meta)
    trix{1} = params(sessix).trialid{2};
    trix{2} = params(sessix).trialid{3};

    for i = 1:numel(trix)
        null{i} = squeeze(mean(rez(sessix).N_null(:,trix{i},:),2));
        potent{i} = squeeze(mean(rez(sessix).N_potent(:,trix{i},:),2));
    end

    selnull = ((null{1} - null{2}).^2);
    selpotent = ((potent{1} - potent{2}).^2);

    sel.null(:,sessix) = mean(selnull,2);
    sel.potent(:,sessix) = mean(selpotent,2);

end

%%
close all
lw = 2;
alph = 0.2;

f = figure;
ax = gca;
hold on;

means = mean(sel.null,2);
errs = std(sel.null,[],2) ./ sqrt(size(sel.null,2));
% CI95 = tinv([0.025 0.975], size(sel.null,2)-1);   
% errs = means*CI95(2);
shadedErrorBar(obj(1).time,means,errs,{'Color',cols.null,'LineWidth',lw},alph,ax)

means = mean(sel.potent,2);
errs = std(sel.potent,[],2) ./ sqrt(size(sel.potent,2));
shadedErrorBar(obj(1).time,means,errs,{'Color',cols.potent,'LineWidth',lw},alph,ax)

sample = mode(obj(sessix).bp.ev.sample) - 2.5;
delay = mode(obj(sessix).bp.ev.delay) - 2.5;
xline(0,'k--')
xline(sample,'k--')
xline(delay,'k--')
xlim([obj(sessix).time(5) 0])
% h = breakyaxis([1 10.5]);



