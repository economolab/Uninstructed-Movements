function nullSpaceCD(rez,obj,params,times)

normRange = [0 1];

% concatenate latents, find mean and stderror
cdEarly{1} = normalizeInRange(rez(1).cd.null.cdEarly_latent(:,params(1).trialid{2}),normRange);
cdEarly{2} = normalizeInRange(rez(1).cd.null.cdEarly_latent(:,params(1).trialid{3}),normRange);
cdLate{1} = normalizeInRange(rez(1).cd.null.cdLate_latent(:,params(1).trialid{2}),normRange);
cdLate{2} = normalizeInRange(rez(1).cd.null.cdLate_latent(:,params(1).trialid{3}),normRange);
cdGo{1} = normalizeInRange(rez(1).cd.null.cdGo_latent(:,params(1).trialid{2}),normRange);
cdGo{2} = normalizeInRange(rez(1).cd.null.cdGo_latent(:,params(1).trialid{3}),normRange);
for i = 2:numel(rez)
    for j = 1:2
        cdEarly{j} = normalizeInRange(cat(2,cdEarly{j},rez(i).cd.null.cdEarly_latent(:,params(i).trialid{j+1})),normRange);
        cdLate{j} = normalizeInRange(cat(2,cdLate{j},rez(i).cd.null.cdLate_latent(:,params(i).trialid{j+1})),normRange);
        cdGo{j} = normalizeInRange(cat(2,cdGo{j},rez(i).cd.null.cdGo_latent(:,params(i).trialid{j+1})),normRange);
    end
end

cdEarly_latent_mean = [nanmean(cdEarly{1},2) nanmean(cdEarly{2},2)];
cdLate_latent_mean = [nanmean(cdLate{1},2) nanmean(cdLate{2},2)];
cdGo_latent_mean = [nanmean(cdGo{1},2) nanmean(cdGo{2},2)];

cdEarly_latent_error = [nanstd(cdEarly{1},[],2) nanstd(cdEarly{2},[],2)] ./ (numel(rez)); % std error
cdLate_latent_error = [nanstd(cdLate{1},[],2) nanstd(cdLate{2},[],2)] ./ numel(rez); % std error
cdGo_latent_error = [nanstd(cdGo{1},[],2) nanstd(cdGo{2},[],2)] ./ numel(rez); % std error


fns = patternMatchCellArray(fieldnames(rez(1).cd.null),{'mode'},'all');

% ve
% 
% ve = zeros(numel(fns),numel(rez));
% for i = 1:numel(rez)
%     ve(1,i) = rez(i).cd.null.ve.early;
%     ve(2,i) = rez(i).cd.null.ve.late;
%     ve(3,i) = rez(i).cd.null.ve.go;
% end
% meanve = mean(ve,2);


close all
clrs = getColors();
lw = 6;
alph = 0.5;

sm = 1;

sample = mode(obj(1).bp.ev.sample - obj(1).bp.ev.(params(1).alignEvent));
delay = mode(obj(1).bp.ev.delay - obj(1).bp.ev.(params(1).alignEvent));


sav = 0;
for i = 1:numel(fns)
    f(i) = figure; f(i).Position = [222   711-(320*(i-1))   468   274];
    ax = axes(f(i)); hold on;
    tempmean = mySmooth(eval([fns{i}(1:end-5) '_latent_mean']),sm);
    temperror = mySmooth(eval([fns{i}(1:end-5) '_latent_error']),sm);
    shadedErrorBar(obj(1).time,tempmean(:,1),temperror(:,1),{'Color',clrs.rhit,'LineWidth',lw},alph, ax)
    shadedErrorBar(obj(1).time,tempmean(:,2),temperror(:,2),{'Color',clrs.lhit,'LineWidth',lw},alph, ax)
    
    xlim([obj(1).time(10);obj(1).time(end)])
%     ylims = [min(min(tempmean)), max(max(tempmean))+5];
%     ylim(ylims);
    
%     title([fns{i} '| MeanVE=' num2str(meanve(i))],'Interpreter','none')
    xlabel('Time (s) from go cue')
    ylabel('Activity (a.u.)')
    ax = gca;
    ax.FontSize = 20;
    
    xline(sample,'k--','LineWidth',2)
    xline(delay,'k--','LineWidth',2)
    xline(0,'k--','LineWidth',2)
    
    curmodename = fns{i};
    timefns = fieldnames(times);
    mask = strfind(timefns,lower(curmodename(3:end-5)));
    ix = cellfun(@(x) isempty(x),mask,'UniformOutput',false);
    ix = ~cell2mat(ix);
    shadetimes = obj(1).time(times.(timefns{ix}));
    x = [shadetimes(1)  shadetimes(end) shadetimes(end) shadetimes(1)];
%     y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
    y = [-1 -1 1 1];
    fl = fill(x,y,'r','FaceColor',[93, 121, 148]./255);
    fl.FaceAlpha = 0.3;
    fl.EdgeColor = 'none';
    
    ylim(normRange);
    
    
    if sav
        pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/kinNullSpace/null/cd';
        fn = [fns{i}];
        mysavefig(f(i),pth,fn);
    end
    
end

