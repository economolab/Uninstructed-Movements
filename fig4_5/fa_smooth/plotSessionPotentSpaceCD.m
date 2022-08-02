function plotSessionPotentSpaceCD(rez,obj,params,times)

fns = patternMatchCellArray(fieldnames(rez.cd.potent),{'mode'},'all');


clrs = getColors();
lw = 1;
alph = 0.3;

sm = 31;

sample = mode(obj(1).bp.ev.sample - obj(1).bp.ev.(params(1).alignEvent));
delay = mode(obj(1).bp.ev.delay - obj(1).bp.ev.(params(1).alignEvent));


sav = 0;
f = figure;
f.Position = [760    96   423   875];
for i = 1:numel(fns)
    ax = subplot(numel(fns),1,i); hold on;
    temp = rez.cd.potent.([fns{i}(1:end-5) '_latent']);
    for j = 1:numel(params.trialid{2})
        toplot = temp(:,params.trialid{2}(j));
        toplot = mySmooth(toplot,sm);
        patchline(obj.time,toplot,'EdgeColor',clrs.rhit,'EdgeAlpha',alph,'LineWidth',lw)
    end
    for j = 1:numel(params.trialid{3})
        toplot = temp(:,params.trialid{3}(j));
        toplot = mySmooth(toplot,sm);
        patchline(obj.time,toplot,'EdgeColor',clrs.lhit,'EdgeAlpha',alph,'LineWidth',lw)
    end

%     plot(obj.time,temp(:,params.trialid{2}),'Color',clrs.rhit,'LineWidth',lw)
%     plot(obj.time,temp(:,params.trialid{3}),'Color',clrs.lhit,'LineWidth',lw)
    
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
    y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
%     y = [-1 -1 1 1];
    fl = fill(x,y,'r','FaceColor',[93, 121, 148]./255);
    fl.FaceAlpha = 0.3;
    fl.EdgeColor = 'none';
    
%     ylim(normRange);
    
    
    if sav
        pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/kinNullSpace/null/cd';
        fn = [fns{i}];
        mysavefig(f(i),pth,fn);
    end
    
end


% trial averaged
sav = 0;
f = figure;
f.Position = [760    96   423   875];
for i = 1:numel(fns)
    ax = subplot(numel(fns),1,i); hold on;
    temp = rez.cd.potent.([fns{i}(1:end-5) '_latent']);

    toplot = temp(:,params.trialid{2});
    toplot = mean(toplot,2);
    toplot = mySmooth(toplot,sm);
    plot(obj.time,toplot,'Color',clrs.rhit,'LineWidth',3)

    toplot = temp(:,params.trialid{3});
    toplot = mean(toplot,2);
    toplot = mySmooth(toplot,sm);
    plot(obj.time,toplot,'Color',clrs.lhit,'LineWidth',3)

%     toplot = temp(:,params.trialid{4});
%     toplot = mean(toplot,2);
%     toplot = mySmooth(toplot,sm+31);
%     plot(obj.time,toplot,'Color',clrs.rmiss,'LineWidth',2)
% 
%     toplot = temp(:,params.trialid{5});
%     toplot = mean(toplot,2);
%     toplot = mySmooth(toplot,sm+31);
%     plot(obj.time,toplot,'Color',clrs.lmiss,'LineWidth',2)

    
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
    y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
%     y = [-1 -1 1 1];
    fl = fill(x,y,'r','FaceColor',[93, 121, 148]./255);
    fl.FaceAlpha = 0.3;
    fl.EdgeColor = 'none';
    
%     ylim(normRange);    
    
    if sav
        pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/kinNullSpace/null/cd';
        fn = [fns{i}];
        mysavefig(f(i),pth,fn);
    end
    
end





% 'a'
end