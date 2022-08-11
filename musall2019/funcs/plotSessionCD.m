function plotSessionCD(rez,obj,params,times,meta)

fns = patternMatchCellArray(fieldnames(rez.cd.null),{'mode'},'all');


clrs = getColors();
lw = 1;
alph = 0.3;

sm = 51;

sample = mode(obj(1).bp.ev.sample - obj(1).bp.ev.(params(1).alignEvent));
delay = mode(obj(1).bp.ev.delay - obj(1).bp.ev.(params(1).alignEvent));

% % single trials
% sav = 0;
% f = figure;
% f.Position = [334    99   423   875];
% for i = 1:numel(fns)
%     ax = subplot(numel(fns),1,i); hold on;
%     temp = rez.cd.null.([fns{i}(1:end-5) '_latent']);
%     for j = 1:numel(params.trialid{2})
%         toplot = temp(:,params.trialid{2}(j));
%         toplot = mySmooth(toplot,sm);
%         patchline(obj.time,toplot,'EdgeColor',clrs.rhit,'EdgeAlpha',alph,'LineWidth',lw)
%     end
%     for j = 1:numel(params.trialid{3})
%         toplot = temp(:,params.trialid{3}(j));
%         toplot = mySmooth(toplot,sm);
%         patchline(obj.time,toplot,'EdgeColor',clrs.lhit,'EdgeAlpha',alph,'LineWidth',lw)
%     end
%     
% 
% %     plot(obj.time,temp(:,params.trialid{2}),'Color',clrs.rhit,'LineWidth',lw)
% %     plot(obj.time,temp(:,params.trialid{3}),'Color',clrs.lhit,'LineWidth',lw)
%     
%     xlim([obj(1).time(10);obj(1).time(end)])
% %     ylims = [min(min(tempmean)), max(max(tempmean))+5];
% %     ylim(ylims);
%     
% %     title([fns{i} '| MeanVE=' num2str(meanve(i))],'Interpreter','none')
%     xlabel('Time (s) from go cue')
%     ylabel('Activity (a.u.)')
%     ax = gca;
%     ax.FontSize = 20;
%     
%     xline(sample,'k--','LineWidth',2)
%     xline(delay,'k--','LineWidth',2)
%     xline(0,'k--','LineWidth',2)
%     
%     curmodename = fns{i};
%     timefns = fieldnames(times);
%     mask = strfind(timefns,lower(curmodename(3:end-5)));
%     ix = cellfun(@(x) isempty(x),mask,'UniformOutput',false);
%     ix = ~cell2mat(ix);
%     shadetimes = obj(1).time(times.(timefns{ix}));
%     x = [shadetimes(1)  shadetimes(end) shadetimes(end) shadetimes(1)];
%     y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
% %     y = [-1 -1 1 1];
%     fl = fill(x,y,'r','FaceColor',[93, 121, 148]./255);
%     fl.FaceAlpha = 0.3;
%     fl.EdgeColor = 'none';
%     
% %     ylim(normRange);    
%     
%     if sav
%         pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/kinNullSpace/null/cd';
%         fn = [fns{i}];
%         mysavefig(f(i),pth,fn);
%     end
%     
% end


% trial averaged
sav = 0;
f(1) = figure;
f(1).Position = [334    99   423   875];
f(2) = figure;
f(2).Position = [760    96   423   875];
for i = 1:numel(fns)
    set(0,'CurrentFigure',f(1));
    ax1 = subplot(numel(fns),1,i); hold on;
    temp = rez.cd.null.([fns{i}(1:end-5) '_latent']);

    toplot = temp(:,params.trialid{2});
    toplot = mean(toplot,2);
    toplot = mySmooth(toplot,sm);
    plot(obj.time,toplot,'Color',clrs.rhit,'LineWidth',3)

    toplot = temp(:,params.trialid{3});
    toplot = mean(toplot,2);
    toplot = mySmooth(toplot,sm);
    plot(obj.time,toplot,'Color',clrs.lhit,'LineWidth',3)
    
    xlim([obj(1).time(20);obj(1).time(end)])
%     ylims = [min(min(tempmean)), max(max(tempmean))+5];
%     ylim(ylims);
    
%     title([fns{i} '| MeanVE=' num2str(meanve(i))],'Interpreter','none')
    sgtitle(['Null | ' meta.anm ' ' meta.date])
    xlabel('Time (s) from go cue')
    ylabel('Activity (a.u.)')
    ax1.FontSize = 20;
    
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
    y = [ax1.YLim(1) ax1.YLim(1) ax1.YLim(2) ax1.YLim(2)];
%     y = [-1 -1 1 1];
    fl = fill(x,y,'r','FaceColor',[93, 121, 148]./255);
    fl.FaceAlpha = 0.3;
    fl.EdgeColor = 'none';

    set(0,'CurrentFigure',f(2));
    ax2 = subplot(numel(fns),1,i); hold on;
    temp = rez.cd.potent.([fns{i}(1:end-5) '_latent']);

    toplot = temp(:,params.trialid{2});
    toplot = mean(toplot,2);
    toplot = mySmooth(toplot,sm);
    plot(obj.time,toplot,'Color',clrs.rhit,'LineWidth',3)

    toplot = temp(:,params.trialid{3});
    toplot = mean(toplot,2);
    toplot = mySmooth(toplot,sm);
    plot(obj.time,toplot,'Color',clrs.lhit,'LineWidth',3)
    
    xlim([obj(1).time(20);obj(1).time(end)])
%     ylims = [min(min(tempmean)), max(max(tempmean))+5];
%     ylim(ylims);
    
%     title([fns{i} '| MeanVE=' num2str(meanve(i))],'Interpreter','none')
    sgtitle(['Potent | ' meta.anm ' ' meta.date])
    xlabel('Time (s) from go cue')
    ylabel('Activity (a.u.)')
    ax2.FontSize = 20;
    
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
    y = [ax2.YLim(1) ax2.YLim(1) ax2.YLim(2) ax2.YLim(2)];
%     y = [-1 -1 1 1];
    fl = fill(x,y,'r','FaceColor',[93, 121, 148]./255);
    fl.FaceAlpha = 0.3;
    fl.EdgeColor = 'none';


%     'a'

    ylims = [ax1.YLim; ax2.YLim];
    [~,ix] = max(abs(ylims(:,1)));
    newylims(1) = ylims(ix,1);
    [~,ix] = max(abs(ylims(:,2)));
    newylims(2) = ylims(ix,2);
    ax1.YLim = newylims;
    ax2.YLim = newylims;

    
end




% 'a'
end