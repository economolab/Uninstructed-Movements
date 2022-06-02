function optimization_plots(rez,obj,dat,params)


% for lfads data
rhit = params.trialid{2};
lhit = params.trialid{3};

trials = [rhit; lhit];
rmask = 1:numel(rhit);
lmask = (numel(rhit)+1):(numel(rhit)+1 + numel(lhit) - 1);

align = mode(obj.bp.ev.(params.alignEvent));
sample = mode(obj.bp.ev.sample) - align;
delay = mode(obj.bp.ev.delay) - align;

try
    N_null = rez.optim.N_null;
    N_potent = rez.optim.N_potent;
catch
    N_null = rez.N_null;
    N_potent = rez.N_potent;
end

lw = 6;
clrs = getColors();

f = figure; clf;  sgtitle(['Null Space Projections | %VE=' num2str(rez.varexp_null*100)])
f.Position = [-856  -124   431   937];
for i = 1:size(N_null,2) % num null dims
    ax = subplot(size(N_null,2),1,i); hold on
    plot(obj.time,mean(squeeze(N_null(:,i,rmask)),2),'Color',clrs.rhit,'LineWidth',lw)
    plot(obj.time,mean(squeeze(N_null(:,i,lmask)),2),'Color',clrs.lhit,'LineWidth',lw)
    
%     for j = 1:numel(rmask)
%         patchline(obj.time,squeeze(N_null(:,i,rmask(j))),'EdgeColor','b','EdgeAlpha',0.4,'LineWidth',1.5);
%     end
%     for j = 1:numel(lmask)
%         patchline(obj.time,squeeze(N_null(:,i,lmask(j))),'EdgeColor','r','EdgeAlpha',0.4,'LineWidth',1.5);
%     end
    
    xline(sample,'k--','LineWidth',2)
    xline(delay,'k--','LineWidth',2)
    xline(0,'k--','LineWidth',2)
    
    
    if i~=size(N_null,2)
        ax.XTick = [];
    else
        ax.XLabel.String = 'Time (s) from go cue';
    end
    
    ax.FontSize = 30;
        
    xlim([obj.time(20) obj.time(end)])
    hold off
end

f = figure;  clf; sgtitle(['Potent Space Projections | %VE=' num2str(rez.varexp_potent*100)])
f.Position = [-856  -124   431   937];
for i = 1:size(N_potent,2) % num null dims
    ax = subplot(size(N_potent,2),1,i); hold on
    
    plot(obj.time,mean(squeeze(N_potent(:,i,rmask)),2),'Color',clrs.rhit,'LineWidth',lw)
    plot(obj.time,mean(squeeze(N_potent(:,i,lmask)),2),'Color',clrs.lhit,'LineWidth',lw)
    
%     for j = 1:numel(rmask)
%         patchline(obj.time,squeeze(N_potent(:,i,rmask(j))),'EdgeColor','b','EdgeAlpha',0.4,'LineWidth',1.5);
%     end
%     for j = 1:numel(lmask)
%         patchline(obj.time,squeeze(N_potent(:,i,lmask(j))),'EdgeColor','r','EdgeAlpha',0.4,'LineWidth',1.5);
%     end
    
    xline(sample,'k--','LineWidth',2)
    xline(delay,'k--','LineWidth',2)
    xline(0,'k--','LineWidth',2)
    
    if i~=size(N_null,2)
        ax.XTick = [];
    else
        ax.XLabel.String = 'Time (s) from go cue';
    end
    
    ax.FontSize = 30;
    
    xlim([obj.time(20) obj.time(end)])
    hold off
end