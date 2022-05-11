function plotNullPotentProjections(obj,dat,rez,params)

% for lfads data
rhit = find(obj.bp.R & obj.bp.hit & ~obj.bp.autowater & ~obj.bp.early);
lhit = find(obj.bp.L & obj.bp.hit & ~obj.bp.autowater & ~obj.bp.early);
mask = ismember(dat.trials,rhit);
rhit = find(mask);
mask = ismember(dat.trials,lhit);
lhit = find(mask);


align = mode(obj.bp.ev.(params.alignEvent));
sample = mode(obj.bp.ev.sample) - align;
delay = mode(obj.bp.ev.delay) - align;

f = figure(1); sgtitle(['Null Space Projections | %VE=' num2str(rez.varexp_null*100)])
f.Position = [-792    58   431   635];
for i = 1:size(rez.N_null,2) % num null dims
    ax = subplot(size(rez.N_null,2),1,i); hold on
    for j = 1:numel(rhit)
        patchline(obj.time,squeeze(rez.N_null(:,i,rhit(j))),'EdgeColor','b','EdgeAlpha',0.4,'LineWidth',1.5);
    end
    for j = 1:numel(lhit)
        patchline(obj.time,squeeze(rez.N_null(:,i,lhit(j))),'EdgeColor','r','EdgeAlpha',0.4,'LineWidth',1.5);
    end
    xlim([obj.time(20) obj.time(end)])
    
    xline(sample,'k--','LineWidth',2)
    xline(delay,'k--','LineWidth',2)
    xline(0,'k--','LineWidth',2)
    
    
    if i~=size(rez.N_null,2)
        ax.XTick = [];
    else
        ax.XLabel.String = 'Time (s) from go cue';
    end
    
    ax.FontSize = 20;
    
    hold off
end

f = figure(2); sgtitle(['Potent Space Projections | %VE=' num2str(rez.varexp_potent*100)])
f.Position = [-792    58   431   635];
for i = 1:size(rez.N_potent,2) % num null dims
    ax = subplot(size(rez.N_potent,2),1,i); hold on
    for j = 1:numel(rhit)
        patchline(obj.time,squeeze(rez.N_potent(:,i,rhit(j))),'EdgeColor','b','EdgeAlpha',0.4,'LineWidth',1.5);
    end
    for j = 1:numel(lhit)
        patchline(obj.time,squeeze(rez.N_potent(:,i,lhit(j))),'EdgeColor','r','EdgeAlpha',0.4,'LineWidth',1.5);
    end
    xlim([obj.time(20) obj.time(end)])

    
    xline(sample,'k--','LineWidth',2)
    xline(delay,'k--','LineWidth',2)
    xline(0,'k--','LineWidth',2)
    
     if i~=size(rez.N_potent,2)
        ax.XTick = [];
    else
        ax.XLabel.String = 'Time (s) from go cue';
    end
    
    ax.FontSize = 20;
    
    hold off
end

% % plot means
% figure(3); sgtitle(['Mean Null Space Projections | %VE=' num2str(rez.varexp_null*100)])
% for i = 1:size(rez.N_null,2) % num null dims
%     subplot(size(rez.N_null,2),1,i); hold on
%     plot(obj.time,mean(squeeze(rez.N_null(:,i,rhit)),2),'LineWidth',2)
%     plot(obj.time,mean(squeeze(rez.N_null(:,i,lhit)),2),'LineWidth',2)
%     xlim([obj.time(20) obj.time(end)])
%     hold off
% end
% 
% figure(4); sgtitle(['Mean Potent Space Projections | %VE=' num2str(rez.varexp_potent*100)])
% for i = 1:size(rez.N_potent,2) % num null dims
%     subplot(size(rez.N_potent,2),1,i); hold on
%     plot(obj.time,mean(squeeze(rez.N_potent(:,i,rhit)),2),'LineWidth',2)
%     plot(obj.time,mean(squeeze(rez.N_potent(:,i,lhit)),2),'LineWidth',2)
%     xlim([obj.time(20) obj.time(end)])
%     hold off
% end


% figure(3); sgtitle('Null Space Projections')
% for i = 1:size(rez.N_null,2) % num null dims
%     subplot(size(rez.N_null,2),1,i);
%     tempdat = squeeze(rez.N_null(:,i,[rhit;lhit]))';
%     imagesc(tempdat)
% end
% 
% figure(4); sgtitle('Potent Space Projections')
% for i = 1:size(rez.N_potent,2) % num null dims
%     subplot(size(rez.N_potent,2),1,i);
%     tempdat = squeeze(rez.N_potent(:,i,[rhit;lhit]))';
%     imagesc(squeeze(rez.N_potent(:,i,:))')
% end
% 
% 
% rmiss = find(obj.bp.R & obj.bp.miss & ~obj.bp.autowater & ~obj.bp.early);
% lmiss = find(obj.bp.L & obj.bp.miss & ~obj.bp.autowater & ~obj.bp.early);
% mask = ismember(dat.trials,rmiss);
% rhit = find(mask);
% mask = ismember(dat.trials,lmiss);
% lhit = find(mask);
% 
% figure(5); sgtitle('Null Space Projections')
% for i = 1:size(rez.N_null,2) % num null dims
%     subplot(size(rez.N_null,2),1,i); hold on
%     for j = 1:numel(rhit)
%         patchline(obj.time,squeeze(rez.N_null(:,i,rhit(j))),'EdgeColor','b','EdgeAlpha',0.4,'LineWidth',1.5);
%     end
%     for j = 1:numel(lhit)
%         patchline(obj.time,squeeze(rez.N_null(:,i,lhit(j))),'EdgeColor','r','EdgeAlpha',0.4,'LineWidth',1.5);
%     end
%     xlim([obj.time(100) obj.time(end)])
%     hold off
% end
% 
% figure(6); sgtitle('Potent Space Projections')
% for i = 1:size(rez.N_potent,2) % num null dims
%     subplot(size(rez.N_potent,2),1,i); hold on
%     for j = 1:numel(rhit)
%         patchline(obj.time,squeeze(rez.N_potent(:,i,rhit(j))),'EdgeColor','b','EdgeAlpha',0.4,'LineWidth',1.5);
%     end
%     for j = 1:numel(lhit)
%         patchline(obj.time,squeeze(rez.N_potent(:,i,lhit(j))),'EdgeColor','r','EdgeAlpha',0.4,'LineWidth',1.5);
%     end
%     xlim([obj.time(100) obj.time(end)])
%     hold off
% end
% 
% 
% rmiss = find(obj.bp.R & obj.bp.hit & obj.bp.autowater & ~obj.bp.early);
% lmiss = find(obj.bp.L & obj.bp.hit & obj.bp.autowater & ~obj.bp.early);
% mask = ismember(dat.trials,rmiss);
% rhit = find(mask);
% mask = ismember(dat.trials,lmiss);
% lhit = find(mask);
% 
% figure(7); sgtitle('Null Space Projections')
% for i = 1:size(rez.N_null,2) % num null dims
%     subplot(size(rez.N_null,2),1,i); hold on
%     for j = 1:numel(rhit)
%         patchline(obj.time,squeeze(rez.N_null(:,i,rhit(j))),'EdgeColor','b','EdgeAlpha',0.4,'LineWidth',1.5);
%     end
%     for j = 1:numel(lhit)
%         patchline(obj.time,squeeze(rez.N_null(:,i,lhit(j))),'EdgeColor','r','EdgeAlpha',0.4,'LineWidth',1.5);
%     end
%     xlim([obj.time(100) obj.time(end)])
%     hold off
% end
% 
% figure(8); sgtitle('Potent Space Projections')
% for i = 1:size(rez.N_potent,2) % num null dims
%     subplot(size(rez.N_potent,2),1,i); hold on
%     for j = 1:numel(rhit)
%         patchline(obj.time,squeeze(rez.N_potent(:,i,rhit(j))),'EdgeColor','b','EdgeAlpha',0.4,'LineWidth',1.5);
%     end
%     for j = 1:numel(lhit)
%         patchline(obj.time,squeeze(rez.N_potent(:,i,lhit(j))),'EdgeColor','r','EdgeAlpha',0.4,'LineWidth',1.5);
%     end
%     xlim([obj.time(100) obj.time(end)])
%     hold off
% end