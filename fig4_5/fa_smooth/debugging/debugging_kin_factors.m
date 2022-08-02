
%% kinematics
six = 1; % session
temp = kinfeats{six}; % time,trials,feats
leg = kin(six).featLeg;

figure(1);
for i = 1:size(temp,2)
    clf
    imagesc(obj(six).time,1:numel(leg),squeeze(temp(:,i,:))');
    set(gca,'YTick',1:numel(leg));
    set(gca,'YTickLabel',leg);
    pause
end


%% kinematics reduced
six = 1; % session
temp = kinfeats_reduced{six}; % time,trials,feats

figure(1);
for i = 1:size(temp,2)
    clf
    imagesc(obj(six).time,1:size(temp,3),squeeze(temp(:,i,:))');
    set(gca,'YTick',1:size(temp,3));
    pause
end



%% neural data factors

six = 1; % session
temp = fa(six).falatents; % time,feats,trials

figure(1);
for i = 1:size(temp,3)
    clf
    imagesc(obj(six).time,1:size(temp,2),squeeze(temp(:,:,i))');
    set(gca,'YTick',1:size(temp,2));
    pause
end


for six = 1:numel(meta)
    % six = 1; % session
    temp = fa(six).falatents; % time,feats,trials
    tmp = cat(3,temp(:,:,params(six).trialid{2}),temp(:,:,params(six).trialid{3}));

    tmpmean{1} = mean(temp(:,:,params(six).trialid{2}),3);
    tmpmean{2} = mean(temp(:,:,params(six).trialid{3}),3);

    figure;
    for i = 1:size(temp,2)
        nexttile;
        imagesc(obj(six).time,1:size(tmp,3),squeeze(tmp(:,i,:))');
        hold on;
        yline(numel(params(six).trialid{2}),'k','LineWidth',2)
    end

    figure;
    for i = 1:size(temp,2)
        nexttile;
        plot(obj(six).time,tmpmean{1}(:,i),'b','LineWidth',2);
        hold on;
        plot(obj(six).time,tmpmean{2}(:,i),'r','LineWidth',2);
    end

    pause
end





%% null/potent single trials

six = 1;
null = rez(six).N_null;

null_rl = cat(2,null(:,params(six).trialid{2},:),null(:,params(six).trialid{3},:));

close all
figure;
for i = 1:size(null,3)
    nexttile;
    imagesc(obj(six).time,1:size(null_rl,2),null_rl(:,:,i)')
    hold on;
    yline(numel(params(six).trialid{2}),'k','LineWidth',2)
    colorbar;
    caxis([-2 2]);
end
sgtitle('null')


six = 1;
potent = rez(six).N_potent;

potent_rl = cat(2,potent(:,params(six).trialid{2},:),potent(:,params(six).trialid{3},:));

figure;
for i = 1:size(potent,3)
    nexttile;
    imagesc(obj(six).time,1:size(potent_rl,2),potent_rl(:,:,i)')
    hold on;
    yline(numel(params(six).trialid{2}),'k','LineWidth',2)
    colorbar;
    caxis([-2 2]);
end
sgtitle('potent')


%% null/potent mean trajectories

six = 1;

close all
figure;
for i = 1:size(rez(six).null_psth,2)
    nexttile;
    plot(obj(six).time,squeeze(rez(six).null_psth(:,i,1)),'b','LineWidth',2)
    hold on
    plot(obj(six).time,squeeze(rez(six).null_psth(:,i,2)),'r','LineWidth',2)
end
sgtitle('null')


figure;
for i = 1:size(rez(six).potent_psth,2)
    nexttile;
    plot(obj(six).time,squeeze(rez(six).potent_psth(:,i,1)),'b','LineWidth',2)
    hold on
    plot(obj(six).time,squeeze(rez(six).potent_psth(:,i,2)),'r','LineWidth',2)
end
sgtitle('potent')



%% cd

close all
clear early late go

six = 1;

early = cat(2,rez(six).cd.null.cdEarly_latent(:,params(six).trialid{2}),rez(six).cd.null.cdEarly_latent(:,params(six).trialid{3}));
late = cat(2,rez(six).cd.null.cdLate_latent(:,params(six).trialid{2}),rez(six).cd.null.cdLate_latent(:,params(six).trialid{3}));
go = cat(2,rez(six).cd.null.cdGo_latent(:,params(six).trialid{2}),rez(six).cd.null.cdGo_latent(:,params(six).trialid{3}));

figure;
nexttile;
imagesc(obj(six).time,1:size(early,2),early')
nexttile; 
imagesc(obj(six).time,1:size(late,2),late')
nexttile;
imagesc(obj(six).time,1:size(go,2),go')

early = cat(2,rez(six).cd.potent.cdEarly_latent(:,params(six).trialid{2}),rez(six).cd.potent.cdEarly_latent(:,params(six).trialid{3}));
late = cat(2,rez(six).cd.potent.cdLate_latent(:,params(six).trialid{2}),rez(six).cd.potent.cdLate_latent(:,params(six).trialid{3}));
go = cat(2,rez(six).cd.potent.cdGo_latent(:,params(six).trialid{2}),rez(six).cd.potent.cdGo_latent(:,params(six).trialid{3}));

figure;
nexttile;
imagesc(obj(six).time,1:size(early,2),early')
nexttile; 
imagesc(obj(six).time,1:size(late,2),late')
nexttile;
imagesc(obj(six).time,1:size(go,2),go')

% ----------------
clear early late go
early{1} = rez(six).cd.null.cdEarly_latent(:,params(six).trialid{2});
early{2} = rez(six).cd.null.cdEarly_latent(:,params(six).trialid{3});
late{1} = rez(six).cd.null.cdLate_latent(:,params(six).trialid{2});
late{2} = rez(six).cd.null.cdLate_latent(:,params(six).trialid{3});
go{1} = rez(six).cd.null.cdGo_latent(:,params(six).trialid{2});
go{2} = rez(six).cd.null.cdGo_latent(:,params(six).trialid{3});


alph = 0.5;
figure;
nexttile; hold on;
for i = 1:size(early{1},2)
    patchline(obj(six).time,early{1}(:,i),'EdgeColor','b','EdgeAlpha',alph)
end
for i = 1:size(early{2},2)
    patchline(obj(six).time,early{2}(:,i),'EdgeColor','r','EdgeAlpha',alph)
end
nexttile; hold on;
for i = 1:size(late{1},2)
    patchline(obj(six).time,late{1}(:,i),'EdgeColor','b','EdgeAlpha',alph)
end
for i = 1:size(late{2},2)
    patchline(obj(six).time,late{2}(:,i),'EdgeColor','r','EdgeAlpha',alph)
end
nexttile; hold on;
for i = 1:size(go{1},2)
    patchline(obj(six).time,go{1}(:,i),'EdgeColor','b','EdgeAlpha',alph)
end
for i = 1:size(go{2},2)
    patchline(obj(six).time,go{2}(:,i),'EdgeColor','r','EdgeAlpha',alph)
end

clear early late go
early{1} = rez(six).cd.potent.cdEarly_latent(:,params(six).trialid{2});
early{2} = rez(six).cd.potent.cdEarly_latent(:,params(six).trialid{3});
late{1} = rez(six).cd.potent.cdLate_latent(:,params(six).trialid{2});
late{2} = rez(six).cd.potent.cdLate_latent(:,params(six).trialid{3});
go{1} = rez(six).cd.potent.cdGo_latent(:,params(six).trialid{2});
go{2} = rez(six).cd.potent.cdGo_latent(:,params(six).trialid{3});

figure;
nexttile; hold on;
for i = 1:size(early{1},2)
    patchline(obj(six).time,early{1}(:,i),'EdgeColor','b','EdgeAlpha',alph)
end
for i = 1:size(early{2},2)
    patchline(obj(six).time,early{2}(:,i),'EdgeColor','r','EdgeAlpha',alph)
end
nexttile; hold on;
for i = 1:size(late{1},2)
    patchline(obj(six).time,late{1}(:,i),'EdgeColor','b','EdgeAlpha',alph)
end
for i = 1:size(late{2},2)
    patchline(obj(six).time,late{2}(:,i),'EdgeColor','r','EdgeAlpha',alph)
end
nexttile; hold on;
for i = 1:size(go{1},2)
    patchline(obj(six).time,go{1}(:,i),'EdgeColor','b','EdgeAlpha',alph)
end
for i = 1:size(go{2},2)
    patchline(obj(six).time,go{2}(:,i),'EdgeColor','r','EdgeAlpha',alph)
end




