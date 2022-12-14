
% get me, null, potent, aligned to first lick and last lick

%%

params_lastLick.alignEvent          = 'lastLick'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params_lastLick.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params_lastLick.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params_lastLick.lowFR               = 0; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params_lastLick.condition(1)     = {'(hit|miss|no)'};                             % all trials
params_lastLick.condition(end+1) = {'R&hit&~stim.enable&~autowater'};             % right hits, no stim, aw off
params_lastLick.condition(end+1) = {'L&hit&~stim.enable&~autowater'};             % left hits, no stim, aw off
params_lastLick.condition(end+1) = {'R&miss&~stim.enable&~autowater'};            % error right, no stim, aw off
params_lastLick.condition(end+1) = {'L&miss&~stim.enable&~autowater'};            % error left, no stim, aw off

params_lastLick.tmin = -3.5;
params_lastLick.tmax = 1.5;
params_lastLick.dt = 1/100;

% smooth with causal gaussian kernel
params_lastLick.smooth = 15;

% cluster qualities to use
params_lastLick.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality
% params_lastLick.quality = {'Excellent','Great','Good','Fair','Multi'};


params_lastLick.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

params_lastLick.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

params_lastLick.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

params_lastLick.advance_movement = 0;

params_lastLick.probe = {meta.probe};


% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% params_lastLick (struct array) - one entry per session
% ----------------------------------------------
[obj_lastLick,params_lastLick] = loadSessionData(meta,params_lastLick);

% ------------------------------------------
% -- Motion Energy --
% me (struct array) - one entry per session
% ------------------------------------------
for sessix = 1:numel(meta)
    me_lastLick(sessix) = loadMotionEnergy(obj_lastLick(sessix), meta(sessix), params_lastLick(sessix), datapth);
end

% project data onto null and potent spaces found from firstLick
rez_lastLick = rez;
for sessix = 1:numel(meta)
    %     clu = params(sessix).cluid;
    %     clu_lastLick = params_lastLick(sessix).cluid;
    %     [numel(clu) numel(clu_lastLick)]
    trialdat_zscored = zscore_singleTrialNeuralData(obj_lastLick(sessix).trialdat);
    rez_lastLick(sessix) = projectNP(params_lastLick(sessix).trialid([1]),trialdat_zscored,rez_lastLick(sessix));
end



%%

close all
for sessix = 11%1:numel(meta)

    temprez = rez(sessix);
    tempme = me(sessix);

    %     trix = 1:obj(sessix).bp.Ntrials;
    trix = find(~obj(sessix).bp.no);

    potent = temprez.N_potent;
    potent = potent(:,trix,:);
    potent = sum(potent.^2,3);
    dims = size(potent);
    potent = normalize(potent(:),'range',[0 1]);
    potent = reshape(potent,dims(1),dims(2));

    plotme = tempme.data(:,trix);
    dims = size(plotme);
    plotme = normalize(plotme(:),'range',[0 1]);
    plotme = reshape(plotme,dims(1),dims(2));


    % sort
    ixs = 180:245;
    [~,sortix] = sort(mean(plotme(ixs,:),1),'descend');


    xlims = [-0.5 0.5];

    t = obj(sessix).time;
    trials = 1:size(plotme,2);

    f = figure;
    f.Position = [717   407   371   538];

%         subplot(2,2,1)
    imagesc(t,trials,plotme(:,sortix)')
    colormap(linspecer)
    xlim(xlims)
    colorbar;
    title('me first lick')

    f = figure;
    f.Position = [717   407   371   538]; 
%     subplot(2,2,2)
    imagesc(t,trials,potent(:,sortix)')
    colormap(linspecer)
    xlim(xlims)
    colorbar;
    title('potent first lick')

    % last lick
    temprez = rez_lastLick(sessix);
    tempme = me_lastLick(sessix);

    potent = temprez.N_potent;
    potent = potent(:,trix,:);
    potent = sum(potent.^2,3);
    dims = size(potent);
    potent = normalize(potent(:),'range',[0 1]);
    potent = reshape(potent,dims(1),dims(2));

    plotme = tempme.data(:,trix);
    dims = size(plotme);
    plotme = plotme(:);
    plotme(plotme==0) = nan;
    plotme = fillmissing(plotme,"movmedian",10000);
    plotme = normalize(plotme,'range',[0 1]);
    plotme = reshape(plotme,dims(1),dims(2));


    t = obj_lastLick(sessix).time;
    trials = 1:size(plotme,2);


%         subplot(2,2,3)
    f = figure;
    f.Position = [717   407   371   538];
    imagesc(t,trials,plotme(:,sortix)')
    colormap(linspecer)
    xlim(xlims)
    colorbar;
    title('me last lick')

    f = figure;
    f.Position = [717   407   371   538]; 
%     subplot(2,2,4)
    imagesc(t,trials,potent(:,sortix)')
    colormap(linspecer)
    xlim(xlims)
    colorbar;
    title('potent last lick')

    %     subplot(2,3,6)
    %     imagesc(t,trials,null(:,sortix)')
    %     colormap(linspecer)
    %     xlim(xlims)


end

%% quantify correlations b/w me and potent space for diff alignments

clear firstLickCC
clear lastLickCC

close all
for sessix = 1:numel(meta)-1

    temprez = rez(sessix);
    tempme = me(sessix);

    %     trix = 1:obj(sessix).bp.Ntrials;
    trix = find(~obj(sessix).bp.no);

    potent = temprez.N_potent;
    potent = potent(:,trix,:);
    potent = sum(potent.^2,3);

    plotme = tempme.data(:,trix);

    % sort
    ixs = 180:245;
    [~,sortix] = sort(mean(plotme(ixs,:),1),'descend');

    xlims = [-0.5 0.5];

    t = obj(sessix).time;
    [~,ix1] = min(abs(t - xlims(1)));
    [~,ix2] = min(abs(t - xlims(2)));

    plotme = plotme(ix1:ix2,:);
    potent = potent(ix1:ix2,:);
    cc = corrcoef(plotme(:),potent(:));
    firstLickCC(sessix) = cc(1,2);

    % last lick
    temprez = rez_lastLick(sessix);
    tempme = me_lastLick(sessix);

    potent = temprez.N_potent;
    potent = potent(:,trix,:);
    potent = sum(potent.^2,3);

    plotme = tempme.data(:,trix);

    t = obj_lastLick(sessix).time;
    [~,ix1] = min(abs(t - xlims(1)));
    [~,ix2] = min(abs(t - xlims(2)));

    plotme = plotme(ix1:ix2,:);
    potent = potent(ix1:ix2,:);
    cc = corrcoef(plotme(:),potent(:));
    lastLickCC(sessix) = cc(1,2);

end

%%

figure;
hold on;
scatter(lastLickCC,firstLickCC,100,'MarkerFaceColor',[0 0 0] ./ 255,'MarkerEdgeColor','w')
xlim([0 0.45])
ylim([0 0.45])

ax = gca;
plot(ax.XLim,ax.YLim,'k--','LineWidth',1)
ax.FontSize = 14;

xlabel('$R^2 (ME,Potent_{lastLick})$', 'Interpreter','latex','FontWeight','bold','FontSize',16)
ylabel('$R^2 (ME,Potent_{firstLick})$', 'Interpreter','latex','FontWeight','bold','FontSize',16);



%% quantify correlations b/w me and potent space for diff alignments

clear firstLickCC
clear lastLickCC
clear sessme_firstLick
clear sesspotent_firstLick
clear sessme_lastLick
clear sesspotent_lastLick

close all
for sessix = 1:numel(meta)-1

    temprez = rez(sessix);
    tempme = me(sessix);

    %     trix = 1:obj(sessix).bp.Ntrials;
    trix = find(obj(sessix).bp.hit);

    potent = temprez.N_potent;
    potent = potent(:,trix,:);
    potent = sum(potent.^2,3);

    plotme = tempme.data(:,trix);



    xlims = [-2 2];

    t = obj(sessix).time;
    [~,ix1] = min(abs(t - xlims(1)));
    [~,ix2] = min(abs(t - xlims(2)));

    plotme = plotme(ix1:ix2,:);
    potent = potent(ix1:ix2,:);

    binSize = 15; % ms
    dt = floor(binSize / (params(sessix).dt*1000)); % samples
    tm = obj(sessix).time(1:dt:numel(obj(sessix).time));
    numT = numel(tm);

    ix = 1;
    for i = 1:dt:(size(potent,1)-dt) % each timepoint
        ixs = i:(i+dt-1);
%         cc = corrcoef(plotme(ixs,:),potent(ixs,:));
%         firstLickCC(i,sessix) = cc(1,2);
        cc = corr(plotme(ixs,:)',potent(ixs,:)','type','Spearman');
        firstLickCC(i,sessix) = cc;
    end
    sessme_firstLick(:,sessix) = mean(plotme,2);
    sesspotent_firstLick(:,sessix) = mean(potent,2);

    % last lick
    temprez = rez_lastLick(sessix);
    tempme = me_lastLick(sessix);

    potent = temprez.N_potent;
    potent = potent(:,trix,:);
    potent = sum(potent.^2,3);

    plotme = tempme.data(:,trix);

    t = obj_lastLick(sessix).time;
    [~,ix1] = min(abs(t - xlims(1)));
    [~,ix2] = min(abs(t - xlims(2)));

    plotme = plotme(ix1:ix2,:);
    potent = potent(ix1:ix2,:);

    dt = floor(binSize / (params_lastLick(sessix).dt*1000)); % samples
    tm = obj_lastLick(sessix).time(1:dt:numel(obj_lastLick(sessix).time));
    numT = numel(tm);

    ix = 1;
    for i = 1:dt:(size(potent,1)-dt) % each timepoint
        ixs = i:(i+dt-1);
%         cc = corrcoef(plotme(ixs,:),potent(ixs,:));
%         lastLickCC(i,sessix) = cc(1,2);
        cc = corr(plotme(ixs,:)',potent(ixs,:)','type','Spearman');
        lastLickCC(i,sessix) = cc;
    end
    sessme_lastLick(:,sessix) = mean(plotme,2);
    sesspotent_lastLick(:,sessix) = mean(potent,2);

end

% sesspotent_lastLick = normalize(sesspotent_lastLick,'range',[0 1]);
% sessme_lastLick = normalize(sessme_lastLick,'range',[0 1]);
%

%%%%

cols = getColors;

alph = 0.2;

figure;
ax = gca;
hold on;
% shadedErrorBar(linspace(xlims(1),xlims(2),size(firstLickCC,1)),mean(firstLickCC,2),std(firstLickCC,[],2)./sqrt(numel(meta)),{'Color',[0 0 0],'LineWidth',2},0.1,ax)
shadedErrorBar(linspace(xlims(1),xlims(2),size(lastLickCC,1)),nanmedian(lastLickCC,2),nanstd(lastLickCC,[],2)./sqrt(numel(meta)),{'Color',cols.potent,'LineWidth',2},alph,ax)

%%%%
figure;
ax = gca;
hold on
yyaxis(ax,'left')
shadedErrorBar(linspace(xlims(1),xlims(2),size(sesspotent_lastLick,1)),nanmean(sesspotent_lastLick,2),nanstd(sesspotent_lastLick,[],2)./sqrt(numel(meta)),{'Color',cols.potent,'LineWidth',2,'LineStyle','-'},alph,ax)
ax.YColor = [0.1500    0.1500    0.1500]./255;


yyaxis(ax,'right');
shadedErrorBar(linspace(xlims(1),xlims(2),size(sessme_lastLick,1)),nanmean(sessme_lastLick,2),nanstd(sessme_lastLick,[],2)./sqrt(numel(meta)),{'Color',[0 0 0],'LineWidth',2,'LineStyle','-'},alph,ax)
yline(10,'k:')
ax.YColor = [0.1500    0.1500    0.1500]./255;
ylim([5 ax.YLim(2)])

xlabel('Time (s) from last lick')

% %%%%
%
% figure;
% ax = gca;
% hold on
% yyaxis(ax,'left')
% plot(linspace(xlims(1),xlims(2),size(sesspotent_lastLick,1)),sesspotent_lastLick,'Color',cols.potent,'LineStyle','-','Marker','none')
% 
% yyaxis(ax,'right');
% plot(linspace(xlims(1),xlims(2),size(sessme_lastLick,1)),sessme_lastLick,'k','LineStyle','-','Marker','none')
% yyaxis(ax,'right')
% yline(10,'k:')




