clear,clc,close all

addpath(genpath(pwd))

% finds cd early, late, go as defined in economo 2018

%% SET RUN PARAMS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to 

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'R&hit&~stim.enable&~autowater&~early'};         % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};         % left hits, no stim, aw off
params.condition(end+1) = {'R&hit&~stim.enable&autowater&~early'};          % right hits, no stim, aw on
params.condition(end+1) = {'L&hit&~stim.enable&autowater&~early'};          % left hits, no stim, aw on
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};        % error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};        % error left, no stim, aw off
params.condition(end+1) = {'~hit&~miss&~stim.enable&~autowater&~early'};    % ignore
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};           % hit 2afc
params.condition(end+1) = {'hit&~stim.enable&autowater&~early'};            % hit aw


% set conditions used for finding activity modes
aw = '2'; % 1-on, 2-off
stim = '0'; % 0-off
params.modecondition(1) = {['R&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};  %rhit, aw off 
params.modecondition(2) = {['L&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};  %lhit, aw off 
params.modecondition(3) = {['R&miss&autowater.nums==' aw '&stim.num==' stim '&~early']}; %rmiss, aw off 
params.modecondition(4) = {['L&miss&autowater.nums==' aw '&stim.num==' stim '&~early']}; %lmiss, aw off 
params.modecondition(5) = {['hit&autowater.nums==' aw '&stim.num==' stim '&~early']};    % hit, aw off 

params.tmin = -2.5;
params.tmax = 4.5;
params.dt = 1/100;

% smooth with causal gaussian kernel
params.smooth = 51;

% cluster qualities to use
% params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality
params.quality = {'Poor,','Fair','Good','Great','Excellent','single'}; 

%% SET METADATA

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];
% meta = loadJEB4_ALMVideo(meta,datapth); % done
% meta = loadJEB5_ALMVideo(meta,datapth); % done
meta = loadJEB6_ALMVideo(meta,datapth); % done
meta = loadJEB7_ALMVideo(meta,datapth); % done
meta = loadEKH1_ALMVideo(meta,datapth); % done
meta = loadEKH3_ALMVideo(meta,datapth); % done
meta = loadJGR2_ALMVideo(meta,datapth); % done
meta = loadJGR3_ALMVideo(meta,datapth); % done


params.probe = [meta.probe]; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written

%% LOAD AND PROCESS DATA

objs = loadObjs(meta);


for metaix = 1:numel(meta)
    obj = objs{metaix};
    disp('______________________________________________________')
    disp(['Processing data for session ' [meta(metaix).anm '_' meta(metaix).date]])
    disp(' ')
    [sessparams{metaix},sessobj{metaix}] = processData(obj,params,params.probe(metaix));
end

% clean up sessparams and sessobj
for metaix = 1:numel(meta)
    params.trialid{metaix} = sessparams{metaix}.trialid;
    params.cluid{metaix} = sessparams{metaix}.cluid{params.probe(metaix)};
    
    objs{metaix} = sessobj{metaix};
    objs{metaix}.psth = objs{metaix}.psth{params.probe(metaix)};
    objs{metaix}.trialdat = objs{metaix}.trialdat{params.probe(metaix)};
    objs{metaix}.presampleFR = objs{metaix}.presampleFR{params.probe(metaix)};
    objs{metaix}.presampleSigma = objs{metaix}.presampleSigma{params.probe(metaix)};
end

disp(' ')
disp('DATA LOADED AND PROCESSED')
disp(' ')

%% Remove unwanted sessions

% remove sessions if:
% 1) less than 40 trials of rhit and lhit each (same as
% 2) atleast 10 clusters of quality that is not noise
use = false(size(objs));
for i = 1:numel(use)
    check(1) = numel(params.trialid{i}{1}) > 40;
    check(2) = numel(params.trialid{i}{2}) > 40;
    check(3) = numel(params.cluid{i}) >= 10;
    if all(check)
        use(i) = true;
    end
end

meta = meta(use);
objs = objs(use);
params.probe = params.probe(use);
params.trialid = params.trialid(use);
params.cluid = params.cluid(use);

%% CODING DIMENSIONS

clearvars -except objs meta params 

for sessix = 1:numel(objs)
    
    
    ev.sample = objs{sessix}.bp.ev.sample;
    ev.delay = objs{sessix}.bp.ev.delay;
    ev.goCue = objs{sessix}.bp.ev.goCue;
    ev.(params.alignEvent) = objs{sessix}.bp.ev.(params.alignEvent);
        
    rez(sessix).time = objs{sessix}.time;
    rez(sessix).psth = objs{sessix}.psth;
    rez(sessix).psth = standardizePSTH(objs{sessix});
    rez(sessix).condition = params.condition;
    rez(sessix).alignEvent = params.alignEvent;
    rez(sessix).ev = ev;
    
    %% cd early mode 
    
    cond{1} = params.modecondition{1};
    cond{2} = params.modecondition{2};
    epoch = 'sample';
    
    e1 = mode(ev.delay) - 0.5 - mode(ev.(params.alignEvent));
    e2 = mode(ev.delay) - 0.1 - mode(ev.(params.alignEvent));
%     e1 = mode(ev.sample) + 0.4 - mode(ev.(params.alignEvent));
%     e2 = mode(ev.sample) + 0.8 - mode(ev.(params.alignEvent));
    
    times.early = rez(sessix).time>e1 & rez(sessix).time<e2;
    tempdat = rez(sessix).psth(:,:,[1,2]);
    mu = squeeze(mean(tempdat(times.early,:,:),1));
    sd = squeeze(std(tempdat(times.early,:,:),[],1));
    cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
    cd(isnan(cd)) = 0;
    cd = cd./sum(abs(cd)); % (ncells,1)
    rez(sessix).cdEarly_mode = cd;
    
    clear cond
    
    %% cd late mode
    
    
    cond{1} = params.modecondition{1};
    cond{2} = params.modecondition{2};
    epoch = 'latedelay';
    
    
%     times.late = rez(sessix).time>-0.5 & rez(sessix).time<-0.1;
    e1 = mode(ev.goCue) - 0.5 - mode(ev.(params.alignEvent));
    e2 = mode(ev.goCue) - 0.1 - mode(ev.(params.alignEvent));
    
    times.late = rez(sessix).time>e1 & rez(sessix).time<e2;
    tempdat = rez(sessix).psth(:,:,[1,2]);
    mu = squeeze(mean(tempdat(times.late,:,:),1));
    sd = squeeze(std(tempdat(times.late,:,:),[],1));
    cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
    cd(isnan(cd)) = 0;
    cd = cd./sum(abs(cd)); % (ncells,1)
    rez(sessix).cdLate_mode = cd;
    
%     rez(sessix).cdLate_mode = calcCD(objs{sessix},params,cond,epoch,rez(sessix).alignEvent,sessix);
    
    clear cond
    
    %% cd go mode
    cond{1} = params.modecondition{1};
    cond{2} = params.modecondition{2};
    epoch = 'go';
    
    e1 = mode(ev.goCue) + 0.02 - mode(ev.(params.alignEvent));
    e2 = mode(ev.goCue) + 0.42 - mode(ev.(params.alignEvent));
    
    times.go = rez(sessix).time>e1 & rez(sessix).time<e2;
    
%     times.go = rez(sessix).time>0.02 & rez(sessix).time<0.42;
    tempdat = rez(sessix).psth(:,:,[1,2]);
    mu = squeeze(mean(tempdat(times.go,:,:),1));
    sd = squeeze(std(tempdat(times.go,:,:),[],1));
    cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
    cd(isnan(cd)) = 0;
    cd = cd./sum(abs(cd)); % (ncells,1)
    rez(sessix).cdGo_mode = cd;
    
%     rez(sessix).cdGo_mode = calcCD(objs{sessix},params,cond,epoch,rez(sessix).alignEvent,sessix);
    
    clear cond
    
     %% cd after mode
    cond{1} = params.modecondition{1};
    cond{2} = params.modecondition{2};
    epoch = 'go';
    
    e1 = mode(ev.goCue) + 1.5 - mode(ev.(params.alignEvent));
    e2 = mode(ev.goCue) + 2.0 - mode(ev.(params.alignEvent));
    
    times.after = rez(sessix).time>e1 & rez(sessix).time<e2;
    
%     times.go = rez(sessix).time>0.02 & rez(sessix).time<0.42;
    tempdat = rez(sessix).psth(:,:,[1,2]);
    mu = squeeze(mean(tempdat(times.go,:,:),1));
    sd = squeeze(std(tempdat(times.go,:,:),[],1));
    cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
    cd(isnan(cd)) = 0;
    cd = cd./sum(abs(cd)); % (ncells,1)
    rez(sessix).cdAfter_mode = cd;
    
%     rez(sessix).cdGo_mode = calcCD(objs{sessix},params,cond,epoch,rez(sessix).alignEvent,sessix);
    
    clear cond
    
    %% orthogonalize
    
    [fns,~] = patternMatchCellArray(fieldnames(rez(sessix)),{'mode'},'all');
    modes = zeros(numel(rez(sessix).cdLate_mode),numel(fns));
    for i = 1:numel(fns)
        modes(:,i) = rez(sessix).(fns{i});
    end
    
    orthModes = gschmidt(modes);
    
    for i = 1:numel(fns) 
        rez(sessix).(fns{i}) = orthModes(:,i);
    end
    
    % lastly, othogonalize early mode to go mode
%     tempmodes = [rez(sessix).cdGo_mode rez(sessix).cdEarly_mode];
%     orthModes = gschmidt(tempmodes);
%     rez(sessix).cdEarly_mode = orthModes(:,2);
    
    %% projections and normalize
    % when pooling trajectories across sessions from hidehikos ppn paper:
    % CD_late projections normalized by mean activity just before go cue
    % (-0.1<t<t_go)
    % CD_go projections normalized by mean activity after go cue
    % (t_go<t<0.4)
    
%     normTimes{1} = rez(sessix).time>e1 & rez(sessix).time<e2; % sample
%     normTimes{2} = rez(sessix).time>-0.4 & rez(sessix).time<0; % delay
%     normTimes{3} = rez(sessix).time>0 & rez(sessix).time<0.4; % go
    
    cond = [1 2];
    for i = 1:numel(fns)
        tempmode = rez(sessix).(fns{i});
        for j = 1:numel(cond)
            c = cond(j);
            
            tempdat = rez(sessix).psth(:,:,c)*rez(sessix).(fns{i});
            
%             normfactor = abs(nanmean(tempdat(normTimes{i})));
            normfactor = 1;
%             normfactor = max(tempdat);
            
            rez(sessix).([fns{i}(1:end-5) '_latent'])(:,j) = tempdat ./ normfactor;
        end
    end
    
    clear cond
    
    %% variance explained
    
    for i = 1:numel(fns)
        psth = rez(sessix).psth;
        datacov = cov([psth(:,:,1) ; psth(:,:,2)]);
        datacov(isnan(datacov)) = 0;
        eigsum = sum(eig(datacov));
        rez(sessix).varexp.(fns{i}(1:end-5)) = var_proj(rez(sessix).(fns{i}), datacov, eigsum);
    end

    
    
end

%% concatenate latents, find mean and stderror

cdEarly{1} = rez(1).cdEarly_latent(:,1);
cdEarly{2} = rez(1).cdEarly_latent(:,2);
cdLate{1} = rez(1).cdLate_latent(:,1);
cdLate{2} = rez(1).cdLate_latent(:,2);
cdGo{1} = rez(1).cdGo_latent(:,1);
cdGo{2} = rez(1).cdGo_latent(:,2);
cdAfter{1} = rez(1).cdGo_latent(:,1);
cdAfter{2} = rez(1).cdGo_latent(:,2);
for i = 2:numel(rez)
    for j = 1:2
        cdEarly{j} = cat(2,cdEarly{j},rez(i).cdEarly_latent(:,j));
        cdLate{j} = cat(2,cdLate{j},rez(i).cdLate_latent(:,j));
        cdGo{j} = cat(2,cdGo{j},rez(i).cdGo_latent(:,j));
        cdGolat{j} = cat(2,cdAfter{j},rez(i).cdAfter_latent(:,j));
    end
end

cdEarly_latent_mean = [nanmean(cdEarly{1},2) nanmean(cdEarly{2},2)];
cdLate_latent_mean = [nanmean(cdLate{1},2) nanmean(cdLate{2},2)];
cdGo_latent_mean = [nanmean(cdGo{1},2) nanmean(cdGo{2},2)];
cdAfter_latent_mean = [nanmean(cdAfter{1},2) nanmean(cdAfter{2},2)];

cdEarly_latent_error = [nanstd(cdEarly{1},[],2) nanstd(cdEarly{2},[],2)] ./ numel(rez); % std error
cdLate_latent_error = [nanstd(cdLate{1},[],2) nanstd(cdLate{2},[],2)] ./ numel(rez); % std error
cdGo_latent_error = [nanstd(cdGo{1},[],2) nanstd(cdGo{2},[],2)] ./ numel(rez); % std error
cdAfter_latent_error = [nanstd(cdAfter{1},[],2) nanstd(cdAfter{2},[],2)] ./ numel(rez); % std error

%%
close all
clrs = getColors();
lw = 6;
alph = 0.5;

sample = mode(rez(1).ev.sample - rez(1).ev.(params.alignEvent));
delay = mode(rez(1).ev.delay - rez(1).ev.(params.alignEvent));

sav = 0;
for i = 1:numel(fns)
    f(i) = figure; ax = axes(f(i)); hold on
    tempmean = eval([fns{i}(1:end-5) '_latent_mean']);
    temperror = eval([fns{i}(1:end-5) '_latent_error']);
    shadedErrorBar(rez(1).time,tempmean(:,1),temperror(:,1),{'Color',clrs.rhit,'LineWidth',lw},alph, ax)
    shadedErrorBar(rez(1).time,tempmean(:,2),temperror(:,2),{'Color',clrs.lhit,'LineWidth',lw},alph, ax)
    
    xlim([rez(1).time(15);rez(1).time(end)])
    ylims = [min(min(tempmean))-5, max(max(tempmean))+5];
    ylim(ylims);
    
    title(fns{i},'Interpreter','none')
    xlabel('Time (s) from go cue')
    ylabel('Activity (a.u.)')
    ax = gca;
    ax.FontSize = 40;
    
    xline(sample,'k--','LineWidth',2)
    xline(delay,'k--','LineWidth',2)
    xline(0,'k--','LineWidth',2)
    
    curmodename = fns{i};
    timefns = fieldnames(times);
    mask = strfind(timefns,lower(curmodename(3:end-5)));
    ix = cellfun(@(x) isempty(x),mask,'UniformOutput',false);
    ix = ~cell2mat(ix);
    shadetimes = objs{1}.time(times.(timefns{ix}));
    x = [shadetimes(1)  shadetimes(end) shadetimes(end) shadetimes(1)];
    y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
    fl = fill(x,y,'r','FaceColor',[93, 121, 148]./255);
    fl.FaceAlpha = 0.3;
    fl.EdgeColor = 'none';
    
    ylim(ylims);
    
    
    if sav
        pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig1/figs/cd';
        fn = [fns{i} '_anmList1_sessionList1_w_excludedsessions_sm_' num2str(params.smooth)];
        mysavefig(f(i),pth,fn);
    end
    
end

%% orthogonality of modes
% close all
% dotProductModes(rez,orthModes,'Orthogonality of Coding Directions')
% f=gcf;
% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig1/figs/cd';
% fn = ['orthogonality_anmList1_sessionList1_w_excludedsessions_sm_' num2str(params.smooth)];
% mysavefig(f,pth,fn);


%% selectivity
clear selectivity

sumsqselectivity.total = zeros(numel(rez(1).time),numel(rez));
sumsqselectivity.early = zeros(numel(rez(1).time),numel(rez));
sumsqselectivity.late = zeros(numel(rez(1).time),numel(rez));
sumsqselectivity.go = zeros(numel(rez(1).time),numel(rez));
for i = 1:numel(rez)
    selectivity.total = rez(i).psth(:,:,1) - rez(i).psth(:,:,2);
    selectivity.early = rez(i).cdEarly_latent(:,1) - rez(i).cdEarly_latent(:,2);
    selectivity.late  = rez(i).cdLate_latent(:,1) - rez(i).cdLate_latent(:,2);
    selectivity.go    = rez(i).cdGo_latent(:,1) - rez(i).cdGo_latent(:,2);
    selfns = fieldnames(selectivity);
    for j = 1:numel(selfns)
        selectivity.(selfns{j})(isnan(selectivity.(selfns{j}))) = 0;
        sumsqselectivity.(selfns{j})(:,i) = sum(selectivity.(selfns{j}).^2,2);
    end
end

% selectivity explained
% calculate variance of selectivity matrix within each of the modes
clear selectivity
% only calculate variance explained between start of cdEarly and end of
% cdGo- all other selectivity is variance that we are not attempting to
% explain with the three CDs
ix1 = find(times.early,1,'first');
ix2 = find(times.go,1,'last');
for i = 1:numel(rez)
    selectivity(i).total = rez(i).psth(ix1:ix2,:,1) - rez(i).psth(ix1:ix2,:,2);
    selectivity(i).early = rez(i).cdEarly_latent(ix1:ix2,1) - rez(i).cdEarly_latent(ix1:ix2,2);
    selectivity(i).late  = rez(i).cdLate_latent(ix1:ix2,1) - rez(i).cdLate_latent(ix1:ix2,2);
    selectivity(i).go    = rez(i).cdGo_latent(ix1:ix2,1) - rez(i).cdGo_latent(ix1:ix2,2);
    
    covtot = cov(selectivity(i).total);
    eigs = sort(eig(covtot),'descend');
    eigsum = sum(eigs(1:3));
    selexp.early(i) = var_proj(rez(i).cdEarly_mode,covtot,eigsum);
    selexp.late(i) = var_proj(rez(i).cdLate_mode,covtot,eigsum);
    selexp.go(i) = var_proj(rez(i).cdGo_mode,covtot,eigsum);
    
%     varianceInTotalSelectivity = trace(cov(selectivity(i).total));
%     selexp.early(i) = trace(cov(selectivity(i).early)) ./ varianceInTotalSelectivity;
%     selexp.late(i) = trace(cov(selectivity(i).late)) ./ varianceInTotalSelectivity;
%     selexp.go(i) = trace(cov(selectivity(i).go)) ./ varianceInTotalSelectivity;
    selexp.sum(i) = selexp.early(i) + selexp.late(i)  + selexp.go(i);
end
selectivityExplained = [selexp.early' selexp.late' selexp.go' selexp.sum'];

%%
close all
sample = mode(rez(1).ev.sample) - mode(rez(1).ev.(params.alignEvent));
delay  = mode(rez(1).ev.delay) - mode(rez(1).ev.(params.alignEvent));

clrs = {'k','b','g','m','c'};

lw = 4;
alph = 0.5;
f = figure; ax = axes(f); hold on;
summean = zeros(numel(rez(1).time),1);
sumstd = zeros(numel(rez(1).time),1);
for i = 1:numel(selfns)
    temp = sumsqselectivity.(selfns{i});
    tempmean = nanmean(temp,2);
    tempstd = nanstd(temp,[],2);
    shadedErrorBar(rez(1).time,tempmean,...
                   tempstd./(numel(rez)),...
                   {'Color',clrs{i},'LineWidth',lw},alph,ax);
   if ~strcmpi(selfns{i},'total')
       summean = summean + tempmean;
       sumstd = sumstd + tempstd;
   end
end
shadedErrorBar(rez(1).time,summean,...
                   sumstd./(numel(rez)),...
                   {'Color',clrs{end},'LineWidth',lw},alph,ax);

xline(sample,'k--','LineWidth',2)
xline(delay,'k--','LineWidth',2)
xline(0,'k--','LineWidth',2)

xlabel('Time (s) from go cue')
ylabel('Squared Sum Selectivity')
% legend('Total selectivity','early','late','go','early + late + go')
xlim([rez(1).time(1)+0.2,rez(1).time(end)])
ax = gca;
ax.FontSize = 20;
 
% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig1/figs/selectivity';
% fn = 'sqsumselectivity_hits_anmList1_sessionList1_psthsm_51_dt_0_005_lowFR_1';
% mysavefig(f,pth,fn);

% variance in selectivity explained
violincols = [70, 70, 235; 70, 235, 81; 224, 70, 235; 70, 235, 235] ./ 255;

varfns = fieldnames(selexp);
f = figure; ax = axes(f);
vs = violinplot(selectivityExplained,{'early','late','go','sum'},...
    'EdgeColor',[1 1 1], 'ViolinAlpha',{0.35,1}, 'ViolinColor',violincols);
ylabel('Selectivity Explained')
ylim([0,1])
ax = gca;
ax.FontSize = 25;

pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig1/figs/selectivity';
fn = 'selexp';
mysavefig(f,pth,fn);


%% correlation pop selectivity vector

selpsth = objs{1}.psth;
for sessix = 2:numel(objs)
    obj = objs{sessix};    
    selpsth = cat(2,selpsth,obj.psth);
end
selpsth(isnan(psth)) = 0;
selpsth(isinf(psth)) = 0;
sm = 31;
psth_selectivity = mySmooth(selpsth(:,:,1) - selpsth(:,:,2),sm);

corr_matrix_selectivity = zeros(size(psth_selectivity,1),size(psth_selectivity,1));

for i = 1:size(corr_matrix_selectivity,1)
    for j = 1:size(corr_matrix_selectivity,1)
        temp = corrcoef(psth_selectivity(i,:),psth_selectivity(j,:));
        corr_matrix_selectivity(i,j) = temp(1,2);
    end
end

%% plot selectivity correlation matrix
close all
f = figure; hold on;
imagesc(obj.time,obj.time,corr_matrix_selectivity);
colorbar; caxis([0 max(max(corr_matrix_selectivity))]);

lw = 4;
ls = '--';
col = [88, 245, 112] ./ 255;
xline(sample,ls,'Color',col,'LineWidth',lw); yline(sample,ls,'Color',col,'LineWidth',lw)
xline(delay,ls,'Color',col,'LineWidth',lw); yline(delay,ls,'Color',col,'LineWidth',lw)
xline(0,ls,'Color',col,'LineWidth',lw); yline(0,ls,'Color',col,'LineWidth',lw)

xlim([rez(1).time(1)+0.2,rez(1).time(end)]);
ylim([rez(1).time(1)+0.2,rez(1).time(end)])
xlabel('Time (s) from go cue')
ylabel('Time (s) from go cue')
ax = gca;
hold off
colormap(hot)
a = colorbar;
a.Label.String = 'Correlation of population selectivity vector';
ax.FontSize = 25;

% 
% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig1/figs/selectivity';
% fn = 'popselectivitycorr_hits_anmList1_sessionList1_psthsm_51_dt_0_005_lowFR_1';
% mysavefig(f,pth,fn);



