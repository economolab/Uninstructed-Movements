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
params.tmax = 2.5;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 51;

% cluster qualities to use
% params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality
params.quality = {'all'}; 

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

% Remove unwanted sessions

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
    

    inite1 = 0;
    inite2 = 0.4;
    dt = 0.4;
    
    for tt = 1:10
        e1 = inite1 + (tt-1)*dt;
        e2 = inite2 + (tt-1)*dt;
        
        e1  = mode(ev.sample) + e1 - mode(ev.(params.alignEvent));
        e2  = mode(ev.sample) + e2 - mode(ev.(params.alignEvent));
        
        times.(['mode' num2str(tt)]) = rez(sessix).time>e1 & rez(sessix).time<e2;
        rez(sessix).(['cd' num2str(tt) '_mode']) = calcCD(rez(sessix),times.(['mode' num2str(tt)]));
        
    end
    
  
    
    % orthogonalize
    
    [fns,~] = patternMatchCellArray(fieldnames(rez(sessix)),{'mode'},'all');
    modes = zeros(numel(rez(sessix).cd1_mode),numel(fns));
    for i = 1:numel(fns)
        modes(:,i) = rez(sessix).(fns{i});
    end
    
    orthModes = gschmidt(modes);
    
    for i = 1:numel(fns)
        rez(sessix).(fns{i}) = orthModes(:,i);
    end
    
    
    % projections and normalize
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
    
    % variance explained
    
    for i = 1:numel(fns)
        psth = rez(sessix).psth;
        datacov = cov([psth(:,:,1) ; psth(:,:,2)]);
        datacov(isnan(datacov)) = 0;
        eigsum = sum(eig(datacov));
        rez(sessix).varexp.(fns{i}(1:end-5)) = var_proj(rez(sessix).(fns{i}), datacov, eigsum);
    end

    
    
end

% concatenate latents from all sessions, find mean and stderror across sessions

fns = patternMatchCellArray(fieldnames(rez(1)),{'latent'},'all');

for fnix = 1:numel(fns)
    fn = fns{fnix};
    % preallocate
    cd.(fn){1} = [];
    cd.(fn){2} = [];
    % concatenate
    for i = 1:numel(rez)
        for j = 1:2
            cd.(fn){j} = [cd.(fn){j} rez(i).(fn)(:,j)];
        end
    end
    % mean and stderror across session (just for plotting)[p
    cd.mean.(fn) = [nanmean(cd.(fn){1},2) nanmean(cd.(fn){2},2)];
    cd.stderr.(fn) = [nanstd(cd.(fn){1},[],2) nanstd(cd.(fn){2},[],2)] ./ numel(rez); % std error
end

%% plot projections onto coding directions 
close all
clrs = getColors();
lw = 6;
alph = 0.5;

sample = mode(rez(1).ev.sample - rez(1).ev.(params.alignEvent));
delay = mode(rez(1).ev.delay - rez(1).ev.(params.alignEvent));

fns = fieldnames(cd.mean);

sav = 0;
f = figure; 
for i = 1:numel(fns)
    ax = nexttile; hold on
    tempmean = cd.mean.(fns{i});
    temperror = cd.stderr.(fns{i});
    shadedErrorBar(rez(1).time,tempmean(:,1),temperror(:,1),{'Color',clrs.rhit,'LineWidth',lw},alph, ax)
    shadedErrorBar(rez(1).time,tempmean(:,2),temperror(:,2),{'Color',clrs.lhit,'LineWidth',lw},alph, ax)
    
    xlim([rez(1).time(15);rez(1).time(end)])
    ylims = [min(min(tempmean))-5, max(max(tempmean))+5];
    ylim(ylims);
    
    title(['$ CD_{' lower(fns{i}(3:end-7)) '}$'],'Interpreter','latex')
%     xlabel('Time (s) from go cue')
%     ylabel('Activity (a.u.)')
    ax = gca;
    ax.FontSize = 12;
    ax.XTick = [];
    ax.YTick = [];
    
    xline(sample,'k--','LineWidth',2)
    xline(delay,'k--','LineWidth',2)
    xline(0,'k--','LineWidth',2)
    
    curmodename = fns{i};
    timefieldname = ['mode' lower(curmodename(3:end-7))];
    shadetimes = objs{1}.time(times.(timefieldname));
    x = [shadetimes(1)  shadetimes(end) shadetimes(end) shadetimes(1)];
    y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
    fl = fill(x,y,'r','FaceColor',[93, 121, 148]./255);
    fl.FaceAlpha = 0.3;
    fl.EdgeColor = 'none';
    
    ylim(ylims);
    
    
%     if sav
%         pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig1/figs/cd';
%         fn = [fns{i} '_anmList1_sessionList1_w_excludedsessions_sm_' num2str(params.smooth)];
%         mysavefig(f(i),pth,fn);
%     end

    hold off
    
end

%% orthogonality of modes
% close all
% dotProductModes(rez,orthModes,'Orthogonality of Coding Directions')
% f=gcf;
% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig1/figs/cd';
% fn = ['orthogonality_anmList1_sessionList1_w_excludedsessions_sm_' num2str(params.smooth)];
% mysavefig(f,pth,fn);


%% Helper Functions


function cd = calcCD(rez,times)
    tempdat = rez.psth(:,:,[1,2]);
    mu = squeeze(mean(tempdat(times,:,:),1));
    sd = squeeze(std(tempdat(times,:,:),[],1));
    cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
    cd(isnan(cd)) = 0;
    cd = cd./sum(abs(cd)); % (ncells,1)
end


