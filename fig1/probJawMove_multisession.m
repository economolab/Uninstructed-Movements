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
params.smooth = 15;

% cluster qualities to use
% params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality
params.quality = {'Fair','Good','Great','Excellent','single'}; 

%% SET METADATA

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];
meta = loadJEB4_ALMVideo(meta,datapth); % done
meta = loadJEB5_ALMVideo(meta,datapth); % done
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

%% COMBINE DATA ACROSS SESSIONS

% TODO: Don't normalize here, and see how this changes cds and selectivity.
% currently, total selectivity looks weird with normalized data (no
% ramping)

psth = normalizePSTH(objs{1});
% psth = objs{1}.psth;
ev.sample = objs{1}.bp.ev.sample;
ev.delay = objs{1}.bp.ev.delay;
ev.goCue = objs{1}.bp.ev.goCue;
ev.(params.alignEvent) = objs{1}.bp.ev.(params.alignEvent);
for sessix = 2:numel(objs)
    obj = objs{sessix};
    
    temppsth = normalizePSTH(obj);
%     temppsth = obj.psth;
    
    psth = cat(2,psth,temppsth);
    
    ev.sample = [ev.sample ; obj.bp.ev.sample];
    ev.delay = [ev.delay ; obj.bp.ev.delay];
    ev.goCue = [ev.goCue ; obj.bp.ev.goCue];
    ev.(params.alignEvent) = [ev.(params.alignEvent) ; obj.bp.ev.(params.alignEvent)];
    
end

psth(isnan(psth)) = 0;
psth(isinf(psth)) = 0;

%% probability of jaw movement

thresh = 0.2;
smth = 21;
view = 1; % side cam

conditions = [1,2];
jawSelectivity = nan(numel(objs{1}.time),numel(objs));
% figure;
for i = 1:numel(objs)
    obj = objs{i};
    jawix = findDLCFeatIndex(objs{i}.traj,view,{'jaw'});
    jawprob = jawProbSessionAvg(obj,view,params.trialid{i},jawix,thresh,conditions,params.dt,params.alignEvent);

    jawSelectivity(:,i) = jawprob{1} - jawprob{2};
%     plot(objs{1}.time,jawSelectivity(:,i));
%     hold on
%     pause
end

jawSelectivity = mySmooth(jawSelectivity,smth);
% figure; plot(jawSelectivity)

%% correlation

corr_matrix_selectivity = zeros(size(jawSelectivity,1),size(jawSelectivity,1));

for i = 1:size(corr_matrix_selectivity,1)
    for j = 1:size(corr_matrix_selectivity,1)
        temp = corrcoef(jawSelectivity(i,:),jawSelectivity(j,:));
        corr_matrix_selectivity(i,j) = temp(1,2);
    end
end

%% plot selectivity correlation matrix

sample = mode(ev.sample) - mode(ev.(params.alignEvent));
delay  = mode(ev.delay) - mode(ev.(params.alignEvent));


figure; hold on;
imagesc(objs{1}.time,objs{1}.time,corr_matrix_selectivity);
colorbar; caxis([0 max(max(corr_matrix_selectivity))]);

lw = 4;
ls = 'g--';
xline(sample,ls,'LineWidth',lw); yline(sample,ls,'LineWidth',lw)
xline(delay,ls,'LineWidth',lw); yline(delay,ls,'LineWidth',lw)
xline(0,ls,'LineWidth',lw); yline(0,ls,'LineWidth',lw)

xlim([objs{1}.time(1)+0.2,objs{1}.time(end)]);
ylim([objs{1}.time(1)+0.2,objs{1}.time(end)])
title('Correlation for selectivity of jaw move prob')
ax = gca;
ax.FontSize = 20;
hold off
colormap(hot)




%% Helper functions

function jawprob = jawProbSessionAvg(obj,view,trialid,jawix,thresh,conditions,dt,alignEv)

jawprob = cell(1,numel(conditions));
for i = 1:numel(conditions)
    cond = conditions(i);
    trialsToUse = trialid{cond};

    derivthresh = thresh;  %If the velocity of the jaw crosses this threshold, the jaw is considered to be moving
    edges = (obj.time(1):dt:obj.time(end)) + mode(obj.bp.ev.(alignEv));
    Ntrials = length(trialsToUse);

    good_ix = NaN(Ntrials, 1);

    traj = obj.traj{view};

    %Find the trials that have numerical values for frameTimes (aka no NaNs)
    for ii = 1:Ntrials
        ix = trialsToUse(ii);
        if ~any(isnan(traj(ix).frameTimes))
            good_ix(ii,1) = ix;       %If the trial is okay, note the trial number
        end
    end

    %Get trial numbers of usable trials
    use_trials = find(~isnan(good_ix))';
    trialsToUse = trialsToUse(use_trials);

    jaw = NaN(numel(edges), length(trialsToUse));

    for ii = 1:length(trialsToUse)
        q = trialsToUse(ii);
        
        ts = mySmooth(traj(q).ts(:, 2, jawix), 21);
        tsinterp = interp1(traj(q).frameTimes-0.5, ts, edges);   %Linear interpolation of jaw position to keep number of time points consistent across trials
        basederiv = median(diff(tsinterp),"omitnan");                   %Find the median jaw velocity (aka baseline)

        %Find when the difference between the jaw velocity and the
        %baseline jaw
        %velocity is above a given threshold (when is jaw moving?)
        jaw(2:end, ii) = abs(diff(tsinterp)-basederiv)>derivthresh;% | abs(tsinterp(2:end)-basepos)>posthresh;

    end
    
%     nojawmovix = ~ismember(trialsToUse, find(obj.earlyMoveix));
%     jawdat_noearly = mean(jaw(:, nojawmovix), 2, 'omitnan');
%     jawdat_noearly = medfilt1(jawdat_noearly, 10);
    
    
    jawdat = mean(jaw,2,'omitnan');         % Take the average probability of jaw movement across trials
    jawdat = medfilt1(jawdat,10);           % Smooth the avg probability of jaw movement 
    jawprob{i} = jawdat;                % Store the avg prob of jaw movement for the current condition
    
%     figure; imagesc(jaw)
%     figure; imagesc(jaw(:, nojawmovix))
end
end % jawProbSessionAvg













