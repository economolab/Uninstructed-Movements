clear,clc,close all

addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Data-Loading-Scripts'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements'));
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
params.smooth = 31;

% cluster qualities to use
% params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality
params.quality = {'excellent','great','good','multi','fair','poor'};

%% SET METADATA

meta = [];
meta = loadJEB6_ALMVideo(meta);
meta = loadJEB7_ALMVideo(meta);
meta = loadEKH1_ALMVideo(meta);
meta = loadEKH3_ALMVideo(meta);
meta = loadJGR2_ALMVideo(meta);
meta = loadJGR3_ALMVideo(meta);

params.probe = [meta.probe]; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written

%% LOAD AND PROCESS DATA

objs = loadObjs(meta);

for metaix = 1:numel(meta)
    obj = objs{metaix};
    disp('______________________________________________________')
    disp(['Processing data for session ' [meta(metaix).anm '_' meta(metaix).date]])             % Display progress
    disp(' ')
    [sessparams{metaix},sessobj{metaix}] = processData(obj,params,params.probe(metaix));        % Function for processing all of the data objs
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

    % cd before mode

    %     e1 = mode(ev.sample) - mode(ev.(params.alignEvent));
    %     e2 = mode(ev.sample) + 0.4 - mode(ev.(params.alignEvent));
    %
    %     times.before = rez(sessix).time>e1 & rez(sessix).time<e2;
    %     rez(sessix).cdBefore_mode = calcCD(rez(sessix),times.before);


    % cd early mode

    e1 = mode(ev.delay) - 0.5 - mode(ev.(params.alignEvent));
    e2 = mode(ev.delay) - 0.1 - mode(ev.(params.alignEvent));
    %     e1 = mode(ev.sample) + 0.4 - mode(ev.(params.alignEvent));
    %     e2 = mode(ev.sample) + 0.8 - mode(ev.(params.alignEvent));

    times.early = rez(sessix).time>e1 & rez(sessix).time<e2;
    rez(sessix).cdEarly_mode = calcCD(rez(sessix),times.early);


    % cd late mode

    e1 = mode(ev.goCue) - 0.5 - mode(ev.(params.alignEvent));
    e2 = mode(ev.goCue) - 0.1 - mode(ev.(params.alignEvent));

    times.late = rez(sessix).time>e1 & rez(sessix).time<e2;
    rez(sessix).cdLate_mode = calcCD(rez(sessix),times.late);


    % cd go mode

    e1 = mode(ev.goCue) + 0.02 - mode(ev.(params.alignEvent));
    e2 = mode(ev.goCue) + 0.42 - mode(ev.(params.alignEvent));

    times.go = rez(sessix).time>e1 & rez(sessix).time<e2;
    rez(sessix).cdGo_mode = calcCD(rez(sessix),times.go);

    % cd after mode

    %     e1 = mode(ev.goCue) + 1.5 - mode(ev.(params.alignEvent));
    %     e2 = mode(ev.goCue) + 2.0 - mode(ev.(params.alignEvent));
    %
    %     times.after = rez(sessix).time>e1 & rez(sessix).time<e2;
    %     rez(sessix).cdAfter_mode = calcCD(rez(sessix),times.after);

    % orthogonalize

    [fns,~] = patternMatchCellArray(fieldnames(rez(sessix)),{'mode'},'all');
    modes = zeros(numel(rez(sessix).cdLate_mode),numel(fns));
    for i = 1:numel(fns)
        modes(:,i) = rez(sessix).(fns{i});
    end

    orthModes = gschmidt(modes);

    for i = 1:numel(fns)
        rez(sessix).(fns{i}) = orthModes(:,i);
    end


    % projections and normalize (not normalizing currently since I
    % standardize PSTHs beforehand)
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
%% Adjust for older functions
for i = 1:numel(objs)
    meta(i).trialid = params.trialid{i};
end
%% Without removing jaw move
for i = 1:numel(objs)
    figure();
    obj = objs{i};
    met = meta(i);
    sesstitle = strcat(meta(i).anm, meta(i).date);

    go = mode(obj.bp.ev.goCue)-mode(obj.bp.ev.goCue);
    del = go-0.9;
    samp = del-1.3;

    conditions = {1,2};
    colors = {[0 0 1],[1 0 0]};

    subplot(4,1,2)
    plot(rez(i).time,rez(i).cdEarly_latent(:,1),'blue','LineWidth',2); hold on; plot(rez(i).time,rez(i).cdEarly_latent(:,2),'red','LineWidth',2);
    addtriallines(go,samp,del)
    legend('Right','Left')
    ylabel('a.u.')
    title('CDearly')
    xlim([-2.3 2.5])
    

    subplot(4,1,3)
    plot(rez(i).time,rez(i).cdLate_latent(:,1),'blue','LineWidth',2); hold on; plot(rez(i).time,rez(i).cdLate_latent(:,2),'red','LineWidth',2);
    addtriallines(go,samp,del)
    legend('Right','Left')
    ylabel('a.u.')
    title('CDlate')
    xlim([-2.3 2.5])

    subplot(4,1,4)
    plot(rez(i).time,rez(i).cdGo_latent(:,1),'blue','LineWidth',2); hold on; plot(rez(i).time,rez(i).cdGo_latent(:,2),'red','LineWidth',2);
    addtriallines(go,samp,del)
    legend('Right','Left')
    ylabel('a.u.')
    title('CDGo')
    xlim([-2.3 2.5])


    subplot(4,1,1)
    taxis = rez(i).time;
    
    jaw_by_cond = findJawVelocity(taxis, obj,conditions,met,'vel',params);    % (1 x conditions cell array)
    % Each cell: (time x trials in that condition)
    jawvel.right = mean(jaw_by_cond{1},2,'omitnan');
    jawvel.left = mean(jaw_by_cond{2},2,'omitnan');
    plot(taxis, jawvel.right,'blue','LineWidth',2); hold on;
    plot(taxis, jawvel.left,'red','LineWidth',2);
    
    ylabel('Jaw velocity')
    addtriallines(go,samp,del)
    legend('Right','Left','Location','best')
    title('Jaw')
    xlim([-2.3 2.5])

    sgtitle(sesstitle)
end
%%  Find jaw mode; orthogonalize other modes to it
rezNoJaw = rez;
for gg = 1:numel(objs)
    obj = objs{gg};
    met = meta(gg);
    sesstitle = strcat(meta(gg).anm, meta(gg).date);

    for c = 1:numel(conditions)
        trix = met.trialid{c};
        obj.trialpsth_cond{c} = obj.trialdat(:,:,trix);
    end

    % Find jaw kinematics
    view = 1; % side
    feat = 4; % jaw
    [~,kin.jawVel] = getFeatureKinematics(taxis,obj,conditions,met,view,feat);

    params.tix = 1:1000;       % time points to use when finding mode
    params.fcut = 50;          % smoothing cutoff frequency
    params.cond = 1:2;         % which conditions to use to find mode
    params.method = 'xcorr';   % 'xcorr' or 'regress' (basically the same)
    params.fa = false;          % if true, reduces neural dimensions to 10 with factor analysis
    
    % get modes based on single trial full neural data (or latents) and kinematic features 
    kinfns = fieldnames(kin);
    for i = 1:numel(kinfns)
        Y = kin.(kinfns{i}); % feature data to use to calculate mode
        [jawmode.(kinfns{i}), dat.(kinfns{i})] = findMode(obj, Y, params,rez);
    end
%     R = rezNoJaw(gg).psth(:,:,1)*jawmode.jawVel; R = R./ normfactor;
%     L = rezNoJaw(gg).psth(:,:,2)*jawmode.jawVel; L = L./ normfactor;
%     figure(); plot(rez(gg).time,R,'blue','LineWidth',2); hold on; plot(rez(gg).time,L,'red','LineWidth',2)
%     title(sesstitle)

    rezNoJaw(gg).cdjawVel_mode = jawmode.jawVel;

    %[fns,~] = patternMatchCellArray(fieldnames(rezNoJaw(sessix)),{'mode'},'all');
    fns = {'cdjawVel_mode';'cdEarly_mode';'cdLate_mode';'cdGo_mode'};
    modes = zeros(numel(rezNoJaw(gg).cdLate_mode),numel(fns));
    for i = 1:numel(fns)
        modes(:,i) = rezNoJaw(gg).(fns{i});
    end

    orthModes = gschmidt(modes);

    for i = 1:numel(fns)
        rezNoJaw(gg).(fns{i}) = orthModes(:,i);
    end

end
%% Project onto the jaw mode and CDs which are orthogonalized to jaw mode
for sessix = 1:numel(objs)
    cond = [1 2];
    for i = 1:numel(fns)
        tempmode = rezNoJaw(sessix).(fns{i});
        for j = 1:numel(cond)
            c = cond(j);

            tempdat = rezNoJaw(sessix).psth(:,:,c)*rezNoJaw(sessix).(fns{i});

            %             normfactor = abs(nanmean(tempdat(normTimes{i})));
            normfactor = 1;
            %             normfactor = max(tempdat);

            rezNoJaw(sessix).([fns{i}(1:end-5) '_latent'])(:,j) = tempdat ./ normfactor;
        end
    end
end
%% Plots with orthogonalizing to jaw mode
for i = 1:numel(objs)
    figure();
    obj = objs{i};
    met = meta(i);
    sesstitle = strcat(meta(i).anm, meta(i).date);

    go = mode(obj.bp.ev.goCue)-mode(obj.bp.ev.goCue);
    del = go-0.9;
    samp = del-1.3;

    conditions = {1,2};
    colors = {[0 0 1],[1 0 0]};

    subplot(4,1,2)
    plot(rezNoJaw(i).time,rezNoJaw(i).cdEarly_latent(:,1),'blue','LineWidth',2); hold on; plot(rezNoJaw(i).time,rezNoJaw(i).cdEarly_latent(:,2),'red','LineWidth',2);
    addtriallines(go,samp,del)
    legend('Right','Left')
    ylabel('a.u.')
    title('CDearly')
    xlim([-2.3 2.5])
    

    subplot(4,1,3)
    plot(rezNoJaw(i).time,rezNoJaw(i).cdLate_latent(:,1),'blue','LineWidth',2); hold on; plot(rezNoJaw(i).time,rezNoJaw(i).cdLate_latent(:,2),'red','LineWidth',2);
    addtriallines(go,samp,del)
    legend('Right','Left')
    ylabel('a.u.')
    title('CDlate')
    xlim([-2.3 2.5])

    subplot(4,1,4)
    plot(rezNoJaw(i).time,rezNoJaw(i).cdGo_latent(:,1),'blue','LineWidth',2); hold on; plot(rezNoJaw(i).time,rezNoJaw(i).cdGo_latent(:,2),'red','LineWidth',2);
    addtriallines(go,samp,del)
    legend('Right','Left')
    ylabel('a.u.')
    title('CDGo')
    xlim([-2.3 2.5])

    subplot(4,1,1)
    plot(rezNoJaw(i).time,rezNoJaw(i).cdjawVel_latent(:,1),'blue','LineWidth',2); hold on; plot(rezNoJaw(i).time,rezNoJaw(i).cdjawVel_latent(:,2),'red','LineWidth',2);
    addtriallines(go,samp,del)
    legend('Right','Left')
    ylabel('a.u.')
    title('CD Jaw Vel')
    xlim([-2.3 2.5])

    sgtitle(sesstitle)
end

%% Helper Functions


function cd = calcCD(rez,times)
tempdat = rez.psth(:,:,[1,2]);
mu = squeeze(mean(tempdat(times,:,:),1));
sd = squeeze(std(tempdat(times,:,:),[],1));
cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)
end

function addtriallines(go,samp,del)
xline(go,'black','LineStyle','--')
xline(del,'black','LineStyle','--')
xline(samp,'black','LineStyle','--')
end


