clear,clc,close all

addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Data-Loading-Scripts'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements'));
% finds cd early, late, go as defined in economo 2018



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

%% SET RUN PARAMS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'
params.kinfind             = 'vel';
% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val
params.moveThresh          = 0.15;

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

% Params for finding kinematic modes
params.fcut = 50;          % smoothing cutoff frequency
params.cond = 1:2;         % which conditions to use to find mode
params.method = 'xcorr';   % 'xcorr' or 'regress' (basically the same)
params.fa = false;         % if true, reduces neural dimensions to 10 with factor analysis


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
%%
% Remove unwanted sessions
[meta,objs,params] = useInclusionCriteria(objs,params,meta);
%% Adjust for older functions
for i = 1:numel(objs)
    meta(i).trialid = params.trialid{i};
    conditions = {1,2};

    for c = 1:numel(conditions)
        trix = meta(i).trialid{c};
        objs{i}.trialpsth_cond{c} = objs{i}.trialdat(:,:,trix);
    end
end
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

    % Calculate all CDs

    [Early, Late, Go] = findAllCDs(rez,sessix,ev,params);
    rez(sessix).cdEarly_mode = Early; rez(sessix).cdLate_mode = Late;  rez(sessix).cdGo_mode = Go;

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

    cond = [1 2];
    for i = 1:numel(fns)
        tempmode = rez(sessix).(fns{i});
        for j = 1:numel(cond)
            c = cond(j);
            tempdat = rez(sessix).psth(:,:,c)*rez(sessix).(fns{i});
            normfactor = 1;
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

    for i = 1:numel(fns)
        temp = rez(sessix).([fns{i}(1:end-5) '_latent']);
        selectivity(sessix).(fns{i}(1:end-5)) = temp(:,1) - temp(:,2);
    end


end
%% Without removing jaw mode
for i = 1:numel(objs)
    figure();
    obj = objs{i};
    met = meta(i);
    sesstitle = strcat(meta(i).anm, meta(i).date);

    plotModesAndJaw(obj,met,rez,i,params,sesstitle)   
end
%%  Find all kin modes; orthogonalize other modes to it
rezNoJaw = rez;
for gg = 1:numel(objs)
    obj = objs{gg};
    met = meta(gg);
    sesstitle = strcat(meta(gg).anm, meta(gg).date);
    taxis = obj.time;

    %%% Calculate all modes (based on late delay) %%%
    psthForProj = cat(3,obj.trialpsth_cond{1},obj.trialpsth_cond{2});       % Concatenated single trial PSTHs (all trials from cond 1 then all trials from cond 2)    

    %%%% REMOVE HALF RIGHT TRIALS AND HALF LEFT TRIALS AT RANDOM
    % Only affects the trials being used to calculate the kinematic modes,
    % not the projections
    removehalf = 'no';
    if strcmp(removehalf,'yes')
        for m = 1:numel(met.trialid)
            nTrials = length(met.trialid{m});
            indx = randsample(nTrials,round(nTrials/2));
            met.trialid{m} = met.trialid{m}(indx);
            obj.trialpsth_cond{m} = obj.trialpsth_cond{m}(:,:,indx);
        end
    end

    %%% GET FEATURE KINEMATICS %%%
    conditions = {1,2};
    kin = struct();
    kin.MEinterp = getME(obj,met,params,taxis,conditions);                  % Find interpolated motion energy
    kinfeat = findAllFeatKin(params, taxis, obj,conditions,met);
    kin.jawVel = kinfeat.jawVel; kin.noseVel = kinfeat.noseVel; kin.tongueVel = kinfeat.tongueVel;
    kin.tongueAngle = findTongueAngle(taxis, obj, met,params,conditions);

    conditions = 1:2;
    kinfns = fieldnames(kin);
    for i = 1:numel(kinfns)
        Y = kin.(kinfns{i}); % feature data to use to calculate mode
        if ~strcmp(kinfns{i},'tongueVel') && ~strcmp(kinfns{i},'tongueAngle')
            e1 = find(taxis>-0.5,1,'first');
            e2 = find(taxis>-0.05,1,'first');
            params.tix = e1:e2;        % time points to use when finding mode (LATE DELAY)
        elseif strcmp(kinfns{i},'tongueVel') || strcmp(kinfns{i},'tongueAngle')
            e1 = find(taxis>0.05,1,'first');
            e2 = find(taxis>0.5,1,'first');
            params.tix = e1:e2;        % time points to use when finding mode (RESPONSE PERIOD)
        end
        [kinmode.(kinfns{i}), dat.(kinfns{i})] = findMode(obj, Y, params,conditions);
        modename = strcat('cd',kinfns{i},'_mode');
        rezNoJaw(gg).(modename) = kinmode.(kinfns{i});
    end

    

    %[fns,~] = patternMatchCellArray(fieldnames(rezNoJaw(sessix)),{'mode'},'all');
    tonguenan = sum(isnan(rezNoJaw(gg).cdtongueAngle_mode));            % Account for the tongue mode being NaNs
    nosenan = sum(isnan(rezNoJaw(gg).cdnoseVel_mode));            % Account for the tongue mode being NaNs
    if tonguenan==0 && nosenan==0
        fns = {'cdMEinterp_mode';'cdjawVel_mode';'cdnoseVel_mode';'cdtongueVel_mode';'cdtongueAngle_mode';'cdEarly_mode';'cdLate_mode';'cdGo_mode'};
    elseif nosenan==0 && tonguenan>0
        fns = {'cdMEinterp_mode';'cdjawVel_mode';'cdnoseVel_mode';'cdEarly_mode';'cdLate_mode';'cdGo_mode'};
    else 
        fns = {'cdMEinterp_mode';'cdjawVel_mode';'cdnoseVel_mode';'cdEarly_mode';'cdLate_mode';'cdGo_mode'};
    end

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

    % Variance explained
    for i = 1:numel(fns)
        psth = rezNoJaw(sessix).psth;
        datacov = cov([psth(:,:,1) ; psth(:,:,2)]);
        datacov(isnan(datacov)) = 0;
        eigsum = sum(eig(datacov));
        rezNoJaw(sessix).varexp.(fns{i}(1:end-5)) = var_proj(rezNoJaw(sessix).(fns{i}), datacov, eigsum);
    end

     for i = 1:numel(fns)
        temp = rezNoJaw(sessix).([fns{i}(1:end-5) '_latent']);
        selectivityNoJaw(sessix).(fns{i}(1:end-5)) = temp(:,1) - temp(:,2);
    end
end
%%
%Plots with orthogonalizing to jaw mode
for i = 1:numel(objs)
    figure();rez(gg).varexp
    obj = objs{i};
    met = meta(i);
    sesstitle = strcat(meta(i).anm, meta(i).date);

    plotJawModeandCDs(obj,rezNoJaw,i,sesstitle)
end
%%
VE.early = [];  VEnoJaw.early = [];
VE.late = [];  VEnoJaw.late = [];
VE.go = [];  VEnoJaw.go = [];
for gg = 1:length(meta)
    VE.early = [VE.early, rez(gg).varexp.cdEarly];  VEnoJaw.early = [VEnoJaw.early, rezNoJaw(gg).varexp.cdEarly];
    VE.late = [VE.late, rez(gg).varexp.cdLate];  VEnoJaw.late = [VEnoJaw.late, rezNoJaw(gg).varexp.cdLate];
    VE.go = [VE.go, rez(gg).varexp.cdGo];  VEnoJaw.go = [VEnoJaw.go, rezNoJaw(gg).varexp.cdGo];
end
%%
% figure();
% x = [1,2,4,5,7,8];
% y = [mean(VE.early),mean(VEnoJaw.early,'omitnan'),mean(VE.late),mean(VEnoJaw.late,'omitnan'),mean(VE.go),mean(VEnoJaw.go,'omitnan')];
% b= bar(x,y);
% b.FaceColor = [0.75 0.75 0.75]; hold on;
% 
% scatter(1,VE.early,'black','filled'); scatter(2,VEnoJaw.early,'magenta','filled')
% scatter(4,VE.late,'black','filled');  scatter(5,VEnoJaw.late,'magenta','filled')
% scatter(7,VE.go,'black','filled');  scatter(8,VEnoJaw.go,'magenta','filled')
% yy = [VE.early',VEnoJaw.early']; plot(x(1:2),yy(:,1:2),'Color','black')
% yy = [VE.late',VEnoJaw.late']; plot(x(3:4),yy(:,1:2),'Color','black')
% yy = [VE.go',VEnoJaw.go']; plot(x(5:6),yy(:,1:2),'Color','black')
% 
% xlim([0 9])
% ylabel('Variance explained (out of 1)')
% title('Difference in VE after orthogonalizing to KinModes')
%%
sel.early = [];  selnoJaw.early = [];
sel.late = [];  selnoJaw.late = [];
sel.go = [];  selnoJaw.go = [];
for gg = 1:length(meta)
    sel.early = [sel.early, selectivity(gg).cdEarly];  selnoJaw.early = [selnoJaw.early, selectivityNoJaw(gg).cdEarly];
    sel.late = [sel.late, selectivity(gg).cdLate];  selnoJaw.late = [selnoJaw.late, selectivityNoJaw(gg).cdLate];
    sel.go = [sel.go, selectivity(gg).cdGo];  selnoJaw.go = [selnoJaw.go, selectivityNoJaw(gg).cdGo];
end

e1 = find(taxis>-1.5,1,'first');
e2 = find(taxis>-0.8,1,'first');
sel.early = abs(mean(sel.early(e1:e2,:),1));  selnoJaw.early = abs(mean(selnoJaw.early(e1:e2,:),1));  

e1 = find(taxis>-0.5,1,'first');
e2 = find(taxis>-0.05,1,'first');
sel.late = abs(mean(sel.late(e1:e2,:),1));  selnoJaw.late = abs(mean(selnoJaw.late(e1:e2,:),1));  

e1 = find(taxis>0.05,1,'first');
e2 = find(taxis>0.5,1,'first');
sel.go = abs(mean(sel.go(e1:e2,:),1));  selnoJaw.go = abs(mean(selnoJaw.go(e1:e2,:),1));  
%%
% figure();
% x = [1,2,4,5,7,8];
% y = [mean(sel.early),mean(selnoJaw.early,'omitnan'),mean(sel.late),mean(selnoJaw.late,'omitnan'),mean(sel.go),mean(selnoJaw.go,'omitnan')];
% b= bar(x,y);
% b.FaceColor = [0.75 0.75 0.75]; hold on;
% 
% scatter(1,sel.early,'black','filled'); scatter(2,selnoJaw.early,'blue','filled')
% scatter(4,sel.late,'black','filled');  scatter(5,selnoJaw.late,'blue','filled')
% scatter(7,sel.go,'black','filled');  scatter(8,selnoJaw.go,'blue','filled')
% yy = [sel.early',selnoJaw.early']; plot(x(1:2),yy(:,1:2),'Color','black')
% yy = [sel.late',selnoJaw.late']; plot(x(3:4),yy(:,1:2),'Color','black')
% yy = [sel.go',selnoJaw.go']; plot(x(5:6),yy(:,1:2),'Color','black')
% 
% xlim([0 9])
% ylabel('Avg selectivity (a.u)')
% title('Difference in Selectivity after orthogonalizing to KinModes')
%% Helper Functions

function MEinterp = getME(obj,met,params,taxis,conditions)
[met,mov,me] = assignEarlyTrials(obj,met,params);
[temp,~] = findInterpME(taxis,conditions, met,mov,me,params,obj);
MEinterp = [];
for c = 1:numel(conditions)
    MEinterp = [MEinterp, temp{c}];
end
end

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

function [cdEarly_mode, cdLate_mode, cdGo_mode] = findAllCDs(rez,sessix,ev,params)

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
    cdEarly_mode = calcCD(rez(sessix),times.early);


    % cd late mode

    e1 = mode(ev.goCue) - 0.5 - mode(ev.(params.alignEvent));
    e2 = mode(ev.goCue) - 0.1 - mode(ev.(params.alignEvent));

    times.late = rez(sessix).time>e1 & rez(sessix).time<e2;
    cdLate_mode = calcCD(rez(sessix),times.late);


    % cd go mode

    e1 = mode(ev.goCue) + 0.02 - mode(ev.(params.alignEvent));
    e2 = mode(ev.goCue) + 0.42 - mode(ev.(params.alignEvent));

    times.go = rez(sessix).time>e1 & rez(sessix).time<e2;
    cdGo_mode = calcCD(rez(sessix),times.go);

    % cd after mode

    %     e1 = mode(ev.goCue) + 1.5 - mode(ev.(params.alignEvent));
    %     e2 = mode(ev.goCue) + 2.0 - mode(ev.(params.alignEvent));
    %
    %     times.after = rez(sessix).time>e1 & rez(sessix).time<e2;
    %     rez(sessix).cdAfter_mode = calcCD(rez(sessix),times.after);
end
%%
function plotModesAndJaw(obj,met,rez,i,params,sesstitle)

go = mode(obj.bp.ev.goCue)-mode(obj.bp.ev.goCue);
    del = go-0.9;
    samp = del-1.3;

    conditions = {1,2};
    colors = {[0 0 1],[1 0 0]};

    subplot(3,1,1)
    plot(rez(i).time,rez(i).cdEarly_latent(:,1),'blue','LineWidth',2); hold on; plot(rez(i).time,rez(i).cdEarly_latent(:,2),'red','LineWidth',2);
    addtriallines(go,samp,del)
    legend('Right','Left')
    ylabel('a.u.')
    currVE = num2str(rez(i).varexp.cdEarly);
    plottitle = strcat('CDearly; VE = ',{' '},currVE);
    title(plottitle)
    xlim([-2.3 2.5])


    subplot(3,1,2)
    plot(rez(i).time,rez(i).cdLate_latent(:,1),'blue','LineWidth',2); hold on; plot(rez(i).time,rez(i).cdLate_latent(:,2),'red','LineWidth',2);
    addtriallines(go,samp,del)
    legend('Right','Left')
    ylabel('a.u.')
    currVE = num2str(rez(i).varexp.cdLate);
    plottitle = strcat('CDlate; VE = ',{' '},currVE);
    title(plottitle)
    xlim([-2.3 2.5])

    subplot(3,1,3)
    plot(rez(i).time,rez(i).cdGo_latent(:,1),'blue','LineWidth',2); hold on; plot(rez(i).time,rez(i).cdGo_latent(:,2),'red','LineWidth',2);
    addtriallines(go,samp,del)
    legend('Right','Left')
    ylabel('a.u.')
    currVE = num2str(rez(i).varexp.cdGo);
    plottitle = strcat('CDgo; VE = ',{' '},currVE);
    title(plottitle)
    xlim([-2.3 2.5])


%     subplot(4,1,1)
%     taxis = rez(i).time;
% 
%     jaw_by_cond = findJawVelocity(taxis, obj,conditions,met,'vel',params);    % (1 x conditions cell array)
%     % Each cell: (time x trials in that condition)
%     jawvel.right = mean(jaw_by_cond{1},2,'omitnan');
%     jawvel.left = mean(jaw_by_cond{2},2,'omitnan');
%     plot(taxis, jawvel.right,'blue','LineWidth',2); hold on;
%     plot(taxis, jawvel.left,'red','LineWidth',2);
% 
%     ylabel('Jaw velocity')
%     addtriallines(go,samp,del)
%     legend('Right','Left','Location','best')
%     title('Jaw')
%     xlim([-2.3 2.5])

    sgtitle(sesstitle)
end

function plotJawModeandCDs(obj,rezNoJaw,i,sesstitle)
go = mode(obj.bp.ev.goCue)-mode(obj.bp.ev.goCue);
    del = go-0.9;
    samp = del-1.3;

    conditions = {1,2};
    colors = {[0.5 0.5 0.8],[0.8 0.5 0.5]};

    subplot(3,1,1)
    plot(rezNoJaw(i).time,rezNoJaw(i).cdEarly_latent(:,1),'Color',colors{1},'LineWidth',2); hold on; plot(rezNoJaw(i).time,rezNoJaw(i).cdEarly_latent(:,2),'Color',colors{2},'LineWidth',2);
    addtriallines(go,samp,del)
    legend('Right','Left')
    ylabel('a.u.')
    currVE = num2str(rezNoJaw(i).varexp.cdEarly);
    plottitle = strcat('CDearly; VE = ',{' '},currVE);
    title(plottitle)
    xlim([-2.3 2.5])


    subplot(3,1,2)
    plot(rezNoJaw(i).time,rezNoJaw(i).cdLate_latent(:,1),'Color',colors{1},'LineWidth',2); hold on; plot(rezNoJaw(i).time,rezNoJaw(i).cdLate_latent(:,2),'Color',colors{2},'LineWidth',2);
    addtriallines(go,samp,del)
    legend('Right','Left')
    ylabel('a.u.')
    currVE = num2str(rezNoJaw(i).varexp.cdLate);
    plottitle = strcat('CDlate; VE = ',{' '},currVE);
    title(plottitle)
    xlim([-2.3 2.5])

    subplot(3,1,3)
    plot(rezNoJaw(i).time,rezNoJaw(i).cdGo_latent(:,1),'Color',colors{1},'LineWidth',2); hold on; plot(rezNoJaw(i).time,rezNoJaw(i).cdGo_latent(:,2),'Color',colors{2},'LineWidth',2);
    addtriallines(go,samp,del)
    legend('Right','Left')
    ylabel('a.u.')
    currVE = num2str(rezNoJaw(i).varexp.cdGo);
    plottitle = strcat('CDgo; VE = ',{' '},currVE);
    title(plottitle)
    xlim([-2.3 2.5])

%     subplot(4,1,1)
%     plot(rezNoJaw(i).time,rezNoJaw(i).cdjawVel_latent(:,1),'Color',colors{1},'LineWidth',2); hold on; plot(rezNoJaw(i).time,rezNoJaw(i).cdjawVel_latent(:,2),'Color',colors{2},'LineWidth',2);
%     addtriallines(go,samp,del)
%     legend('Right','Left')
%     ylabel('a.u.')
%     currVE = num2str(rezNoJaw(i).varexp.cdjawVel);
%     plottitle = strcat('JawVel Mode; VE = ',{' '},currVE);
%     title(plottitle)
%     xlim([-2.3 2.5])

    sgtitle(sesstitle)
end