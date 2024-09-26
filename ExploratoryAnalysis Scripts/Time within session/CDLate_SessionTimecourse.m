% DECODING CDlate FROM ALL KINEMATIC FEATURES
clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\NullPotent_Context';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'fig1')));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Decoding Analysis'));
% rmpath(genpath(fullfile(utilspth,'fig1/')));
% rmpath(genpath(fullfile(utilspth,'mc_stim/')));

%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials
params.condition(end+1) = {'R&hit&~stim.enable&~autowater'};             % R 2AFC hits, no stim
params.condition(end+1) = {'L&hit&~stim.enable&~autowater'};             % L 2AFC hits, no stim
params.condition(end+1) = {'R&miss&~stim.enable&~autowater'};            % R error 2AFC, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater'};            % L error 2AFC, no stim

params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

params.advance_movement = 0;
%% SPECIFY DATA TO LOAD

datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';

meta = [];

% --- ALM ---
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
%meta = loadJEB13_ALMVideo(meta,datapth);
meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);

% --- M1TJ ---
% meta = loadJEB14_M1TJVideo(meta,datapth);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA

% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[obj,params] = loadSessionData(meta,params);

% ------------------------------------------
% -- Motion Energy --
% me (struct array) - one entry per session
% ------------------------------------------
for sessix = 1:numel(meta)
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
end
%% Calculate all CDs and find single trial projections
clearvars -except obj meta params me sav kin

disp('----Calculating coding dimensions----')
cond2use = [2,3];
cond2proj = [2,3];
regr = getCodingDimensions(obj,params,cond2use,cond2proj);

disp('----Projecting single trials onto CDlate----')
cd = 'early';
regr = getSingleTrialProjs(regr,obj,cd);
%% Divide CD projections on all trials into groups of chronological order (i.e. the first 20 trials, second 20 trials, third...)
ngroups = 4;
colorsR = {[0 0 1],[0.25 0.35 1],[0.5 0.5 1],[0.85 0.75 1]};
colorsL = {[1 0 0],[1 0.35 0.25],[1 0.5 0.5],[1 0.75 0.85]};
for sessix = 1:length(meta)
    singleproj_quartile = cell(1,length(cond2use));
    for c = 1:length(cond2use)
        cond = cond2use(c);
        condtrix = params(sessix).trialid{cond};
        temp = regr(sessix).singleProj(:,condtrix);

        nTrialsGroup = floor(length(condtrix)/ngroups);
        cnt = 1;
        groupavg = NaN(length(obj(1).time),ngroups);
        for g = 1:ngroups
           start = cnt;
           stop = cnt+nTrialsGroup;
           if stop > length(condtrix)
               stop = length(condtrix);
           end
           groupavg(:,g) = mySmooth(mean(temp(:,start:stop),2,'omitnan'),21);
           cnt = cnt+nTrialsGroup+1;
        end
        singleproj_quartile{c} = groupavg;
    end
        
%     figure();
%     for g = 1:ngroups
%         for c = 1:length(cond2use)
%             if c==1
%                 col = colorsR{g};
%             else
%                 col = colorsL{g};
%             end
%             plot(singleproj_quartile{c}(:,g),'Color',col,'LineWidth',1.5); hold on
%         end
%     end
%     legend('R first','L first','R second','L second','R third','L third','R last','L last')
temp = NaN(length(obj(1).time),ngroups);
for g = 1:ngroups
    temp(:,g) = singleproj_quartile{1}(:,g)-singleproj_quartile{2}(:,g);
end
selectivity(sessix).group = abs(temp);
end

%% Plot selectivity in CDs in each of the four quartiles
colors = {[0 0 0],[0.5 0.5 0.5], [0.67 0.67 0.67],[0.85 0.85 0.85]};
% for sessix = 1:length(meta)
%     figure()
%     for g = 1:ngroups
%         temp = selectivity(sessix).group(:,g);
%         plot(obj(1).time, temp,'Color',colors{g}); hold on;
%     end
%     legend('First','Second','Third','Last')
% end
%% Average the selectivity in CDs across all sessions
across_sess = cell(1,ngroups);
for g = 1:ngroups
    sel = [];
    for sessix = 1:length(meta)
        temp = selectivity(sessix).group(:,g);
        sel = [sel, temp];
    end
    across_sess{g} = sel;
end
%% Plot the average selectivity in CDs across all sessions
figure();
alpha = 0.2;
for g = 1:ngroups
    toplot = mean(across_sess{g},2,'omitnan');
    err = 1.96*(std(across_sess{g},0,2,'omitnan')/length(meta));
    ax = gca;
    shadedErrorBar(obj(1).time, toplot, err ,{'Color',colors{g},'LineWidth',2}, alpha, ax); hold on;
    xlabel('Time from go cue (s)')
    ylabel('Selectivity in CDlate (a.u.)')
    xlim([-2.5 2.5])
end

%%
function rez = getCodingDimensions(obj,params,cond2use,cond2proj)

cd_labels = {'early','late','go'};
cd_epochs = {'delay','goCue','goCue'};
cd_times = {[-0.42 -0.02], [-0.42 -0.02], [0.02 0.42]}; % in seconds, relative to respective epochs

for sessix = 1:numel(obj)

    %-------------------------------------------
    % --setup results struct--
    % ------------------------------------------
    rez(sessix).time = obj(sessix).time;
    rez(sessix).psth = standardizePSTH(obj(sessix));
    rez(sessix).condition = params(sessix).condition;
    rez(sessix).trialid = params(sessix).trialid;
    rez(sessix).alignEvent = params(sessix).alignEvent;
    rez(sessix).align = mode(obj(sessix).bp.ev.(rez(sessix).alignEvent));
    rez(sessix).ev.sample = obj(sessix).bp.ev.sample;
    rez(sessix).ev.delay = obj(sessix).bp.ev.delay;
    rez(sessix).ev.goCue = obj(sessix).bp.ev.goCue;


    % ------------------------------------------
    % --get coding directions--
    % ------------------------------------------
    rez(sessix).cd_mode = zeros(size(rez(sessix).psth,2),numel(cd_labels)); % (neurons,numCDs)
    for ix = 1:numel(cd_labels)
        % find time points to use
        e1 = mode(rez(sessix).ev.(cd_epochs{ix})) + cd_times{ix}(1) - rez(sessix).align;
        e2 = mode(rez(sessix).ev.(cd_epochs{ix})) + cd_times{ix}(2) - rez(sessix).align;
        times.(cd_labels{ix}) = rez(sessix).time>e1 & rez(sessix).time<e2;
        % calculate coding direction
        rez(sessix).cd_mode(:,ix) = calcCD(rez(sessix).psth,times.(cd_labels{ix}),cond2use);
    end


    % ------------------------------------------
    % --orthogonalize coding directions--
    % ------------------------------------------
    rez(sessix).cd_mode_orth = gschmidt(rez(sessix).cd_mode);


    % ------------------------------------------
    % --project neural population on CDs--
    % ------------------------------------------
    temp = permute(rez(sessix).psth(:,:,cond2proj),[1 3 2]); % (time,cond,neurons), permuting to use tensorprod() on next line for the projection
    rez(sessix).cd_proj = tensorprod(temp,rez(sessix).cd_mode_orth,3,1); % (time,cond,cd), cond is in same order as con2use variable defined at the top of this function


    % ------------------------------------------
    % --variance explained--
    % ------------------------------------------
    psth = rez(sessix).psth(:,:,cond2use);
    datacov = cov(cat(1,psth(:,:,1),psth(:,:,2)));
    datacov(isnan(datacov)) = 0;
    eigsum = sum(eig(datacov));
    for i = 1:numel(cd_labels)
        % whole trial
        rez(sessix).cd_varexp(i) = var_proj(rez(sessix).cd_mode_orth(:,i), datacov, eigsum);
        % respective epoch
        epoch_psth = rez(sessix).psth(times.(cd_labels{i}),:,cond2use);
        epoch_datacov = cov(cat(1,epoch_psth(:,:,1),epoch_psth(:,:,2)));
        epoch_datacov(isnan(epoch_datacov)) = 0;
        epoch_eigsum = sum(eig(epoch_datacov));
        rez(sessix).cd_varexp_epoch(i) = var_proj(rez(sessix).cd_mode_orth(:,i), epoch_datacov, epoch_eigsum);
    end

    % ------------------------------------------
    % --selectivity--
    % ------------------------------------------
    % coding directions
    rez(sessix).selectivity_squared = squeeze(rez(sessix).cd_proj(:,1,:) - rez(sessix).cd_proj(:,2,:)).^2; 
    % sum of coding directions
    rez(sessix).selectivity_squared(:,4) = sum(rez(sessix).selectivity_squared,2); 
    % full neural pop
    temp = rez(sessix).psth(:,:,cond2use);
    temp = (temp(:,:,1) - temp(:,:,2)).^2;
    rez(sessix).selectivity_squared(:,5) = sum(temp,2); % full neural pop

    % ------------------------------------------
    % --selectivity explained--
    % ------------------------------------------
    full = rez(sessix).selectivity_squared(:,5);
    rez(sessix).selexp = zeros(numel(rez(sessix).time),4); % (time,nCDs+1), +1 b/c sum of CDs
    for i = 1:4
        % whole trial
        rez(sessix).selexp(:,i) = 1 - ((full - rez(sessix).selectivity_squared(:,i)) ./ full);
    end


    % set some more rez variables to keep track of
    rez(sessix).cd_times = times;
    rez(sessix).cd_labels = cd_labels;
    rez(sessix).cd_epochs = cd_epochs;
    rez(sessix).cd_times_epoch = cd_times; % relative to respective epochs


end

end