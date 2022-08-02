clear,clc,close all


addpath(genpath(pwd))

%% PARAMETERS

% --SPECIFY WHICH ANIMAL AND SESSION TO LOAD
meta(1).anm = 'JEB7'; 
meta(1).date = '2021-04-29'; 

meta(end+1).anm = 'JEB7'; 
meta(end).date = '2021-04-30'; 

meta(end+1).anm = 'JEB6'; 
meta(end).date = '2021-04-18'; 

meta(end+1).anm = 'JGR2'; 
meta(end).date = '2021-11-16';

meta(end+1).anm = 'JGR2'; 
meta(end).date = '2021-11-17';

meta(end+1).anm = 'JGR3'; 
meta(end).date = '2021-11-18';

meta(end+1).anm = 'EKH3'; 
meta(end).date = '2021-08-11';

% meta(end+1).anm = 'EKH3'; 
% meta(end).date = '2021-08-07';

meta(end+1).anm = 'EKH1'; 
meta(end).date = '2021-08-07';









meta = assignDataPath(meta);

% the data should be stored in the following structure:
% /params.datapth/
%  --- /DataObjects/
%  --- --- /animal_name_number/ (contains data_structure*.mat)
%  --- /lfads/
%  --- --- /input/ (contains anm_name_session_date*.mat/.h5)
%  --- --- /output/ (contains model_train/valid_anm_name_session_date*.h5)


% --SPECIFY PREPROCESSING PARAMETERS YOU WANT TO CHANGE BELOW THE FUNCTION CALL
% if loading lfads data, these default params will be overwritten by params
% used to generate lfads input data
dfparams = getDefaultParams();
% params.probe = 2; % change default params.probe from '1' to '2'


% --SPECIFY METHOD OF DENOISING SINGLE TRIALS
% we need denoised, smooth single trial neural activity. We have two
% options:
% 1) load data that's ALREADY been passed through an lfads model
% 2) perform Factor Analysis on binned single trial data, followed by
%    smoothing
dfparams.lfads_or_fa = 'lfads'; % 'lfads' or 'fa'
dfparams.lfads_run = 'run1'; % 'run3' , leave empty to use most recent run
dfparams.fcut_post_fa = 31; % if performing FA, cutoff freq to smooth rates and factors with a butterworth filter
dfparams.feat_varToExplain = 75; % num factors for dim reduction of video features should explain this much variance
dfparams.full_or_reduced = 'reduced'; % 'full'  or 'reduced' -- which data to use in regression
% using the full data will require another method, the system seems to be
% highly overfit as is
assert(strcmpi(dfparams.full_or_reduced,'reduced'),'method to use full dimensional data doesnt exist, params.full_or_reduced should be set to `reduced`')

% --SPECIFY TIME POINTS AND LAG TO USE
dfparams.prep = [-2.5 -0.05]; % initial and final time points (seconds) defining prep epoch, relative to alignevent
dfparams.move = [-2.5 1.5];   % initial and final time points (seconds) defining move epoch, relative to alignevent
dfparams.advance_movement = 0.025; % seconds, amount of time to advance movement data relative to neural data

%% NEURAL ACTIVITY

use = true(size(meta));
for i = 1:numel(meta)
    disp(['Loading data for ' meta(i).anm ' ' meta(i).date]);
%     [meta(i),params(i),obj(i),dat(i)] = getNeuralActivity(meta(i),dfparams);
    fa(i) = getFAData(meta(i),'run2');
    params(i) = fa(i).params;
    if isfield(fa(i).obj,'meta')
        fa(i).obj = rmfield(fa(i).obj,'meta');
    end
%     obj(i) = fa(i).obj;
%     me(i) = loadMotionEnergy(obj(i),meta(i),params(i),1:obj(i).bp.Ntrials); 
%     if ~me(i).use
%         use(i) = false;
%     end
    disp('DONE');
    disp(' ');
end

fa = fa(use);
params = params(use);
meta = meta(use);
obj = obj(use);
me = me(use);

%% SET RUN PARAMS
clear params

params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to 

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                                % all trials
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};         % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};         % left hits, no stim, aw off
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};        % error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};        % error left, no stim, aw off


params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 0.02;

% smooth with causal gaussian kernel
params.smooth = 0;

% cluster qualities to use
% params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality
params.quality = {'all'}; 

%% SET METADATA

% datapth = '/Users/Munib/Documents/Economo-Lab/data/';
datapth = 'M:/Economo-Lab/data/';

meta = [];
% meta = loadJEB4_ALMVideo(meta,datapth); % done
% meta = loadJEB5_ALMVideo(meta,datapth); % done
meta = loadJEB7_ALMVideo(meta,datapth); % done
meta = loadJEB6_ALMVideo(meta,datapth); % done
meta = loadJGR2_ALMVideo(meta,datapth); % done
meta = loadJGR3_ALMVideo(meta,datapth); % done
meta = loadEKH3_ALMVideo(meta,datapth); % done
meta = loadEKH1_ALMVideo(meta,datapth); % done


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

fa = fa(use);
meta = meta(use);
objs = objs(use);
params.probe = params.probe(use);
params.trialid = params.trialid(use);
params.cluid = params.cluid(use);

%% Motion Energy

use = true(size(meta));
for i = 1:numel(meta)
    me(i) = loadMotionEnergy(objs{i},meta(i),params,1:objs{i}.bp.Ntrials); 
    if ~me(i).use
        use(i) = false;
    end
    disp('DONE');
    disp(' ');
end


params.probe = params.probe(use);
params.trialid = params.trialid(use);
params.cluid = params.cluid(use);
meta = meta(use);
objs = objs(use);
me = me(use);
fa = fa(use);


%%

clearvars -except meta params obj objs dat fa dfparams me kin kinfeats kinfeats_reduced


%% VIDEO FEATURES

% getKinematicsFromVideo() returns 2 variables
% - kin: feature matrix of size (time,trials,features). Features defined in
%         params.traj_features
% - featLeg: legend corresponding to features in kin (for 2nd dimension)
for i = 1:numel(meta)
    [kin(i).dat,kin(i).featLeg] = getKinematicsFromVideo(objs{i},dfparams,params,1:objs{i}.bp.Ntrials);
    
%     % TONGUE ANGLE AND LENGTH % TODO
    [ang,len] = getLickAngleAndLength(kin(i).featLeg,kin(i).dat);
    kin(i).featLeg{end+1} = 'tongue_angle';
    kin(i).featLeg{end+1} = 'tongue_length';
    
    kinfeats{i} = kin(i).dat;
%     % create feature matrix, feats, and assign to a field in dat
%     kinfeats{i} = cat(3,kin(i).dat,reshape(ang,size(ang,1),size(ang,2),1));
%     kinfeats{i} = cat(3,kinfeats{i},reshape(len,size(len,1),size(len,2),1));
    
    
    % MOTION ENERGY
    if me(i).use
        kinfeats{i} = cat(3,kinfeats{i},reshape(me(i).data,size(me(i).data,1),size(me(i).data,2),1));
        kin(i).featLeg{end+1} = 'motion_energy';
    end
    % To generate a motionEnergy*.mat file for a session, see https://github.com/economolab/videoAnalysisScripts/blob/main/motionEnergy.m
    
    
    % STANDARDIZE FEATURES (ZERO MEAN, UNIT VARIANCE)
    kinfeats{i} = standardizeFeatures(kinfeats{i});
    
    
    % DIMENSIONALITY REDUCTION
    % many of the video features will be highly correlated, so we will perform PCA/FA
    % on the matrix of features to reduce the dimensionality to a set of factors that
    % best explain the movement captured by the video recordings
    kinfeats_reduced{i} = reduceDimensionVideoFeatures(kinfeats{i},dfparams.feat_varToExplain);
    
end



disp('DONE CREATING FEATURE MATRIX AND REDUCED DIM FEATURE MATRIX')



%% REGRESSION

% here we pose the problem
%           
%               V = NW
%
% where V is a matrix of size (numVideoFeatures , time*trials)'
% and   N is a matrix of size (numFactors/numClusters , time*trials)'
% 
% We will estimate W with ridge regression (regularized linear regression).
% It is the linear transformation that takes neural activity to video
% activity
% The null space of W forms subspace that we will call the         NULL SPACE   of N
% The column space of W forms the subspace that we will call the   POTENT SPACE of N

for i = 1:numel(meta)
    clear newrez
    input_data.factors = fa(i).falatents;
    input_data.feats = kinfeats_reduced{i};
    newrez = estimateW(objs{i},input_data,dfparams,objs{i}.time,me(i)); % N,V are zscored neural activity and feature matrix
    
    
    
    % NULL AND POTENT SPACE OF W
    
    [newrez.W_null,newrez.W_potent,newrez.N_null,newrez.N_potent,newrez.dPrep,newrez.dMove] = ... 
        getNullPotentSpaces_SVD(newrez.W',newrez,input_data); % transposing W b/c we found W by solving V'=N'*W, rather than V = WN
    
    % correlation between V and V_hat = rez.N * rez.W
    % rez.corrcoef = calcCorrCoef(rez.V,rez.N,rez.W);
    
    % variance explained by null and potent space
    % (this is how much variance is explained out of max variance to be
    % explained from total number of null and potent dims)
    % TODO: var explained of trial averaged data, not single trials b/c the
    % variability messes things up...
    verez = calVarExp(newrez);
    
    rezfns = fieldnames(verez);
    for j = 1:numel(rezfns)
        rez(i).(rezfns{j}) = verez.(rezfns{j});
    end
end

%% variance explained plots
close all

null_total = zeros(size(rez));
potent_total = zeros(size(rez));
for i = 1:numel(rez)
    null_total(i) = rez(i).varexp_null;
    potent_total(i) = rez(i).varexp_potent;
end

violincols = [50, 168, 82; 168, 50, 142] ./ 225;
varexp_full = [null_total ; potent_total]';
f = figure; ax = axes(f);
vs = violinplot(varexp_full,{'Null','Potent'},...
    'EdgeColor',[1 1 1], 'ViolinAlpha',{0.2,1}, 'ViolinColor', violincols);
ylabel('Normalized Variance Explained (Whole Trial)')
ylim([0,1])
ax = gca;
ax.FontSize = 20;
% % 
% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/kinNullSpace';
% fn = 've_total';
% mysavefig(f,pth,fn);


%% plot projections

plotProjections(params,obj,rez)


%%

sessix = 1;
potenttemp = rez(sessix).N_potent;
kintemp = kinfeats_reduced{sessix};

metemp = me(sessix);

close all

nTrials = 30;
trix = randsample(size(metemp.data,2),nTrials);

medata = metemp.data(:,trix);
potentdim = 1;
potentdata = potenttemp(:,trix,potentdim);

kindim = 3;
kindata = kintemp(:,trix,kindim);

medata = medata(:);

potentdata = potentdata(:);

kindata = kindata(:);

newtime = (1:numel(potentdata))./200;
figure; hold on;
patchline(newtime,medata,'EdgeColor','k','EdgeAlpha',0.35,'LineWidth',2);
ix = medata > metemp.moveThresh;
z = medata;
z(~ix) = nan;
plot(newtime,z,'r','LineWidth',2)
% plot(newtime,potentdata*45,'b','LineWidth',1);
plot(newtime,kindata*3,'b','LineWidth',1);

%%


trialOffset = 0;
f = figure; hold on;
% f.Position = [-1323        -145         574         968];
for i = 1:30:size(metemp.data,2)
    temp = mySmooth(metemp.data(:,i),21);
    
    ix = metemp.data(:,i)>(metemp.moveThresh);
    z = temp;
    z(~ix) = nan;
    ztime = obj(sessix).time;
    patchline(obj(sessix).time,trialOffset + temp,'EdgeColor','k','EdgeAlpha',0.35,'LineWidth',4);
    plot(obj(sessix).time,trialOffset + z,'r','LineWidth',2)

    plot(obj(sessix).time, trialOffset + potenttemp(:,i,1)*45, 'b', 'LineWidth', 0.5)


    trialOffset = trialOffset + 90;
end
xlabel('Time (s) from go cue')
ylabel('Motion Energy')
xlim([obj(sessix).time(15),obj(sessix).time(end)]);
ax = gca;
% ax.YTick = [];
ax.FontSize = 35;
% 
% align = mode(obj.bp.ev.(params.alignEvent));
% sample = mode(obj.bp.ev.sample) - align;
% delay = mode(obj.bp.ev.delay) - align;
% xline(sample,'k--','LineWidth',2)
% xline(delay,'k--','LineWidth',2)
% xline(0,'k--','LineWidth',2)

hold off

%% activity modes

clear cdrez times

for i = 1:numel(rez)
    [cdrez(i),times] = cdNullSpace_elsayed(rez(i),objs{i},params(i));
end

rez = cdrez;


%% null space cds

nullSpaceCD(rez,obj,params,times)


%% potent space cds

potentSpaceCD(rez,obj,params,times)





