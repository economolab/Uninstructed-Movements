clear,clc,close all

% finds null/potent spaces for a single session by :
% 1. split single trials into move and non-move time points
% 2. perform pca on non-move time points. Top PCs make up null space
% 3. subtract out null space from original data
% $. perform pca on residual data. Top PCs make up potent space


% * (most functions associated with this script are in the kinNullSpace directory) *

addpath(genpath(pwd))

%% PARAMETERS

% --SPECIFY WHICH ANIMAL AND SESSION TO LOAD
meta.anm = 'JEB6'; % 'JEB7'  'JEB6'  'JGR2'
meta.date = '2021-04-18'; % '2021-04-29'  '2021-04-18'  '2021-11-17'


% --SPECIFY PATH TO DATA HERE
if ispc
    meta.datapth = 'C:\Code\uninstructedMovements-Munib/data';
else
    meta.datapth = '/Users/Munib/Documents/Economo-Lab/data/';
end
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
params = getDefaultParams();
% params.probe = 2; % change default params.probe from '1' to '2'


% --SPECIFY METHOD OF DENOISING SINGLE TRIALS
% we need denoised, smooth single trial neural activity. We have two
% options:
% 1) load data that's ALREADY been passed through an lfads model
% 2) perform Factor Analysis on binned single trial data, followed by
%    smoothing
params.lfads_or_fa = 'lfads'; % 'lfads' or 'fa'
params.lfads_run = ''; % 'run3' , leave empty to use most recent run
params.fcut_post_fa = 31; % if performing FA, cutoff freq to smooth rates and factors with a butterworth filter
params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance
params.full_or_reduced = 'reduced'; % 'full'  or 'reduced' -- which data to use in regression
% using the full data will require another method, the system seems to be
% highly overfit as is
assert(strcmpi(params.full_or_reduced,'reduced'),'method to use full dimensional data doesnt exist, params.full_or_reduced should be set to `reduced`')

% --SPECIFY TIME POINTS AND LAG TO USE
params.prep = [-2.5 -0.05]; % initial and final time points (seconds) defining prep epoch, relative to alignevent
params.move = [-2.5 1.5];   % initial and final time points (seconds) defining move epoch, relative to alignevent
params.advance_movement = 0.025; % seconds, amount of time to advance movement data relative to neural data

%% NEURAL ACTIVITY

% getNeuralActivity() returns 4 main variables
% - dat: contains lfads/fa smoothed firing rates, factors, and trial numbers
%        dat.factors and dat.rates are size (time,factors/clusters,trials)
% - meta: session meta data
% - params: parameters used for preprocessing lfads input data
% - obj: preprocessed data obj
[meta,params,obj,dat] = getNeuralActivity(meta,params);
dat.factors = standardizeFactors(dat.factors);

% for i = 1:size(dat.factors,2)
%     figure(i)
% %     plot(obj.time,squeeze(dat.factors(:,i,1:numel(params.trialid{1}))))
%     imagesc(obj.time,1:numel(dat.trials),squeeze(dat.factors(:,i,:))');
% end
% 
% 
% % factor 8 shows differences within a condition, first half and second half
% figure;
% subplot(2,1,1)
% imagesc(obj.time,1:numel(params.trialid{1}),squeeze(dat.factors(:,8,1:numel(params.trialid{1})))');
% title('JEB7 4-29, Factor 8, Right Hit Trials')
% subplot(2,1,2)
% plot(obj.time,mean(temp(:,1:37),2),'LineWidth',2);
% hold on;
% plot(obj.time,mean(temp(:,38:75),2),'LineWidth',2)
% legend('First Half, avg right hit trials','Second Half, avg right hit trials')


%% MOTION ENERGY
% me is a struct with fields
% - data:       cell array (ntrials,1), with motion energy for each time point
% - moveThresh: a scalar - below this value of motion energy, animal is
%               said to be stationary
% - use:        0 or 1, indicating whether to use me.data. 1 if
%               a motionEnergy*.mat file is found, 0 if file not found
me = loadMotionEnergy(obj,meta,params,dat.trials); 

% To generate a motionEnergy*.mat file for a session, see https://github.com/economolab/videoAnalysisScripts/blob/main/motionEnergy.m

%%

close all

trialOffset = 0;
f = figure; hold on;
f.Position = [-1323        -145         574         968];
for i = 1:30:size(me.data,2)
    temp = mySmooth(me.data(:,i),21);
    
    ix = me.data(:,i)>(me.moveThresh+5);
    z = temp;
    z(~ix) = nan;
    ztime = obj.time;
    patchline(obj.time,trialOffset + temp,'EdgeColor','k','EdgeAlpha',0.35,'LineWidth',4);
    plot(obj.time,trialOffset + z,'r','LineWidth',2)
    trialOffset = trialOffset + 90;
end
xlabel('Time (s) from go cue')
ylabel('Motion Energy')
xlim([obj.time(15),obj.time(end)]);
ax = gca;
ax.YTick = [];
ax.FontSize = 35;

align = mode(obj.bp.ev.(params.alignEvent));
sample = mode(obj.bp.ev.sample) - align;
delay = mode(obj.bp.ev.delay) - align;
xline(sample,'k--','LineWidth',2)
xline(delay,'k--','LineWidth',2)
xline(0,'k--','LineWidth',2)

hold off

% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4/figs/motionEnergy/';
% fn = ['motionEnergy_JEB7_2021-04-29_threshincreasedby5'];
% mysavefig(f,pth,fn);


%%
close all

data = me.data;

z = data - (me.moveThresh+5);
z(z<=0) = -max(max(z));

f = figure; hold on;
f.Position = [-1311        -179         798        1003];
imagesc(obj.time,1:size(z,2),z')

xline(sample,'k--','LineWidth',3)
xline(delay,'k--','LineWidth',3)
xline(0,'k--','LineWidth',3)

hold off

ylim([1,size(z,2)])

a = colorbar;
a.Label.String = '(Motion Energy - Threshold) > 0';
xlabel('Time (s) from go cue')
ylabel('Trials')
ax = gca;
ax.FontSize = 30;

colormap cool
% 
% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4/figs/motionEnergy/';
% fn = ['motionEnergyMap_JEB7_2021-04-29_threshincreasedby5'];
% mysavefig(f,pth,fn);












