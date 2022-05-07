clear,clc,close all


addpath(genpath(pwd))

%% PARAMETERS

% --SPECIFY WHICH ANIMAL AND SESSION TO LOAD
meta.anm = 'JEB7'; % 'JEB7'  'EKH3'  'JGR2'
meta.date = '2021-04-29'; % '2021-04-29'  '2021-08-11'  '2021-11-16'


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
params.lfads_run = 'run3'; % 'run3' , leave empty to use most recent run
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

%%

lfads = dat.rates;
pop = obj.trialdat(:,:,dat.trials);

sm = 31;

% f = figure;
corrs = zeros(size(pop,2),1);
for i = 1:size(pop,2)
%     clf(f); hold on;
%     
%     plot(obj.time,mean(lfads(:,i,:),3))
%     plot(obj.time,mySmooth(mean(pop(:,i,:),3),51))
%     
%     pause; hold off;

    lfadspsth = mean(lfads(:,i,:),3);
    poppsth = mySmooth(mean(pop(:,i,:),3),sm);
    
    temp = corrcoef(lfadspsth,poppsth);
    corrs(i) = temp(1,2);

end

%%
close all
figure; 
nbins = round(size(pop,2) / 3);
h = histogram(corrs,nbins,'EdgeColor','none');
hold on;
xline(mean(corrs),'k','LineWidth',3)
xlabel('Correlation')
ylabel('Count')
ax = gca;
ax.FontSize = 20;

% create smaller axes in top right, and plot on it
axes('Position',[.2 .6 .3 .25])
box on
h2 = cdfplot(corrs);
h2.LineWidth = 2;
ax = h2.Parent;
h2.Parent.XLabel.String = '';
h2.Parent.YLabel.String = '';
h2.Parent.Title.String = 'CDF';
h2.Parent.FontSize = 15;










