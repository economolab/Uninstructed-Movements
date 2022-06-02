clear,clc,close all


addpath(genpath(pwd))

%% PARAMETERS

% --SPECIFY WHICH ANIMAL AND SESSION TO LOAD
meta(1).anm = 'JEB7'; % 'JEB7'  'EKH3'  'JGR2'
meta(1).date = '2021-04-29'; % '2021-04-29'  '2021-08-11'  '2021-11-16'

meta(end+1).anm = 'JEB6'; % 'JEB7'  'EKH3'  'JGR2'
meta(end).date = '2021-04-18'; % '2021-04-29'  '2021-08-11'  '2021-11-16'

meta(end+1).anm = 'JGR2'; % 'JEB7'  'EKH3'  'JGR2'
meta(end).date = '2021-11-17'; % '2021-04-29'  '2021-08-11'  '2021-11-16'

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
dfparams.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance
dfparams.full_or_reduced = 'reduced'; % 'full'  or 'reduced' -- which data to use in regression
% using the full data will require another method, the system seems to be
% highly overfit as is
assert(strcmpi(dfparams.full_or_reduced,'reduced'),'method to use full dimensional data doesnt exist, params.full_or_reduced should be set to `reduced`')

% --SPECIFY TIME POINTS AND LAG TO USE
dfparams.prep = [-2.5 -0.05]; % initial and final time points (seconds) defining prep epoch, relative to alignevent
dfparams.move = [-2.5 1.5];   % initial and final time points (seconds) defining move epoch, relative to alignevent
dfparams.advance_movement = 0.025; % seconds, amount of time to advance movement data relative to neural data

%% NEURAL ACTIVITY

% getNeuralActivity() returns 4 main variables
% - dat: contains lfads/fa smoothed firing rates, factors, and trial numbers
%        dat.factors and dat.rates are size (time,factors/clusters,trials)
% - meta: session meta data
% - params: parameters used for preprocessing lfads input data
% - obj: preprocessed data obj

for i = 1:numel(meta)
    disp(['Loading data for ' meta(i).anm ' ' meta(i).date]);
%     [meta,params,obj,dat] = getNeuralActivity(meta(i),params(i));
    [meta(i),params(i),obj(i),dat(i)] = getNeuralActivity(meta(i),dfparams);
end


%% correlate lfads psths and gaussian smoothed psths

corrs = cell(numel(meta),1);

for sessix = 1:numel(meta)
    
    lfads = dat(sessix).rates;
    
    gauss = obj(sessix).trialdat;
    
    N = 20;
    for i = 1:size(gauss,2) % for each clu
        gauss(:,i,:) = mySmooth(squeeze(gauss(:,i,:)),N); % causal gaussian filter
    end
    
    
    corrs{sessix} = zeros(size(gauss,2),1);
    for i = 1:size(gauss,2)
        
        
        lfadspsth = mean(lfads(:,i,:),3);
        poppsth = mySmooth(mean(gauss(:,i,:),3),N);
        
        temp = corrcoef(lfadspsth,poppsth);
        corrs{sessix}(i) = temp(1,2);
        
    end
    
end

corrs = cell2mat(corrs);


%%

close all
f = figure; 
nbins = round(numel(corrs) / 3);
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

% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig3/figs/lfads_psth_corr';
% fn = [meta.anm '_' meta.date '_' params.lfads_run ];
% mysavefig(f,pth,fn)


%% neural reconstruction

% To compare our ability to reconstruct single-trial responses using Gaussian smoothing and LFADS, 
% we first computed peri-event time histograms (PETHs) within condition using all training trials (excluding one test trial).
% We then computed the correlation between the firing rates of each test trial (smoothed with a Gaussian kernel or 
% reconstructed with LFADS) with the PETH of the corresponding condition averaged across the training trails 
% (Figure 2—figure supplement 1A). We repeated this procedure with a different trial left out for each condition.
% We report the difference in correlation coefficient obtained after LFADS processing and Gaussian smoothing (Figure 2—figure supplement 1B).


lfads_recon = cell(numel(meta),1);
gauss_recon = cell(numel(meta),1);
for sessix = 1:numel(meta)
    
    lfads = dat(sessix).rates;
    
    gauss = obj(sessix).trialdat;
    
    rhits = params(sessix).trialid{2};
    lhits = params(sessix).trialid{3};
    
    trials = [rhits;lhits]; % using right and left hit trials only
    
    lfads = lfads(:,:,trials);
    
    gauss = gauss(:,:,trials);
    
    sm = 31; sigma = 31;
    for i = 1:size(gauss,2) % for each clu
        gauss(:,i,:) = mySmooth(squeeze(gauss(:,i,:)),N); % causal gaussian filter
    end
    
    
    
    lfads_recon{sessix} = nan(numel(trials),1);
    gauss_recon{sessix} = nan(numel(trials),1);
    for i = 1:numel(trials)
        test = i;
        train = trials~=trials(i);
        
        lfads_train = lfads(:,:,train);
        lfads_test = lfads(:,:,test);
        
        gauss_train = gauss(:,:,train);
        gauss_test = gauss(:,:,test);
        
        psth = mean(gauss_train,3);
        
        lfads_recon{sessix}(i) = corr(lfads_test(:),psth(:));
        
        gauss_recon{sessix}(i) = corr(gauss_test(:),psth(:));
        
    end
    
end

lfads_recon = cell2mat(lfads_recon);
gauss_recon = cell2mat(gauss_recon);

%%

f = figure;
cols = [130, 176, 120 ; 106, 115, 171] ./ 255;
vs = violinplot([lfads_recon, gauss_recon],{'LFADS','Gaussian 30 ms'},...
                'ViolinColor',cols,'EdgeColor',[1 1 1]);
ylabel('Neural Reconstruction, CC')
ylim([0,1])
ax = gca;
ax.FontSize = 20;

% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig3/figs/reconstruction';
% fn = [meta.anm '_' meta.date '_' params.lfads_run ];
% mysavefig(f,pth,fn)











