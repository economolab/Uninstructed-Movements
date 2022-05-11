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



%% conditions

rhit = find(obj.bp.R & obj.bp.hit & ~obj.bp.early & ~obj.bp.autowater);
lhit = find(obj.bp.L & obj.bp.hit & ~obj.bp.early & ~obj.bp.autowater);
rmiss = find(obj.bp.R & obj.bp.miss & ~obj.bp.early & ~obj.bp.autowater);
lmiss = find(obj.bp.L & obj.bp.miss & ~obj.bp.early & ~obj.bp.autowater);

rhit_lfads = ismember(dat.trials,rhit);
lhit_lfads = ismember(dat.trials,lhit);
rmiss = ismember(dat.trials,rmiss);
lmiss = ismember(dat.trials,lmiss);

rhit = find(rhit_lfads);
lhit = find(lhit_lfads);


%% choice decoding

nTrials = min(numel(rhit),numel(lhit));


k = 10;

for j = 1:k
    
    disp(['Bootstrap Iteration (' num2str(j) '/' num2str(k) ')']);
    
    trials.rlfads = randsample(rhit,nTrials);
    trials.llfads = randsample(lhit,nTrials);
    trials.rpop = dat.trials(trials.rlfads);
    trials.lpop = dat.trials(trials.llfads);
    
    trials.nums = cat(1,trials.rlfads,trials.llfads);
    
    X = dat.factors(:,:,trials.nums);
    
    y = cat(1,ones(nTrials,1),2*ones(nTrials,1));
    
    %     y = zeros(size(rhit));
    %     y(rhit) = 1;
    %     y(lhit) = 2;
    %     ix = find(y==0);
    %     y(ix) = [];
    %     X(:,:,trials.nums) = [];
    
    ce = zeros(size(X,1),1);
    costs_lfads = zeros(size(X,1),k);
    acc_lfads = costs_pca;
    for i = 1:size(X,1) % each timepoint
        x = squeeze(X(i,:,:)); % (factors,trials) for timepoint i
        [Mdl_lfads,FitInfo] = fitclinear(x',y);
        %         [Mdl_lfads] = fitclinear(x',y,'KFold',5,...
        %             'Learner','logistic');
        %         ce(i) = kfoldLoss(Mdl_lfads);
        costs_lfads(i,j) = FitInfo.Objective;
        acc_lfads(i,j) = sum(Mdl_pca.predict(x') == y) / numel(y);
    end
    
    
    
    
    % compare against smoothed single trials projected onto first nFactors PCs
    rdat = obj.trialdat(:,:,trials.rpop);
    ldat = obj.trialdat(:,:,trials.lpop);
    
    temp = cat(3,rdat,ldat);
    temp = permute(temp,[1 3 2]);
    temp = reshape(temp,size(temp,1)*size(temp,2),size(temp,3));
    pcs = pca(temp,'NumComponents',size(X,2));
    
    
    pclatents = zeros(size(X,1),size(X,2),numel(dat.trials));
    for i = 1:numel(dat.trials)
        trialnum = dat.trials(i);
        temp = obj.trialdat(:,:,trialnum);
        pclatents(:,:,i) = temp * pcs;
    end
    
    
    X = pclatents(:,:,trials.nums);
    
    y = cat(1,ones(nTrials,1),2*ones(nTrials,1));
    
    
    %     y = zeros(size(rhit));
    %     y(rhit) = 1;
    %     y(lhit) = 2;
    %     ix = find(y==0);
    %     y(ix) = [];
    %     X(:,:,ix) = [];
    
    ce = zeros(size(X,1),1);
    costs_pca = zeros(size(X,1),k);
    acc_pca = costs_pca;
    for i = 1:size(X,1) % each timepoint
        x = squeeze(X(i,:,:)); % (factors,trials) for timepoint i
        [Mdl_pca,FitInfo] = fitclinear(x',y);
        %         [Mdl_pca] = fitclinear(x',y,'KFold',5,...
        %             'Learner','logistic');
        %         ce(i) = kfoldLoss(Mdl_pca);
        acc_pca(i,j) = sum(Mdl_pca.predict(x') == y) / numel(y);
        costs_pca(i,j) = FitInfo.Objective;
    end
    
end


%% costs
align = mode(obj.bp.ev.(params.alignEvent));
sample = mode(obj.bp.ev.sample) - align;
delay = mode(obj.bp.ev.delay) - align;

figure(10); clf; hold on;

m = mean(costs_pca,2);
CI(:,1) = quantile(costs_pca',0.025);
CI(:,2) = quantile(costs_pca',0.975);
% stdmean = std(costs_pca,[],2);
% CI = m-stdmean;
% CI(CI<0) = 0;
% CI = [CI,m+stdmean];

col = [115, 68, 130]./255;
col2 = [175, 130, 189]./255;
plot_ci(obj.time,[m,CI(:,1),CI(:,2)],'PatchColor', col, 'PatchAlpha', 0.1, ...
    'MainLineWidth', 2, 'MainLineStyle', '-', 'MainLineColor', col, ...
    'LineWidth', 0.5, 'LineStyle','-', 'LineColor', col2);

m = mean(costs_lfads,2);
CI(:,1) = quantile(costs_lfads',0.025);
CI(:,2) = quantile(costs_lfads',0.975);
% stdmean = std(costs_lfads,[],2);
% CI = m-stdmean;
% CI(CI<0) = 0;
% CI = [CI,m+stdmean];
col = [70, 176, 106]./255;
col2 = [104, 186, 132]./255;
plot_ci(obj.time,[m,CI(:,1),CI(:,2)],'PatchColor', col, 'PatchAlpha', 0.1, ...
    'MainLineWidth', 3, 'MainLineStyle', '-', 'MainLineColor', col, ...
    'LineWidth', 0.5, 'LineStyle','-', 'LineColor', col2);


xline(sample,'k--','LineWidth',2)
xline(delay,'k--','LineWidth',2)
xline(0,'k--','LineWidth',2)

title([meta.anm ' ' meta.date ' | LFADS ' params.lfads_run ])
xlabel('Time (s) from go cue')
ylabel('Error')
% legend({'LFADS','PCA'},'Location','best')
ax = gca;
ax.FontSize = 20;
xlim([obj.time(30) obj.time(end)])
hold off

%% accuracy

align = mode(obj.bp.ev.(params.alignEvent));
sample = mode(obj.bp.ev.sample) - align;
delay = mode(obj.bp.ev.delay) - align;

figure(11); clf; hold on;

m = mean(acc_pca,2);
CI = quantile(acc_pca',[0.025 0.975])';
% CI(:,2) = quantile(acc_pca',0.975);
% stdmean = std(acc_pca,[],2);
% CI = m-stdmean;
% CI(CI<0) = 0;
% CI = [CI,m+stdmean];

col = [115, 68, 130]./255;
col2 = [175, 130, 189]./255;
plot_ci(obj.time,[m,CI(:,1),CI(:,2)],'PatchColor', col, 'PatchAlpha', 0.1, ...
    'MainLineWidth', 2, 'MainLineStyle', '-', 'MainLineColor', col, ...
    'LineWidth', 0.5, 'LineStyle','-', 'LineColor', col2);

m = mean(acc_lfads,2);
CI(:,1) = quantile(acc_lfads',0.025);
CI(:,2) = quantile(acc_lfads',0.975);
% stdmean = std(acc_lfads,[],2);
% CI = m-stdmean;
% CI(CI<0) = 0;
% CI = [CI,m+stdmean];
col = [70, 176, 106]./255;
col2 = [104, 186, 132]./255;
plot_ci(obj.time,[m,CI(:,1),CI(:,2)],'PatchColor', col, 'PatchAlpha', 0.1, ...
    'MainLineWidth', 3, 'MainLineStyle', '-', 'MainLineColor', col, ...
    'LineWidth', 0.5, 'LineStyle','-', 'LineColor', col2);


xline(sample,'k--','LineWidth',2)
xline(delay,'k--','LineWidth',2)
xline(0,'k--','LineWidth',2)

title([meta.anm ' ' meta.date ' | LFADS ' params.lfads_run ])
xlabel('Time (s) from go cue')
ylabel('Accuracy')
% legend({'LFADS','PCA'},'Location','best')
ax = gca;
ax.FontSize = 20;
xlim([obj.time(30) obj.time(end)])
hold off









