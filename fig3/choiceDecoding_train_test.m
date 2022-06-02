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
params.lfads_run = 'run1'; % 'run3' , leave empty to use most recent run
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



%% lfads factors, and gaussian smoothed pcs/factors

lfads = dat.factors;

gauss = obj.trialdat;

% project gauss onto first nFactors(lfads) PCs using all data
% using all trials b/c that's what we did with lfads as well
nFactors = size(lfads,2);
temp = permute(gauss,[1 3 2]);
temp_reshaped = reshape(temp,size(temp,1)*size(temp,2),size(temp,3));
[pcs,pcalatents] = pca(temp_reshaped,'NumComponents',nFactors);
[~,~,~,~,falatents] = factoran(temp_reshaped,nFactors);

pcalatents = reshape(pcalatents,size(temp,1),size(temp,2),nFactors);
gauss_pca = permute(pcalatents,[1 3 2]);

falatents = reshape(falatents,size(temp,1),size(temp,2),nFactors);
gauss_fa = permute(falatents,[1 3 2]);

sm = 31; sigma = 31;
for i = 1:size(gauss_pca,2) % for each factor
    gauss_pca(:,i,:) = mySmooth(squeeze(gauss_pca(:,i,:)),sm,sigma); % causal gaussian filter
end


%% choice decoding

k = 10; % number of iterations (bootstrap)

dt = 10; % train/test a model every dt'th time point
mdlTime = obj.time(1:dt:numel(obj.time));
numT = numel(mdlTime);

acc_lfads = zeros(numT,k);
acc_gauss_pca = zeros(numT,k);
acc_gauss_fa = zeros(numT,k);

train = 0.8;
nTrials = min(numel(params.trialid{2}),numel(params.trialid{3}));
nTrain = round(train*nTrials); % per condition
nTest  = nTrials - nTrain;     % per condition

for j = 1:k
    
    disp(['Iteration ' num2str(j) '/' num2str(k)])
    
    rhit_train = randsample(params.trialid{2},nTrain);
    lhit_train = randsample(params.trialid{3},nTrain);

    rhitremain = params.trialid{2}(~ismember(params.trialid{2},rhit_train));
    lhitremain = params.trialid{3}(~ismember(params.trialid{3},lhit_train));
    
    rhit_test = randsample(rhitremain,nTest);
    lhit_test = randsample(lhitremain,nTest);
    
    y_train = cat(1,0*ones(numel(rhit_train),1),1*ones(numel(lhit_train),1));
    
    y_test = cat(1,0*ones(numel(rhit_test),1),1*ones(numel(lhit_test),1));
    
    X_train_lfads = lfads(:,:,[rhit_train;lhit_train]);
    X_test_lfads = lfads(:,:,[rhit_test;lhit_test]);
    
    X_train_gauss_pca = gauss_pca(:,:,[rhit_train;lhit_train]);
    X_test_gauss_pca = gauss_pca(:,:,[rhit_test;lhit_test]);
    
    X_train_gauss_fa = gauss_fa(:,:,[rhit_train;lhit_train]);
    X_test_gauss_fa = gauss_fa(:,:,[rhit_test;lhit_test]);
    

    acc_lfads = train_test_choiceDecoder(acc_lfads,X_train_lfads,X_test_lfads,...
                                         y_train,y_test,j,dt);
                                     
    acc_gauss_pca = train_test_choiceDecoder(acc_gauss_pca,X_train_gauss_pca,X_test_gauss_pca,...
                                         y_train,y_test,j,dt);              
                                     
    acc_gauss_fa = train_test_choiceDecoder(acc_gauss_fa,X_train_gauss_fa,X_test_gauss_fa,...
                                         y_train,y_test,j,dt);      
       
    
end


%% accuracy

close all

align = mode(obj.bp.ev.(params.alignEvent));
sample = mode(obj.bp.ev.sample) - align;
delay = mode(obj.bp.ev.delay) - align;

f = figure(11); ax = axes(f); hold on;
stderr = std(acc_lfads,[],2) ./ k;
means = mean(acc_lfads,2);
shadedErrorBar(mdlTime, means, stderr, {'Color',[70, 176, 106]./255,'LineWidth',3},0.5, ax);

stderr = std(acc_gauss_pca,[],2) ./ k;
means = mean(acc_gauss_pca,2);
shadedErrorBar(mdlTime, means, stderr, {'Color',[115, 68, 130]./255,'LineWidth',3},0.5, ax);

stderr = std(acc_gauss_fa,[],2) ./ k;
means = mean(acc_gauss_fa,2);
shadedErrorBar(mdlTime, means, stderr, {'Color',[222, 111, 142]./255,'LineWidth',3},0.5, ax);

xline(sample,'k--','LineWidth',2)
xline(delay,'k--','LineWidth',2)
xline(0,'k--','LineWidth',2)

title([meta.anm ' ' meta.date ' | LFADS ' params.lfads_run ])
xlabel('Time (s) from go cue')
ylabel('Accuracy')
% legend({'LFADS','PCA'})
ax = gca;
ax.FontSize = 20;
xlim([mdlTime(5) mdlTime(end)])
hold off

% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig3/figs/choiceDecoder';
% fn = [meta.anm '_' meta.date '_' params.lfads_run ];
% mysavefig(f,pth,fn)


%% Helper Functions

function acc = train_test_choiceDecoder(acc,X_train,X_test,y_train,y_test,bootiter,dt)
ix = 1;
for i = 1:dt:size(X_train,1) % each timepoint
    x_train = squeeze(X_train(i,:,:)); % (factors,trials) for timepoint i
    %         mdl = fitclinear(x_train',y_train);
    mdl = fitcsvm(x_train',y_train);
    x_test = squeeze(X_test(i,:,:)); % (factors,trials) for timepoint i
    pred = predict(mdl,x_test');
    acc(ix,bootiter) = sum(pred == y_test) / numel(y_test);
    ix = ix + 1;
end
end

