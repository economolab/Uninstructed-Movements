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
    [meta(i),params(i),obj(i),dat(i)] = getNeuralActivity(meta(i),dfparams);
    gpfa(i) = getGPFAData(meta(i),'run1');
    
    disp('DONE');
    disp(' ');
end




%% lfads factors, and gaussian smoothed pcs/factors

k = 5; % number of iterations (bootstrap)

dt = 4; % train/test a model every dt'th time point
mdlTime = obj(1).time(1:dt:numel(obj(1).time));
numT = numel(mdlTime);

acc_lfads = zeros(numT,k,numel(meta));
acc_gpfa = zeros(numT,k,numel(meta));
acc_gauss_pca = zeros(numT,k,numel(meta));
acc_gauss_fa = zeros(numT,k,numel(meta));

for sessix = 1:numel(meta)
    
    disp(['Session ' num2str(sessix) '/' num2str(numel(meta))])
    
    lfads = dat(sessix).factors;
    
    gpfa_ = gpfa(sessix).gpfalatents; % stored in descending trial order. flip 3rd dim
    
    gauss = obj(sessix).trialdat;
    
    % project gauss onto first nFactors(lfads) PCs using all data
    % using all trials b/c that's what we did with lfads as well
    nFactors = size(lfads,2);
    temp = permute(gauss,[1 3 2]);
    temp_reshaped = reshape(temp,size(temp,1)*size(temp,2),size(temp,3));
%     [pcs,pcalatents] = pca(temp_reshaped,'NumComponents',nFactors);
    [~,~,~,~,falatents] = factoran(temp_reshaped,nFactors);
    
%     pcalatents = reshape(pcalatents,size(temp,1),size(temp,2),nFactors);
%     gauss_pca = permute(pcalatents,[1 3 2]);
    
    falatents = reshape(falatents,size(temp,1),size(temp,2),nFactors);
    gauss_fa = permute(falatents,[1 3 2]);
    
    sm = 14;
    for i = 1:size(gauss_fa,2) % for each factor
%         gauss_pca(:,i,:) = mySmooth(squeeze(gauss_pca(:,i,:)),sm); % causal gaussian filter
        gauss_fa(:,i,:) = mySmooth(squeeze(gauss_fa(:,i,:)),sm); % causal gaussian filter
    end
    
    
    
    
    %% choice decoding
    

    
    train = 0.8;
    nTrials = min(numel(params(sessix).trialid{2}),numel(params(sessix).trialid{3}));
    nTrain = round(train*nTrials); % per condition
    nTest  = nTrials - nTrain;     % per condition
    
    for j = 1:k
        
        disp(['--Iteration ' num2str(j) '/' num2str(k)])
        
        rhit_train = randsample(params(sessix).trialid{2},nTrain);
        lhit_train = randsample(params(sessix).trialid{3},nTrain);
        
        rhitremain = params(sessix).trialid{2}(~ismember(params(sessix).trialid{2},rhit_train));
        lhitremain = params(sessix).trialid{3}(~ismember(params(sessix).trialid{3},lhit_train));
        
        rhit_test = randsample(rhitremain,nTest);
        lhit_test = randsample(lhitremain,nTest);
        
        y_train = cat(1,0*ones(numel(rhit_train),1),1*ones(numel(lhit_train),1));
        
        y_test = cat(1,0*ones(numel(rhit_test),1),1*ones(numel(lhit_test),1));
        
        X_train_lfads = lfads(:,:,[rhit_train;lhit_train]);
        X_test_lfads = lfads(:,:,[rhit_test;lhit_test]);
        
        X_train_gpfa = gpfa_(:,:,[rhit_train;lhit_train]);
        X_test_gpfa = gpfa_(:,:,[rhit_test;lhit_test]);
        
%         X_train_gauss_pca = gauss_pca(:,:,[rhit_train;lhit_train]);
%         X_test_gauss_pca = gauss_pca(:,:,[rhit_test;lhit_test]);
        
        X_train_gauss_fa = gauss_fa(:,:,[rhit_train;lhit_train]);
        X_test_gauss_fa = gauss_fa(:,:,[rhit_test;lhit_test]);
        
        
        acc_lfads = train_test_choiceDecoder(acc_lfads,X_train_lfads,X_test_lfads,...
            y_train,y_test,j,dt,sessix);
        
        acc_gpfa = train_test_choiceDecoder(acc_gpfa,X_train_gpfa,X_test_gpfa,...
            y_train,y_test,j,dt,sessix);
        
%         acc_gauss_pca = train_test_choiceDecoder(acc_gauss_pca,X_train_gauss_pca,X_test_gauss_pca,...
%             y_train,y_test,j,dt,sessix);
        
        acc_gauss_fa = train_test_choiceDecoder(acc_gauss_fa,X_train_gauss_fa,X_test_gauss_fa,...
            y_train,y_test,j,dt,sessix);
        
        
    end
    
end


%% accuracy

close all

align = mode(obj(1).bp.ev.(params(1).alignEvent));
sample = mode(obj(1).bp.ev.sample) - align;
delay = mode(obj(1).bp.ev.delay) - align;

f = figure(11); ax = axes(f); hold on;
stderr = std(std(acc_lfads,[],2),[],3) ./ (k);
means = mean(mean(acc_lfads,2),3);
shadedErrorBar(mdlTime, means, stderr, {'Color',[70, 176, 106]./255,'LineWidth',3},0.5, ax);

stderr = std(std(acc_gpfa,[],2),[],3) ./ (k);
means = mean(mean(acc_gpfa,2),3);
shadedErrorBar(mdlTime, means, stderr, {'Color',[50, 102, 168]./255,'LineWidth',3},0.5, ax);

% stderr = std(std(acc_gauss_pca,[],2),[],3) ./ (k*numel(meta));
% means = mean(mean(acc_gauss_pca,2),3);
% shadedErrorBar(mdlTime, means, stderr, {'Color',[115, 68, 130]./255,'LineWidth',3},0.5, ax);

stderr = std(std(acc_gauss_fa,[],2),[],3) ./ (k);
means = mean(mean(acc_gauss_fa,2),3);
shadedErrorBar(mdlTime, means, stderr, {'Color',[115, 68, 130]./255,'LineWidth',3},0.5, ax);

xline(sample,'k--','LineWidth',2)
xline(delay,'k--','LineWidth',2)
xline(0,'k--','LineWidth',2)

title([meta.anm ' ' meta.date ' | LFADS ' params.lfads_run ])
xlabel('Time (s) from go cue')
ylabel('Choice Decoder Accuracy')
% legend({'LFADS','PCA'})
ax = gca;
ax.FontSize = 20;
xlim([mdlTime(5) mdlTime(end)])
hold off

pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig3/figs/choiceDecoder';
fn = [meta.anm '_' meta.date '_' params.lfads_run '_' 'winsize_' num2str(sm)];
mysavefig(f,pth,fn)


%% Helper Functions

function acc = train_test_choiceDecoder(acc,X_train,X_test,y_train,y_test,bootiter,dt,sessix)
ix = 1;
for i = 1:dt:size(X_train,1) % each timepoint
    x_train = squeeze(X_train(i,:,:)); % (factors,trials) for timepoint i
    %         mdl = fitclinear(x_train',y_train);
    mdl = fitcsvm(x_train',y_train);
    x_test = squeeze(X_test(i,:,:)); % (factors,trials) for timepoint i
    pred = predict(mdl,x_test');
    acc(ix,bootiter,sessix) = sum(pred == y_test) / numel(y_test);
    ix = ix + 1;
end
end

