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
params.lfads_run = 'run12'; % 'run3' , leave empty to use most recent run
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



%% lfads and gaussian smoothed single trials (left and right hits only)

lfads = dat.rates;

gauss = obj.trialdat(:,:,dat.trials);

rhits = params.trialid{2};
lhits = params.trialid{3};
[~,mask] = ismember(dat.trials,rhits);
rhits = find(mask);
[~,mask] = ismember(dat.trials,lhits);
lhits = find(mask);

temp{1} = lfads(:,:,rhits);
temp{2} = lfads(:,:,lhits);

lfads = temp; clear temp

temp{1} = gauss(:,:,rhits);
temp{2} = gauss(:,:,lhits);

gauss = temp; clear temp


sm = 31; winsize = sm*params.dt;
for j = 1:2
    for i = 1:size(gauss{j},2) % for each clu
        gauss{j}(:,i,:) = mySmooth(squeeze(gauss{j}(:,i,:)),sm); % causal gaussian filter
    end
end

nFactors = size(lfads{1},2);
temp = cat(3,gauss{1},gauss{2});
temp2 = permute(temp,[1 3 2]);
temp_reshaped = reshape(temp2,size(temp2,1)*size(temp2,2),size(temp2,3));
[pcs,gausslatents] = pca(temp_reshaped,'NumComponents',nFactors);
% [~,~,~,~,gausslatents] = factoran(temp_reshaped,nFactors);




gausslatents = reshape(gausslatents,size(temp2,1),size(temp2,2),nFactors);
gausslatents = permute(gausslatents,[1 3 2]);

lfadslatents = cat(3,lfads{1},lfads{2});

%% choice decoding

k = 1; % number of iterations (bootstrap)

acc_lfads = zeros(size(lfadslatents,1),k);
acc_pca = zeros(size(gausslatents,1),k);

train = 0.8;
nTrials = min(numel(rhits),numel(lhits));
nTrain = round(train*nTrials); % per condition
nTest  = nTrials - nTrain;     % per condition

for j = 1:k
    
    disp(['Iteration ' num2str(j) '/' num2str(k)])
    
    r = 1:numel(rhits);
    l = (numel(rhits)+1):(numel(rhits) + numel(lhits));
    
    rhit_train = randsample(r,nTrain);
    lhit_train = randsample(l,nTrain);
    
    rhitremain = find(~ismember(r,rhit_train));
    lhitremain = find(~ismember(l,lhit_train));
    
    rhit_test = randsample(rhitremain,nTest);
    lhit_test = randsample(lhitremain,nTest);
    
    y_train = cat(1,0*ones(numel(rhit_train),1),1*ones(numel(lhit_train),1));
    
    y_test = cat(1,0*ones(numel(rhit_test),1),1*ones(numel(lhit_test),1));
    

    
    %%
    
    % LFADS factors
    
    X_train = lfadslatents(:,:,[rhit_train;lhit_train]);
    X_test = lfadslatents(:,:,[rhit_test;lhit_test]);

    for i = 1:size(X_train,1) % each timepoint
        x_train = squeeze(X_train(i,:,:)); % (factors,trials) for timepoint i
%         mdl = fitclinear(x_train',y_train);
        mdl = fitcsvm(x_train',y_train);
        x_test = squeeze(X_test(i,:,:)); % (factors,trials) for timepoint i
        pred = predict(mdl,x_test');
%         unique(pred)
        acc_lfads(i,j) = sum(pred == y_test) / numel(y_test);
    end
    
    % PCA factors
    
    
    X_train = gausslatents(:,:,[rhit_train;lhit_train]);
    X_test = gausslatents(:,:,[rhit_test;lhit_test]);
    
    
    for i = 1:size(X_train,1) % each timepoint
        x_train = squeeze(X_train(i,:,:)); % (factors,trials) for timepoint i
        mdl = fitcsvm(x_train',y_train);
        x_test = squeeze(X_test(i,:,:)); % (factors,trials) for timepoint i
        pred = predict(mdl,x_test');
        acc_pca(i,j) = sum(pred == y_test) / numel(y_test);
    end
    
    
    lfads_r = lfadslatents(:,:,rhit_train);
    lfads_l = lfadslatents(:,:,lhit_train);
    dist_r = zeros(size(lfads_r));
    dist_l = zeros(size(lfads_l));
    for i = 1:size(X_train_lfads,1) % each timepoint
        dist_r(i,:,:) = squeeze(lfads_r(i,:,:));
        dist_l(i,:,:) = squeeze(lfads_l(i,:,:));
        
%         histogram(r(:),'FaceColor','b'); hold on;
%         histogram(l(:),'FaceColor','r');
%         pause
%         hold off
    end
    test = dist_r - dist_l;
    
    
end


%% accuracy
align = mode(obj.bp.ev.(params.alignEvent));
sample = mode(obj.bp.ev.sample) - align;
delay = mode(obj.bp.ev.delay) - align;

f = figure(11); ax = axes(f); hold on;

stderr = mySmooth(std(acc_lfads,[],2) ./ k, 31);
means = mySmooth(mean(acc_lfads,2) , 31);
shadedErrorBar(obj.time, means, stderr, {'Color',[70, 176, 106]./255,'LineWidth',1.5},0.5, ax);

stderr = mySmooth(std(acc_pca,[],2) ./ k, 31);
means = mySmooth(mean(acc_pca,2) , 31);
shadedErrorBar(obj.time, means, stderr, {'Color',[115, 68, 130]./255,'LineWidth',1.5},0.5, ax);

% plot(obj.time,mean(acc_lfads,2),'Color',[70, 176, 106]./255,'LineWidth',3)
% plot(obj.time,mean(acc_pca,2),'Color',[115, 68, 130]./255,'LineWidth',3)

xline(sample,'k--','LineWidth',2)
xline(delay,'k--','LineWidth',2)
xline(0,'k--','LineWidth',2)

title([meta.anm ' ' meta.date ' | LFADS ' params.lfads_run ])
xlabel('Time (s) from go cue')
ylabel('Accuracy')
% legend({'LFADS','PCA'})
ax = gca;
ax.FontSize = 20;
xlim([obj.time(30) obj.time(end)])
hold off

pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig3/figs/choiceDecoder';
fn = [meta.anm '_' meta.date '_' params.lfads_run ];
mysavefig(f,pth,fn)



