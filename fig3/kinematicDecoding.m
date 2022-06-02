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

%% NEURAL ACTIVITY and KINEMATICS

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
    [kin{i}, featLeg] = getKinematicsFromVideo(obj(i),params(i),sort(dat(i).trials));
    gpfa(i) = getGPFAData(meta(i),'run1');
    disp(' ')
    disp('DONE')
end

%% Decode


feats2decode = {'tongue_ydisp_view1','jaw_ydisp_view1','nose_ydisp_view1'};
[~,mask] = patternMatchCellArray(featLeg,feats2decode,'any');
featix = find(mask);

winsize = 14;
kernsd = std(1:winsize)*params(1).dt;

for sessix = 1:numel(meta)
    
    lfadsdat = permute(dat(sessix).factors,[1,3,2]);
    lfadsdat = reshape(lfadsdat,size(lfadsdat,1)*size(lfadsdat,2),size(lfadsdat,3));
    lfadsdat = cat(2,ones(size(lfadsdat,1),1),lfadsdat);
    
    gpfadat = permute(gpfa(sessix).gpfalatents,[1,3,2]);
    gpfadat = reshape(gpfadat,size(gpfadat,1)*size(gpfadat,2),size(gpfadat,3));
    gpfadat = cat(2,ones(size(gpfadat,1),1),gpfadat);
    
    [~, gaussdat_fa] = getGaussData(obj(sessix),size(dat(sessix).factors,2),winsize);
    
    for i = 1:numel(featix)
        ix = featix(i);
        
        kindat = mySmooth(kin{sessix}(:,:,ix),11);
        
        [~,~,r2.lfads{sessix}(i)] = regressData(kindat,lfadsdat);
        [~,~,r2.gpfa{sessix}(i)] = regressData(kindat,gpfadat);
%         [~,~,r2.gausspca{sessix}(i)] = regressData(kindat,gaussdat_pca);
        [~,~,r2.gaussfa{sessix}(i)] = regressData(kindat,gaussdat_fa);

    end
end

%%
close all

cols = {[0 0.4470 0.9010],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880]};

f = figure; ax = axes(f); hold on;
for sessix = 1:numel(meta)
    for i = 1:numel(featix)
        scatter(r2.gaussfa{sessix}(i),r2.gpfa{sessix}(i),90,cols{i},'filled')
    end
end
refline([1 0])
xlabel('R^2 FA')
ylabel('R^2 GPFA')
xlim([0 1])
ylim([0 1])
refl = refline([1 0]);
refl.Color = [0 0 0];
refl.LineStyle = '--';
refl.LineWidth = 2;
title(['Gaussian Sigma = ' num2str(kernsd*1000) ' ms'])
ax.FontSize = 25;
legend({'Tongue','Jaw','Nose'},'Location','southeast','FontSize',15);

% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig3/figs/kindecoding';
% fn = ['run1_gpfa_fa_winsize_' num2str(winsize)];
% mysavefig(f,pth,fn);

f = figure; ax = axes(f); hold on;
for sessix = 1:numel(meta)
    for i = 1:numel(featix)
        scatter(r2.gpfa{sessix}(i),r2.lfads{sessix}(i),90,cols{i},'filled')
    end
end
refline([1 0])
xlabel('R^2 GPFA')
ylabel('R^2 LFADS')
xlim([0 1])
ylim([0 1])
refl = refline([1 0]);
refl.Color = [0 0 0];
refl.LineStyle = '--';
refl.LineWidth = 2;
title(['Gaussian Sigma = ' num2str(kernsd*1000) ' ms'])
ax.FontSize = 25;
legend({'Tongue','Jaw','Nose'},'Location','southeast','FontSize',15);

% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig3/figs/kindecoding';
% fn = ['run1_gpfa_lfads_winsize_' num2str(winsize)];
% mysavefig(f,pth,fn);

f = figure; ax = axes(f); hold on;
for sessix = 1:numel(meta)
    for i = 1:numel(featix)
        scatter(r2.gaussfa{sessix}(i),r2.lfads{sessix}(i),90,cols{i},'filled')
    end
end
xlabel('R^2 FA')
ylabel('R^2 LFADS')
xlim([0 1])
ylim([0 1])
refl = refline([1 0]);
refl.Color = [0 0 0];
refl.LineStyle = '--';
refl.LineWidth = 2;
title(['Gaussian Sigma = ' num2str(kernsd*1000) ' ms'])
ax.FontSize = 25;
legend({'Tongue','Jaw','Nose'},'Location','southeast','FontSize',15);


% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig3/figs/kindecoding';
% fn = ['run1_fa_lfads_winsize_' num2str(winsize)];
% mysavefig(f,pth,fn);



%% Helper Functions

function [gaussdat_pca, gaussdat_fa] = getGaussData(obj,nFactors,winsize)

gaussdat = permute(obj.trialdat,[1,3,2]);
nTime = size(gaussdat,1);
nTrials = size(gaussdat,2);
nClu = size(gaussdat,3);
gaussdat = reshape(gaussdat,nTime*nTrials,nClu);

% pca
[~,gaussdat_pca] = pca(gaussdat,'NumComponents',nFactors);
% fa
[~,~,~,~,gaussdat_fa] = factoran(gaussdat,nFactors);

gaussdat_pca = reshape(gaussdat_pca,nTime,nTrials,nFactors);
gaussdat_fa = reshape(gaussdat_fa,nTime,nTrials,nFactors);

N = winsize;
for i = 1:nFactors
    gaussdat_pca(:,:,i) = mySmooth(gaussdat_pca(:,:,i),N);
    gaussdat_fa(:,:,i) = mySmooth(gaussdat_fa(:,:,i),N);
end

gaussdat_pca = reshape(gaussdat_pca,nTime*nTrials,nFactors);
gaussdat_fa = reshape(gaussdat_fa,nTime*nTrials,nFactors);

gaussdat_pca = cat(2,ones(size(gaussdat_pca,1),1),gaussdat_pca); % add a column of ones
gaussdat_fa = cat(2,ones(size(gaussdat_fa,1),1),gaussdat_fa); % add a column of ones

end


function [mode,latent,r2] = regressData(kindat,neuraldat)


[mode,~,~,~,stats] = regress(kindat(:),neuraldat);
mode = mode./sum(abs(mode));

latent = neuraldat * mode;
latent = reshape(latent,size(kindat,1),size(kindat,2),size(latent,3));

r2 = stats(1);

end












