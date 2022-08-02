clear,clc,close all
addpath(genpath(pwd))

%% PARAMETERS

% --SPECIFY WHICH ANIMALs/SESSIONs TO LOAD
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

% % meta(end+1).anm = 'EKH3';
% % meta(end).date = '2021-08-07';

meta(end+1).anm = 'EKH1';
meta(end).date = '2021-08-07';

meta = assignDataPath(meta);

% ------------------------------------------------------------------ %
% get some parameters that are common to all sessions/analyses
dfparams = getDefaultParams();

%% NEURAL ACTIVITY

% data objects contain
% 1) processed obj
% 2) processing params
% 3) factor analysis results
% 4) params used for factor analysis
% 5) meta data that we don't really need

[obj,fa,params,me,kin,kinfeats,kinfeats_reduced] = loadProcessedData(meta);


%% covariance estimation
% some helpful analysis in https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0271136

sessix = 1;


% get fa or neural data (single trials)
covariates = [1 2];
fulldat = fa(sessix).falatents;
fulldat = permute(fulldat,[1 3 2]); % (time,trials,factors);
% fulldat = obj(sessix).trialdat;
% fulldat = permute(fulldat,[1 3 2]); % (time,trials,factors);

dat{1} = fulldat(:,params(sessix).trialid{2},covariates); % right trials
dat{2} = fulldat(:,params(sessix).trialid{3},covariates); % left trials

dat_reshape = cellfun(@(x) reshape(x,size(x,1)*size(x,2),size(x,3)), dat, 'UniformOutput',false);
cov_true = cellfun(@(x) cov(x), dat_reshape, 'UniformOutput',false);
truecov = cellfun(@(x) x(1,2), cov_true, 'UniformOutput',false);

nIters = 100; % number of iterations

for k = 1:nIters

    % randomly sort trials in dat
    sorted_dat = cell(size(dat));
    for i = 1:numel(dat)
        sortix = randsample(size(dat{i},2),size(dat{i},2));
        sorted_dat{i} = dat{i}(:,sortix,:);
    end

    % select 1:nTrials trials to estimate cov with, plot error for each
    % of 1:nTrials
    for i = 1:numel(sorted_dat)
        temp = sorted_dat{i};
        for j = 1:size(temp,2) % for each trial
            trialdat = temp(:,1:j,:);
            trialdat = reshape(trialdat,size(trialdat,1)*size(trialdat,2),size(trialdat,3));
            datcov = cov(trialdat);
            estcov{i}(j,k) = datcov(1,2);
            cov_error{i}(j,k) = mse(estcov{i}(j,k),truecov{i});
        end
    end


end

mean_error = cellfun(@(x)   mean(x,2)  , cov_error,'UniformOutput',false);
for i = 1:numel(dat)
    temp = cov_error{i};
    for k = 1:size(temp,1) % trials
        ci_low{i}(k) = prctile(temp(k,:),5);
        ci_high{i}(k) = prctile(temp(k,:),95);
    end
end


%%

clrs = getColors();
cols{1} = clrs.rhit;
cols{2} = clrs.lhit;
lw = 2;
alph = 0.5;

close all
f = figure; f.Position = [370         278        1163         670];
ax = axes(f); hold on;
for i = 1:numel(cov_error)
%     plot(1:numel(mean_error{i}),mean_error{i},'Color',cols{i},'LineWidth',lw)
    h(i) = errorbar(1:numel(mean_error{i}),mean_error{i},ci_low{i},ci_high{i}, 'Color',cols{i},'LineWidth',2);
end
legend({['True=' num2str(truecov{1})] , ['True=' num2str(truecov{2})]})
xlabel('# of trials')
ylabel('Error in Covariance Estimation')



    

% figure; 
% hold on
% plot(mean_error{1})
% plot(ci_low{1})
% plot(ci_high{1})














