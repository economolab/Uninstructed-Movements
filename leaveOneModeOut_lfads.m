clear,clc,close all

addpath(genpath(pwd))


% load lfads input and output
datapth = '/Users/Munib/Documents/Economo-Lab/code/lfads/input_data/';
datafn = 'JEB7_2021-04-29_run2.mat';
temp = load(fullfile(datapth,datafn));
datin = temp.lfads; clear temp;

% load lfads output
lfadsdatapth = '/Users/Munib/Documents/Economo-Lab/code/lfads/output';
lfadsfname = 'model_runs_h5_train_posterior_sample_and_average_JEB7-2021-04-29_run2';
datout = HDF5ToStruct(lfadsdatapth,lfadsfname);

% conditions
rhit = find(datin.obj.bp.R & datin.obj.bp.hit & ~datin.obj.bp.early & ~datin.obj.bp.autowater);
lhit = find(datin.obj.bp.L & datin.obj.bp.hit & ~datin.obj.bp.early & ~datin.obj.bp.autowater);
rmiss = find(datin.obj.bp.R & datin.obj.bp.miss & ~datin.obj.bp.early & ~datin.obj.bp.autowater);
lmiss = find(datin.obj.bp.L & datin.obj.bp.miss & ~datin.obj.bp.early & ~datin.obj.bp.autowater);
if contains(lfadsfname,'train')
    trials = datin.train_trials;
    
elseif contains(lfadsfname,'valid')
    trials = datin.valid_trials;
end
rhit = ismember(trials,rhit);
lhit = ismember(trials,lhit);
rmiss = ismember(trials,rmiss);
lmiss = ismember(trials,lmiss);


%% neural data

% lfads
rates_lfads = datout.output_dist_params;
factors_lfads = datout.factors;
rates_lfads = permute(rates_lfads,[2,1,3]);
factors_lfads = permute(factors_lfads,[2,1,3]);

% just smoothing
rates_smooth = datin.obj.trialdat(:,:,datin.train_trials);
for i = 1:size(rates_smooth,3)
    rates_smooth(:,:,i) = mySmooth(rates_smooth(:,:,i),51);
end

factors_smooth = zeros(size(factors_lfads));
for i = 1:size(factors_smooth,3)
    [~,factors_smooth(:,:,i)] = pca(rates_smooth(:,:,i),'numComponents',10);
end

% figure(11); clf(figure(11));
% ax(1) = subplot(121);
% imagesc(squeeze(factors_lfads(:,2,:))')
% ax(2) = subplot(122);
% imagesc(squeeze(factors_smooth(:,2,:))')
% linkaxes(ax,'xy');

%%

traj = datin.obj.traj{1};  % Get the video data

trials = datin.train_trials;
taxis = datin.obj.time;

jawPos = nan(numel(taxis),numel(trials));
jawVel = nan(numel(taxis),numel(trials));
jawAcc = nan(numel(taxis),numel(trials));
for i = 1:numel(trials)                        % For every trial in the condition
    trix = trials(i);
    if isnan(traj(trix).NdroppedFrames )                       % If the video data from this trial isn't good, skip it
        continue;
    end
    
    if ~isnan(traj(trix).frameTimes)                           % If the video data from this trial is good...
        ts = mySmooth(traj(trix).ts(:, 2, 2), 21);                                               % Side-view, up and down position of the jaw, smoothed
        jawPos(:,i) = interp1(traj(trix).frameTimes-0.5-(datin.obj.bp.ev.jawOnset(trix)), ts, taxis);          % Linear interpolation of jaw position to keep number of time points consistent across trials
        jawPos(:,i) = fillmissing(jawPos(:,i),'nearest');
        basederiv = median(myDiff(jawPos(:,i),400),'omitnan');                                         % Find the median jaw velocity (aka baseline)
    end
    %Find the difference between the jaw velocity and the
    %baseline jaw velocity
    jawVel(:, i) = abs(myDiff(jawPos(:,i),400)-basederiv);      % Values > 0 = jaw is moving
    % acc
    jawAcc(:,i) = abs(myDiff(jawVel(:,i),400));
end



% concatenate time and trials to make each 3d matrix a 2d matrix of size
% (time*trials,numFeats)

temp = rates_lfads;
temp = permute(temp,[1 3 2]);
temp = reshape(temp,size(temp,1)*size(temp,2),size(temp,3));
rates_lfads = temp;

temp = factors_lfads;
temp = permute(temp,[1 3 2]);
temp = reshape(temp,size(temp,1)*size(temp,2),size(temp,3));
factors_lfads = temp;

temp = rates_smooth;
temp = permute(temp,[1 3 2]);
temp = reshape(temp,size(temp,1)*size(temp,2),size(temp,3));
rates_smooth = temp;

temp = factors_smooth;
temp = permute(temp,[1 3 2]);
temp = reshape(temp,size(temp,1)*size(temp,2),size(temp,3));
factors_smooth = temp;

jawPos = jawPos(:);
jawVel = jawVel(:);
jawAcc = jawAcc(:);

save('jeb7_decoding_lfads',...
       'rates_lfads','factors_lfads','rates_smooth','factors_smooth',...
       'jawPos','jawVel','jawAcc')



