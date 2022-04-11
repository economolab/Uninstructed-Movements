clear,clc,close all

% finds null/potent spaces for a single session based on the linear equation

%               V = NW, 

% where V is a video feature matrix and N is a matrix of 
% single trial neural activity. W is the transformation
% matrix between video features and neural activity . null(W) is a basis for
% the null space of W and colspace(W) is a basis for the row space of W (which
% is the potent space).

% * (all functions associated with this script are in the kinModes directory) *

addpath(genpath(pwd))

%% PARAMETERS

% --SPECIFY WHICH ANIMAL AND SESSION TO LOAD
meta.anm = 'JEB7'; % 'JEB7'  'EKH3'  'JGR2'
meta.date = '2021-04-29'; % '2021-04-29'  '2021-08-11'  '2021-11-16'


% --SPECIFY PATH TO DATA HERE
meta.datapth = '/Users/Munib/Documents/Economo-Lab/data/';
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
params.feat_varToExplain = 90; % num factors for dim reduction of video features should explain this much variance
params.full_or_reduced = 'reduced'; % 'full'  or 'reduced' -- which data to use in regression
% using the full data will require another method, the system seems to be
% highly overfit as is

% % TODO: parameters for time points to use for regression (and a lag b/w neural activity
% % and video features)
% params.time_to_use = 0;
% params.lag = 0;

%% NEURAL ACTIVITY

% getNeuralActivity() returns 4 main variables
% - dat: contains lfads/fa smoothed firing rates, factors, and trial numbers
%        dat.factors and dat.rates are size (time,factors/clusters,trials)
% - meta: session meta data
% - params: parameters used for preprocessing lfads input data
% - obj: preprocessed data obj
[meta,params,obj,dat] = getNeuralActivity(meta,params);

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

%% VIDEO FEATURES

% getKinematicsFromVideo() returns 2 variables
% - kin: feature matrix of size (time,features,trials). Features defined in
%         params.traj_features
% - featLeg: legend corresponding to features in kin (for 2nd dimension)
[kin,dat.featLeg] = getKinematicsFromVideo(obj,params,dat.trials);


% JAW ANGLE
jawAngle = getJawAngle(obj.time, obj, dat.trials, params.alignEvent); % (time, trials)
dat.featLeg{end+1} = 'jaw_angle';

% create feature matrix, feats, and assign to a field in dat
dat.feats = cat(2,kin,reshape(jawAngle,size(jawAngle,1),1,size(jawAngle,2)));


% MOTION ENERGY
% me is a struct with fields
% - data:       cell array (ntrials,1), with motion energy for each time point
% - moveThresh: a scalar - below this value of motion energy, animal is
%               said to be stationary
% - use:        0 or 1, indicating whether to use me.data. 1 if
%               a motionEnergy*.mat file is found, 0 if file not found
me = loadMotionEnergy(obj,meta,params,dat.trials); 
if me.use
    dat.feats = cat(2,dat.feats,reshape(me.data,size(me.data,1),1,size(me.data,2))); 
    dat.featLeg{end+1} = 'motion_energy';
end
% To generate a motionEnergy*.mat file for a session, see https://github.com/economolab/videoAnalysisScripts/blob/main/motionEnergy.m


% DIMENSIONALITY REDUCTION
% many of the video features will be highly correlated, so we will perform PCA/FA
% on the matrix of features to reduce the dimensionality to a set of factors that
% best explain the movement captured by the video recordings
dat.feats_reduced = reduceDimensionVideoFeatures(dat.feats,params.feat_varToExplain);

disp('DONE CREATING FEATURE MATRIX AND REDUCED DIM FEATURE MATRIX')

% for i = 1:size(dat.feats_reduced,2)
%     figure(i)
% %     plot(obj.time,squeeze(dat.feats_reduced(:,i,1:numel(params.trialid{1}))))
% %     imagesc(obj.time,dat.trials,squeeze(dat.feats(:,i,:))'); title(dat.featLeg{i},'Interpreter','none')
%     imagesc(obj.time,dat.trials,squeeze(dat.feats_reduced(:,i,:))')
% %     colorbar
% %     caxis([-500 200])
% end

%% 

% We now have a struct, dat, that contains:
% - trials:        trial numbers in use
% - factors:       neural data that's been reduced using lfads or factor analysis (time,factors,trials)
% - rates:         denoised single trial neural data from lfads or smoothing (time,cells,trials)
% - featLeg:       cell array of strings describing the features in dat.feats (2nd dim)
% - feats:         displacements, velocity, jaw angle, and motion energy (time,feature,trials)
% - feats_reduced: pca reduced dat.feats (time,factors,trials);


%% REGRESSION

% here we pose the problem
%           
%               V = NW
%
% where V is a matrix of size (time*trials , numVideoFeatures)
% and   N is a matrix of size (time*trials , numFactors/numClusters)
% 
% We will estimate W with ridge regression (regularized linear regression).
% It is the linear transformation that takes neural activity to video
% activity
% The null space of W forms subspace that we will call the         NULL SPACE   of N
% The column space of W forms the subspace that we will call the   POTENT SPACE of N

[W,N,V] = estimateW(dat,params); % N,V are zscored neural activity and feature matrix



% gprMdl = fitrgp(N(1:5000,:),V(1:5000,1));
% V_ = resubPredict(gprMdl);
% 
% V_=predict(gprMdl,N(5001:10000,:));
% L = loss(gprMdl,N(5001:10000,:),V(5001:10000,1))



% V_ = N*W;
% ft = 1;
% figure; plot(V(1:1000,ft)); hold on; plot(V_(1:1000,ft))

% TODO:
% from kaufman2014:
% To accommodate the known lag between motor cortical activity and muscle activity,
% M was advanced by 50 ms relative to N (refs. 34,46); the movement period for the muscle data 
% therefore began at movement onset
% - implement a lag b/w N and V
% - split trials into 2 or 3 sets, train,test,valid 

%% NULL AND POTENT SPACE OF W

% TODO: W is (10,9) or (N,V). rank ends up being 10 or close to it (9
% sometimes). So 1 dimension for null space or no null space...
% will have to figure out good number of dims for factors and feats_reduced

% http://pillowlab.princeton.edu/teaching/statneuro2018/slides/notes03a_SVDandLinSys.pdf

% W_potent = colspace(W);
% W_null = null(W);

% find null and potent spaces of W (reference above)
[u,s,v] = svd(W);
W_potent = u(logical(diag(s)),:)'; % can also use licols to get col space of W or W'
W_null = u(~sum(s'),:)';

N_potent = N * W_potent;
N_null   = N * W_null;

N_potent = reshape(N_potent,size(dat.factors,1),size(dat.factors,3),size(N_potent,2));
N_potent = permute(N_potent,[1,3,2]);

N_null = reshape(N_null,size(dat.factors,1),size(dat.factors,3),size(N_null,2));
N_null = permute(N_null,[1,3,2]);
% 
% % % prep tuning
% gamma = norm((W_null * data.N(data.move_idx,:)'),'fro')^2 / ...
%     norm((W_potent * data.N(data.move_idx,:)'),'fro')^2;
% 
% N_null = (W_null * data.N(data.prep_idx,:)')';
% N_potent = (W_potent * data.N(data.prep_idx,:)')';
% 
% tuning_ratio = (1/gamma) * (norm(N_null,'fro')^2 / norm(N_potent,'fro')^2)


%% dPCA


% [~,mask] = patternMatchCellArray(dat.featLeg,{'jaw','xvel','0'},'all');
% ix = find(mask);
% 
% temp = squeeze(dat.feats(:,ix,1:numel(params.trialid{1})));
% figure; imagesc(temp')
% figure; plot(nanmean(temp,2))
% 
% 
% 
% 
% figure; imagesc(squeeze(N_potent(:,1,:))')
% figure; imagesc(squeeze(N_null(:,1,:))')
% 
% % for fa data
% figure; 
% plot(obj.time,squeeze(N_potent(:,1,1:numel(params.trialid{1}))),'b'); hold on
% trix = (numel(params.trialid{1})+1):((numel(params.trialid{1})+1) + numel(params.trialid{2}));
% plot(obj.time,squeeze(N_potent(:,1,trix)),'r')

% for lfads data
rhit = find(obj.bp.R & obj.bp.hit & ~obj.bp.autowater & ~obj.bp.early);
lhit = find(obj.bp.L & obj.bp.hit & ~obj.bp.autowater & ~obj.bp.early);
mask = ismember(dat.trials,rhit);
rhit = find(mask);
mask = ismember(dat.trials,lhit);
lhit = find(mask);

figure; 
dim = 1;
plot(obj.time,squeeze(N_null(:,dim,rhit)),'b'); hold on;
plot(obj.time,squeeze(N_null(:,dim,lhit)),'r')

figure; 
dim = 1;
plot(obj.time,squeeze(N_potent(:,dim,rhit)),'b'); hold on;
plot(obj.time,squeeze(N_potent(:,dim,lhit)),'r')


nPlot = 20;
gc = round(numel(obj.time)/2);
figure; hold on
trix = randsample(rhit,nPlot);
plot3(squeeze(N_null(101:end,dim,trix)),squeeze(N_potent(101:end,2,trix)),squeeze(N_potent(101:end,1,trix)),'b')
plot3(squeeze(N_null(101,dim,trix)),squeeze(N_potent(101,2,trix)),squeeze(N_potent(101,1,trix)),'b.','MarkerSize',20)
plot3(squeeze(N_null(gc,dim,trix)),squeeze(N_potent(gc,2,trix)),squeeze(N_potent(gc,1,trix)),'g.','MarkerSize',20)

trix = randsample(lhit,nPlot);
plot3(squeeze(N_null(101:end,dim,trix)),squeeze(N_potent(101:end,2,trix)),squeeze(N_potent(101:end,1,trix)),'r')
plot3(squeeze(N_null(101,dim,trix)),squeeze(N_potent(101,2,trix)),squeeze(N_potent(101,1,trix)),'r.','MarkerSize',20)
plot3(squeeze(N_null(gc,dim,trix)),squeeze(N_potent(gc,2,trix)),squeeze(N_potent(gc,1,trix)),'g.','MarkerSize',20)

hold off
xlabel('Null 1')
ylabel('Potent 2')
zlabel('Potent 1')



nPlot = 20;
gc = round(numel(obj.time)/2);
figure; hold on
trix = randsample(rhit,nPlot);
plot3(squeeze(N_null(101:end,dim,trix)),squeeze(N_null(101:end,2,trix)),squeeze(N_potent(101:end,1,trix)),'b')
plot3(squeeze(N_null(101,dim,trix)),squeeze(N_potent(101,2,trix)),squeeze(N_potent(101,1,trix)),'b.','MarkerSize',20)
plot3(squeeze(N_null(gc,dim,trix)),squeeze(N_potent(gc,2,trix)),squeeze(N_potent(gc,1,trix)),'g.','MarkerSize',20)

trix = randsample(lhit,nPlot);
plot3(squeeze(N_null(101:end,dim,trix)),squeeze(N_null(101:end,2,trix)),squeeze(N_potent(101:end,1,trix)),'r')
plot3(squeeze(N_null(101,dim,trix)),squeeze(N_potent(101,2,trix)),squeeze(N_potent(101,1,trix)),'r.','MarkerSize',20)
plot3(squeeze(N_null(gc,dim,trix)),squeeze(N_potent(gc,2,trix)),squeeze(N_potent(gc,1,trix)),'g.','MarkerSize',20)

hold off
xlabel('Null 1')
ylabel('Null 2')
zlabel('Potent 1')










