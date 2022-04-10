clear,clc,close all

% finds kinematic modes for a single session and determines a null and
% potent space based on the linear equation V = NW, where M is a video feature matrix and
% N is a matrix of single trial neural activity. W is the transformation
% matrix between video features and neural activity. null(W) is a basis for
% the null space of W and null(W') is a basis for the row space of W (which
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
params.lfads_or_fa = 'fa'; % 'lfads' or 'fa'
params.lfads_run = ''; % 'run3' , leave empty to use most recent run
params.fcut_post_fa = 31; % if performing FA, cutoff freq to smooth rates and factors with a butterworth filter


% % TODO: parameters for time points to use (and a lag b/w neural activity
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


%% VIDEO FEATURES

% getKinematicsFromVideo() returns 2 variables
% - kin: feature matrix of size (time,features,trials). Features defined in
%         params.traj_features
% - featLeg: legend corresponding to features in kin (for 2nd dimension)
[kin,featLeg] = getKinematicsFromVideo(obj,params,dat.trials);

% JAW ANGLE
jawAngle = getJawAngle(obj.time, obj, dat.trials, params.alignEvent); % (time, trials)
featLeg{end+1} = 'jaw_angle';

% MOTION ENERGY
% me is a struct with fields
% - data:       cell array (ntrials,1), with motion energy for each time point
% - moveThresh: a scalar - below this value of motion energy, animal is
%               said to be stationary
% - use:        0 or 1, indicating whether to use me.data. 1 if
%               a motionEnergy*.mat file is found, 0 if file not found
me = loadMotionEnergy(meta);
if me.use
    featLeg{end+1} = 'motion_energy';
end
% To generate a motionEnergy*.mat file for a session, see https://github.com/economolab/videoAnalysisScripts/blob/main/motionEnergy.m


% dimensionality reduction
% the video features will be highly correlated, so we will perform PCA/FA
% on the matrix of features to reduce the features to a set of factors that
% best explain the movement captured by the video recordings


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
% The null space of W forms subspace that we will call the     	NULL SPACE   of N
% The row space of W forms the subspace that we will call the   POTENT SPACE of N

%% NULL AND POTENT SPACE OF W



%% dPCA







